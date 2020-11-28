package rbsseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.SamLocusIterator; 
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.stat.inference.BinomialTest;

public class ScorePseudouridinePositions {
	//Filtering settings
	private double errorRate = 0.001;
	private double splitThresh = 0.5;
	private float pval = (float)0.05;
	private int minBsCov = 10;
	private int minNbsCov = 10;
	private int minBsDel = 5;
	private double maxNbsFrac = 0.01;
	private double minBsFrac = 0.005;
	private int flankLength = 4;
	private int delDistance = 5;
	private int hpLength = 6; //Length to consider homopolymer
	
	//File settings
	private File biomartFile = null;
	private File repbaseFile = null;
	private File ucscFile = null;
	private File bisulfiteAlignment = null;
	private File nonBisulfiteAlignment = null;
	private File preParsedFile = null;
	private File referenceFile = null;
	private File outputPrefix = null;

	//Data dictionaries
	private HashMap<String,String> refSeq = new HashMap<String,String>();
	private HashMap<String,String[]> biomartDict = new HashMap<String,String[]>();
	private HashMap<String,ArrayList<Feature>> flankDict = new HashMap<String,ArrayList<Feature>>();
	private HashMap<String,ArrayList<Feature>> geneDict= new HashMap<String,ArrayList<Feature>>();
	
	public HashMap<String,String> revComp = new HashMap<String,String>();

	//Counters!
	private int minCovBsFilter = 0;
	private int minCovNbsFilter = 0;
	private int minDelBsFilter = 0;
	private int maxNbsFractionFilter = 0;
	private int minBsFractionFilter = 0;
	private int minObsOK = 0;
	private int totalPositions = 0;
	private int afterCollapse = 0;
	
	private int lowConfidenceFlag = 0;
	private int highBackgroundExactFlag = 0;
	private int highBackgroundProxFlag = 0;
	private int homopolymerFlag = 0;
	private int flankFlag = 0;
	private int qvalue = 0;
	
	private int passed = 0;
	private int annotated = 0;
	private int unannotated = 0;
	
	private int sharedPositions = 0;
	private int sharedUsablePositions = 0;
	
	//Data containers
	private ArrayList<PositionGroup> finalPositionGroups = new ArrayList<PositionGroup>();
	private ArrayList<PositionGroup> originalPositionGroups = new ArrayList<PositionGroup>();
	private HashSet<String> nbsPositions = new HashSet<String>();
	
	public static void main(String[] args) {
		new ScorePseudouridinePositions(args);

	}
	
	public ScorePseudouridinePositions(String[] args) {
		//Setup dictionary
		revComp.put("A", "T");
		revComp.put("T", "A");
		revComp.put("C", "G");
		revComp.put("G", "C");
		revComp.put("N", "N");
		
		//Parse command line arguments
		System.out.println("Parsing command line arguments... ");
		processArgs(args);
		
		//Read in biomart annotations
		if (biomartFile != null) {
			System.out.println("Loading biomart anntotations... ");
			readBiomart();
		}
		
		
		//Read in UCSC annotations
		System.out.println("Loading ucsc annotations... ");
		readGeneTable();
		
		//Read in reference sequence
		System.out.println("Loading reference sequence... ");
		readReferenceSequence();
		
		//Read in repbase annotations
		if (repbaseFile != null) {
			System.out.println("Reading repbase annotations... ");
			readRepBase();
		}
		
		
		//Filter mpileup
		if (preParsedFile != null) {
		    System.out.println("Reading in pre-parsed file");
			parseExisting();
		} else if (nonBisulfiteAlignment != null) {
			System.out.println("Parsing bisulfite and non-bisulfite alignments");
			parseAlignmentFiles();
		} else {
			System.out.println("Parsing bisulfite alignments");
			parseSingleFile();
		}
			
		BinomialTest bt = new BinomialTest();
		
		//Annotating deletions
		System.out.println("Annotating deletions... ");
		annotatePositions(bt);
		
		//Annotating deletions
		System.out.println("Scanning for nearby deletions... ");
		findNearbyNbsDeletions();
		
		//Identifying potential shared sites
		System.out.println("Scanning for shared originating Ts...");
		removedSharedSites();
		
		//Calculating qvalue
		System.out.println("Calculating qvalues... ");
		calculateQvalue();
		
		//Writing out results
		System.out.println("Writing results to file...");
		writeResults();
		
		writeStats();
		System.out.println("Finished!");
	}
	
	private void writeStats() {
		System.out.println("\n\n******* Threshold Filtering ******");
		System.out.println(String.format("%d positions processed",totalPositions));
		System.out.println(String.format("Fewer than %d deletions in bisulfite sample: %d (%.4f%%)",minBsDel,minDelBsFilter,(float)minDelBsFilter/totalPositions*100));
		System.out.println(String.format("Fewer than %d coverage in bisulfite sample: %d (%.4f%%)",minBsCov,minCovBsFilter,(float)minCovBsFilter/totalPositions*100));
		System.out.println(String.format("Deletion rate lower than %.4f in bisulfite sample: %d (%.4f%%)",minBsFrac,minBsFractionFilter,(float)minBsFractionFilter/totalPositions*100));
		System.out.println(String.format("Fewer than %d coverage in nbs sample: %d (%.4f%%)",minNbsCov,minCovNbsFilter,(float)minCovNbsFilter/totalPositions*100));
		System.out.println(String.format("Deletion rate higher than %.4f in nbs sample: %d (%.4f%%)",maxNbsFrac,maxNbsFractionFilter,(float)maxNbsFractionFilter/totalPositions*100));
		System.out.println(String.format("Passing positions: %d (%.4f%%)",minObsOK,(float)minObsOK/totalPositions*100));
		System.out.println(String.format("Positons after collapsing: %d (%.4f%%)",afterCollapse,(float)afterCollapse/minObsOK*100));
		System.out.println("\n\n******* Artifact and Confidence Filtering ******");
		System.out.println(String.format("Low confidence (pval %.4f): %d (%.4f%%)",pval,lowConfidenceFlag,(float)lowConfidenceFlag/afterCollapse*100));
		System.out.println(String.format("High Background Exact (pval %.4f): %d (%.4f%%)",pval,highBackgroundExactFlag,(float)highBackgroundExactFlag/afterCollapse*100));
		System.out.println(String.format("High Background Proximity (%d bp): %d (%.4f%%)",delDistance,highBackgroundProxFlag,(float)highBackgroundProxFlag/afterCollapse*100));
		System.out.println(String.format("Homopolymer (%d bp): %d (%.4f%%)",hpLength,homopolymerFlag,(float)homopolymerFlag/afterCollapse*100));
		System.out.println(String.format("Splice Junction (%d bp): %d (%.4f%%)",flankLength,flankFlag,(float)flankFlag/afterCollapse*100));
		System.out.println(String.format("Annotation Ready: %d (%.4f%%)",passed,(float)passed/afterCollapse*100));
	
		
		System.out.println("\n\n******* Annotation ******");
		System.out.println(String.format("Annotated: %d (%.4f%%)",annotated,(float)annotated/passed*100));
		System.out.println(String.format("Unannotated: %d (%.4f%%)",unannotated,(float)unannotated/passed*100));
		
		System.out.println("\n\n******* Annotation ******");
		System.out.println(String.format("Shared Originating T all: %d (%.4f%%)",sharedPositions,(float)sharedPositions/passed*100));
		System.out.println(String.format("Shared Originating T Usable: %d (%.4f%%)",sharedUsablePositions,(float)sharedUsablePositions/passed*100));
		
		System.out.println("\n\n******* Qvalue ******");
		System.out.println(String.format("Corrected positions (Low confidence + Annotation Ready - Shared Usable): %d",passed + lowConfidenceFlag - sharedUsablePositions));
		System.out.println(String.format("qvalue passed: %d (%.4f%%)",qvalue,(float)qvalue/(passed + lowConfidenceFlag - sharedUsablePositions)*100));
	
	}
	
	private void removedSharedSites() {
		HashMap<String,Double> maxScore = new HashMap<String,Double>();
		for (PositionGroup pg: finalPositionGroups) {
			String key = pg.getMaxPos().getChrom() + ":" + String.valueOf(pg.getMaxPos().getPredPos());
			Double fraction = pg.getMaxPos().getBsPer();
			if (maxScore.containsKey(key)) {
				if (maxScore.get(key) < fraction) {
					maxScore.put(key,fraction);
				}
			}  else {
				maxScore.put(key,fraction);
			}
		}
		
		ArrayList<PositionGroup> sharedRemoved = new ArrayList<PositionGroup>();
		for (PositionGroup pg: finalPositionGroups) {
			String key = pg.getMaxPos().getChrom() + ":" + String.valueOf(pg.getMaxPos().getPredPos());
			Double fraction = pg.getMaxPos().getBsPer();
			if (maxScore.get(key).equals(fraction)) {
				sharedRemoved.add(pg);
			} else {
				sharedPositions += 1;
				if (pg.getFilterFlag().equals("Annotated") || pg.getFilterFlag().equals("Intron/Intergenic") || pg.getFilterFlag().equals("LowConfidence")) {
					sharedUsablePositions += 1;
				} 
			}
		}
		
		finalPositionGroups = sharedRemoved;	
	}
	
	private void annotatePositions(BinomialTest bt) {
		ArrayList<Feature> sortedFlank = null;
		ArrayList<Feature> sortedGene = null;
		
		
		
		Pattern p1 = Pattern.compile("^[CT]*$");
		Pattern p2 = Pattern.compile("^[GA]*$");
		
		for (PositionGroup p: finalPositionGroups) {
			p.determineBaseToReport(refSeq);
			//p.calculateFisher(errorRate, fe);
			p.calculateBinomal(errorRate, bt);
			String filterFlag = "NA";
			String currentChrom = "NA";
			String geneName = "NA";
			String biotype = "NA";
			if (p.getBsPval() > pval) {
				lowConfidenceFlag++;
				
				filterFlag = "LowConfidence";
			} else if (p.getNbsPval() <= pval) {
				highBackgroundExactFlag++;
				filterFlag = "HighBackgroundExact";
			} 
		
			
			String chrom = p.getMaxPos().getChrom();
			int pos = p.getMaxPos().getPos();
			
			//Check to see if deletion is in flanking region.
			if (filterFlag.equals("NA")) {
				if (!chrom.equals(currentChrom) || sortedFlank == null) {
					if (flankDict.containsKey(chrom)) {
						sortedFlank = flankDict.get(chrom);
						Collections.sort(sortedFlank, new FeatureComparator());
					} else {
						sortedFlank = new ArrayList<Feature>();
					}
				}
				
				for (Feature f: sortedFlank) {
					if (f.getStart() > pos) {
						break;
					}
					if (pos > f.getStart() && pos < f.getEnd()) {
						filterFlag = "InFlank";
						flankFlag++;
						break;
					}
				}
			}
			
			
			//Add annotation
			if (filterFlag.equals("NA")) {
				int searchPos =  pos;
				
				String seq = refSeq.get(chrom).substring(searchPos,searchPos+hpLength).toUpperCase();
				Matcher m1 = p1.matcher(seq);
				Matcher m2 = p2.matcher(seq);
								
				if (m1.matches()) {
					filterFlag = "Homopolymer";
					homopolymerFlag++;
				} else if (m2.matches()) {
					filterFlag = "Homopolymer";
					homopolymerFlag++;
				} else {
					if (!chrom.equals(currentChrom) || sortedGene == null) {
						if (geneDict.containsKey(chrom)) {
							sortedGene = geneDict.get(chrom);
							Collections.sort(sortedGene,new FeatureComparator());
						} else {
							sortedGene = new ArrayList<Feature>();
						}
					}
					
					boolean found = false;
					ArrayList<String> biotypes = new ArrayList<String>();
					ArrayList<String> names = new ArrayList<String>();
					for (Feature gene: sortedGene) {
						if (gene.getStart() > pos) {
							break;
						} 
						
						if (pos > gene.getStart() && pos < gene.getEnd()) {
							if (!names.contains(gene.getName())) {
								names.add(gene.getName());
								biotypes.add(gene.getBiotype());
							}
							found = true;
						}
					}
					
					if (!found) {
						filterFlag = "Intron/Intergenic";
						unannotated++;
						passed++;
					} else {
						filterFlag = "Annotated";
						geneName = join(names, ";");
						biotype = join(biotypes, ";");
						annotated++;
						passed++;
					}
				}
			}
			
			p.addAnnotation(geneName,filterFlag,biotype);
			
			currentChrom = chrom;
		}
	}
	
	private String join(ArrayList<String> list, String sep) {
		StringBuilder sb = new StringBuilder("");
		for (String l: list) {
			sb.append(sep + l);
		}
		String finalString = sb.toString().substring(sep.length());
		return finalString;
	}
	
	private void calculateQvalue() {
		
		ArrayList<PositionGroup> sortedPosition = new ArrayList<PositionGroup>();
    	ArrayList<Double> pvalueList = new ArrayList<Double>();
		
		for (PositionGroup p: finalPositionGroups) {
			if (p.getFilterFlag().equals("Annotated") || p.getFilterFlag().equals("Intron/Intergenic") || p.getFilterFlag().equals("LowConfidence")) {
				sortedPosition.add(p);
			} 
		}
		
		Collections.sort(sortedPosition,new PositionComparator());
		
		for (PositionGroup p: sortedPosition) {
			pvalueList.add(p.getBsPval());
		}
		
		double[] pvalue = new double[pvalueList.size()];
		for (int i=0; i< pvalue.length ;i++) {
			pvalue[i] = pvalueList.get(i);
		}
		
		benjaminiHochbergCorrect(pvalue);
		
		for (int i=0; i < pvalue.length; i++) {
			if (pvalue[i] < 0.05) {
				sortedPosition.get(i).addQvalue(pvalue[i], false);
				qvalue++;
			} else {
				sortedPosition.get(i).addQvalue(pvalue[i], true);
			}
		}
		

	}
	
	private void findNearbyNbsDeletions() {
		for (PositionGroup p: finalPositionGroups) {
			boolean found = false;
			for (int i=p.getMaxPos().getPos()-delDistance;i<p.getMaxPos().getPos()+delDistance+1;i++) {
				String key = p.getMaxPos().getChrom() + ":" + String.valueOf(i);
				if (nbsPositions.contains(key)) {
					found = true;
				}
			}
			
			if (found) {
				if (p.getFilterFlag().equals("Annotated")) {
					annotated--;
					passed--;
					p.addAnnotation("NA", "HighBackgroundProx", "NA");
					highBackgroundProxFlag++;
				} else if (p.getFilterFlag().equals("Intron/Intergenic")) {
					unannotated--;
					passed--;
					p.addAnnotation("NA", "HighBackgroundProx", "NA");
					highBackgroundProxFlag++;
				}
			}
		}
	}
	
	private void writeResults() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputPrefix + ".results.txt"));
			
			if (finalPositionGroups.size() > 0) {
				bw.write(finalPositionGroups.get(0).getHeader());
			}
			
			for (PositionGroup p: finalPositionGroups) {
				bw.write(p.outputString(refSeq));
			}
			
			bw.close();
		} catch (IOException ioex) {
			System.out.println("Error writing to output file: " + ioex.getMessage());
			System.exit(1);
		}
	}
	
	private void parseSingleFile() {
		try {
			BufferedWriter bwStats = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".parsed.txt.gz"))));
		
			SamReader srBS = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bisulfiteAlignment);
			SamLocusIterator sliBS = new SamLocusIterator(srBS);
			sliBS.setEmitUncoveredLoci(false);
			
			PositionGroup currentGroup = new PositionGroup(splitThresh, hpLength);
			LocusInfo lBS = null;
			
			int counter = 0;
			while(true) {
				totalPositions++;
				//Get coverage information
				int countBS;
				int covBS;
				int countNBS = 0;
				int covNBS = 100;
				
				//strand
				int reverse;
				int forward;
				
				if (sliBS.hasNext()) {
					lBS = sliBS.next();
				} else {
					break;
				}
				
				//location
				String chrom;
				String position=null;
				countBS = lBS.getDeletionCount();
				covBS = lBS.getLocusCoverage();
				reverse = lBS.getReverseCount();
				forward = lBS.getFowardCount();
				
				chrom = lBS.getSequenceName();
				position = String.valueOf(lBS.getPosition());
				
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(counter + " " + chrom + " " + position);
				}
				counter += 1;
	
				currentGroup = processPosition(currentGroup, covNBS, covBS, countBS, countNBS, forward, reverse, chrom, position, bwStats);	
			}
			
			if (currentGroup.getPosList().size() > 0) {
				originalPositionGroups.add(currentGroup);
			}
			
			for (PositionGroup p: originalPositionGroups) {
				ArrayList<PositionGroup> split = p.splitGroup();
				finalPositionGroups.addAll(split);
			}
			
			afterCollapse = finalPositionGroups.size();
			bwStats.close();
			sliBS.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
	}
	
	/***************
	 * This function processes a bisulfite and non-bisulfite alignment file
	 */
	private void parseAlignmentFiles() {
		try {
			BufferedWriter bwStats = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".parsed.txt.gz"))));
			
			SamReader srBS = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bisulfiteAlignment);
			SamReader srNBS = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(nonBisulfiteAlignment);
			SamLocusIterator sliBS = new SamLocusIterator(srBS);
			sliBS.setEmitUncoveredLoci(false);
			SamLocusIterator sliNBS = new SamLocusIterator(srNBS);
			sliNBS.setEmitUncoveredLoci(false);
			String refBS = null;
			String refNBS = null;
			String refLast = null;
			
			PositionGroup currentGroup = new PositionGroup(splitThresh, hpLength);
			
			LocusInfo lBS = null;
			LocusInfo lNBS = null;
			if (sliBS.hasNext()) {
				lBS = sliBS.next();
				refBS = lBS.getSequenceName();
			}
			if (sliNBS.hasNext()) {
				lNBS = sliNBS.next();
				refNBS = lNBS.getSequenceName();
			}
			
			int counter = 0;
			while(true) {
				totalPositions++;
				//Get coverage information
				int countBS;
				int covBS;
				int countNBS;
				int covNBS;
				
				//strand
				int reverse;
				int forward;
				
				//location
				String chrom;
				String position=null;
				
				if (lBS == null && lNBS == null) {
					break;
				} else if (lBS != null && (lNBS == null || (lBS.getPosition() < lNBS.getPosition() && refNBS == refBS) || refNBS != refLast)) {
					countBS = lBS.getDeletionCount();
					covBS = lBS.getLocusCoverage();
					countNBS = 0;
					covNBS = 0;
					
					reverse = lBS.getReverseCount();
					forward = lBS.getFowardCount();
					
					chrom = lBS.getSequenceName();
					position = String.valueOf(lBS.getPosition());
					
					if (sliBS.hasNext()) {
						lBS = sliBS.next();
						refBS = lBS.getSequenceName();
						if (refBS.equals(refNBS)) {
							refLast = chrom;
						}
					} else {
						lBS = null;
					}
				} else  if (lNBS != null  && (lBS == null || (lNBS.getPosition() < lBS.getPosition() && refNBS == refBS) || refBS != refLast)) {
					countNBS = lNBS.getDeletionCount();
					covNBS = lNBS.getLocusCoverage();
					countBS = 0;
					covBS = 0;
					
					reverse = lNBS.getReverseCount();
					forward = lNBS.getFowardCount();
					
					chrom = lNBS.getSequenceName();
					position = String.valueOf(lNBS.getPosition());
					
					if (sliNBS.hasNext()) {
						lNBS = sliNBS.next();
						refNBS = lNBS.getSequenceName();
						if (refBS.equals(refNBS)) {
							refLast = chrom;
						}
					} else {
						lNBS = null;
					}
				} else {
					countNBS = lNBS.getDeletionCount();
					covNBS = lNBS.getLocusCoverage();
					countBS = lBS.getDeletionCount();
					covBS = lBS.getLocusCoverage();
					
					reverse = lNBS.getReverseCount() + lBS.getReverseCount();
					forward = lNBS.getFowardCount() + lBS.getFowardCount();
					
					chrom = lNBS.getSequenceName();
					position = String.valueOf(lNBS.getPosition());
					
					if (sliNBS.hasNext()) {
						lNBS = sliNBS.next();
						refNBS = lNBS.getSequenceName();
					} else {
						lNBS = null;
					}
					
					if (sliBS.hasNext()) {
						lBS = sliBS.next();
						refBS = lBS.getSequenceName();
					} else {
						lBS = null;
					}
					
					if (refBS.equals(refNBS)) {
						refLast = chrom;
					}
				}
				
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(counter + " " + chrom + " " + position);
				}
				counter += 1;
				
				currentGroup = processPosition(currentGroup, covNBS, covBS, countBS, countNBS, forward, reverse, chrom, position, bwStats);

			}
						
			if (currentGroup.getPosList().size() > 0) {
				originalPositionGroups.add(currentGroup);
			}
			
			
			for (PositionGroup p: originalPositionGroups) {
				ArrayList<PositionGroup> split = p.splitGroup();
				finalPositionGroups.addAll(split);
			}
			originalPositionGroups.clear();
			
			afterCollapse = finalPositionGroups.size();
			bwStats.close();
			sliBS.close();
			sliNBS.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
	} 
			
	private void parseExisting() {
		try {
			
			BufferedWriter bwStats = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".parsed.txt.gz"))));
			PositionGroup currentGroup = new PositionGroup(splitThresh, hpLength);
			BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(preParsedFile))));
			String temp = null;
			
			int counter = 0;
			while((temp = br.readLine()) != null) {
				String[] parts = temp.split("\t");
				
				totalPositions++;
				
				//Get coverage information
				int countBS = Integer.parseInt(parts[2]);
				int covBS = Integer.parseInt(parts[3]);
				int countNBS = Integer.parseInt(parts[4]);
				int covNBS = Integer.parseInt(parts[5]);
				
				//strand
				int forward = Integer.parseInt(parts[6]);
				int reverse = Integer.parseInt(parts[7]);
				
				//location
				String chrom = parts[0];
				String position = parts[1];
				
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(counter + " " + chrom + " " + position);
				}
				counter++;
						
				currentGroup = processPosition(currentGroup, covNBS, covBS, countBS, countNBS, forward, reverse, chrom, position, bwStats);
			}
			
			if (currentGroup.getPosList().size() > 0) {
				originalPositionGroups.add(currentGroup);
			}
			
			
			for (PositionGroup p: originalPositionGroups) {
				ArrayList<PositionGroup> split = p.splitGroup();
				finalPositionGroups.addAll(split);
			}
			
			afterCollapse = finalPositionGroups.size();
			br.close();
			bwStats.close();
		} catch (IOException ioex) {
			System.out.println("Error reading file");
		}
	}
	
	/************************
	 * This method takes basic information about a position and does basic filtering.  If this position passes, it is either added to an existing group
	 * or used to create a new group.
	 * 
	 * @param currentGroup
	 * @param covNBS
	 * @param covBS
	 * @param countBS
	 * @param countNBS
	 * @param forward
	 * @param reverse
	 * @param chrom
	 * @param position
	 * @return PositonGroup
	 */
	private PositionGroup processPosition(PositionGroup currentGroup, int covNBS, int covBS, int countBS, int countNBS, int forward, int reverse, String chrom, String position, BufferedWriter bwParsed) throws IOException{
		
		double fracNBS = 0;
		double fracBS = 0;
		if (covNBS > 0) {
			fracNBS = (double)countNBS / covNBS;
		}
		
		if (covBS >= 0) {
			fracBS = (double)countBS / covBS;
		}
		
		String status = "";
		if (countBS < minBsDel) {
			minDelBsFilter++;
			status += "minDelBs;";
		}
		
		if (fracBS < minBsFrac) {
			minBsFractionFilter++;
			status += "minBsFrac;";
		}
		
		if (covBS < minBsCov) {
			minCovBsFilter++;
			status += "minBsCov;";
		}
		
		if (fracNBS > maxNbsFrac) {
			maxNbsFractionFilter++;
			status += "minNbsFrac;";
			//Store potential NBS deletions
			if (covNBS >= minNbsCov && countNBS >= minBsDel) {
				nbsPositions.add(chrom + ":" + position);
			}
		} 
		
		if (covNBS < minNbsCov) {
			minCovNbsFilter++;
			status += "minNbsCov;";
		}
		
		if (!status.equals("")) {
			bwParsed.write(chrom + "\t" + position + "\t" + countBS + "\t" + covBS + "\t" + countNBS + "\t" + covNBS + "\t" + forward + "\t" + reverse + "\t" + status + "\n");
			return currentGroup;
		}
		
		bwParsed.write(chrom + "\t" + position + "\t" + countBS + "\t" + covBS + "\t" + countNBS + "\t" + covNBS + "\t" + forward + "\t" + reverse + "\tPASSED\n");
		minObsOK++;

		double perF = (double)forward / (forward + reverse);
		double perR = (double)reverse / (forward + reverse);
		
		String strand = "NA";
		if (perF > perR) {
			if (perF < perR * 2) {
				strand = "+*";
			} else {
				strand = "+";
			}
		} else {
			if (perR < perF * 2) {
				strand = "-*";
			} else {
				strand = "-";
			}
		}
		
		Integer pos = Integer.parseInt(position)-1;
		String base = refSeq.get(chrom).substring(pos, pos+1);
		Position p = new Position(chrom,pos,strand,base,covBS,countBS,covNBS,countNBS);
		
		//Collapse
		if ((currentGroup.getLastPos() + 1) == p.getPos()) {
			currentGroup.addPosition(p);
		} else {
			if (currentGroup.getPosList().size() > 0) {
				originalPositionGroups.add(currentGroup);
			}
			currentGroup = new PositionGroup(splitThresh, hpLength);
			currentGroup.addPosition(p);
		}
		return currentGroup;
	}
	
	/* 
	 * Process reference files
	 */
					
	private void readRepBase() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(repbaseFile));
			
			String temp = "";
			while((temp = br.readLine()) != null) {
				String[] parts = temp.split("\t");
				Feature repbaseFeature = new Feature(Integer.parseInt(parts[1]),Integer.parseInt(parts[2]),parts[3],"repbase");
				if (!geneDict.containsKey(parts[0])) {
					geneDict.put(parts[0], new ArrayList<Feature>());
				}
				geneDict.get(parts[0]).add(repbaseFeature);
			}
			br.close();
		} catch (IOException ioex) {
			System.out.println("Error reading repbase file: " + ioex.getMessage());
		}
	}
	
	
	private void readGeneTable() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(ucscFile));
			
			String temp = "";
			while((temp = br.readLine()) != null) {
				String[] parts = temp.split("\t");
				String chrom = parts[2];
				if (!flankDict.containsKey(chrom)) {
					flankDict.put(chrom,new ArrayList<Feature>());
				}
				if (!geneDict.containsKey(chrom)) {
					geneDict.put(chrom,new ArrayList<Feature>());
				}
				
				if (parts[9].endsWith(",")) {
					parts[9] = parts[9].substring(0,parts[9].length());
				}
				if (parts[10].endsWith(",")) {
					parts[10] = parts[10].substring(0,parts[10].length());
 				}
				
				String[] starts = parts[9].split(",");
				String[] ends =  parts[10].split(",");
				
				for (int i=0; i<starts.length; i++) {
					String geneName = "NA";
					String biotype = "NA";
					
					if (biomartDict.containsKey(parts[0])) {
						String[] annotations = biomartDict.get(parts[0]);
						if (annotations[0].equals("")) {
							geneName = parts[0];
						} else {
							geneName = annotations[0];
						}
						
						biotype = annotations[1];
					}
					
					int start = Integer.parseInt(starts[i]);
					int end = Integer.parseInt(ends[i]);
					
					
					Feature geneFeature = new Feature(start,end,geneName,biotype);
					geneDict.get(chrom).add(geneFeature);
					
					if (i != 0) {
						Feature upFlank = new Feature(start-flankLength,start+flankLength,geneName,"NA");
						flankDict.get(chrom).add(upFlank);
					}
					if (i != starts.length-1) {
						Feature downFlank = new Feature(end-flankLength,end+flankLength,geneName,"NA");
						flankDict.get(chrom).add(downFlank);
					}
				}
			}
		
			br.close();
		} catch (IOException ioex) {
			System.out.println("Error reading genetable, exiting: " + ioex.getMessage());
		}
	}

	private void readBiomart() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(biomartFile));
			
			String temp = "";
			while ((temp = br.readLine()) != null) {
				String[] parts = temp.split("\t");
				String[] annotations = new String[2];
				annotations[0] = parts[2];
				annotations[1] = parts[3];
				biomartDict.put(parts[0], annotations);
			}
			
			br.close();
		} catch (IOException ioex) {
			System.out.println("Error reading biomart file, exiting: " + ioex.getMessage());
		}
	}
	
	private void readReferenceSequence() {
		
		
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(referenceFile));
			String temp = null;
			Pattern p = Pattern.compile(">(\\S+).*$");
			String chrom = null;
			StringBuilder seq = null;
			while ((temp = br.readLine()) != null) {
				Matcher m = p.matcher(temp);
				if (m.matches()) {
					if (chrom != null) {
						refSeq.put(chrom, seq.toString());
					}
					chrom = m.group(1);
					seq = new StringBuilder("");
				} else {
					seq.append(temp);
				}
			}
			refSeq.put(chrom,seq.toString());
			br.close();
		} catch (IOException ioex) {
			System.out.println("Error reading reference fasta file: " + ioex.getMessage());
		}
	}
	
	
	private void processArgs(String[] args) {
		CommandLineParser parser = new DefaultParser();
		
		//Create options
		Options options = new Options();
		options.addOption(Option.builder("a").longOpt("bis").desc("Bisulfite alignment file in bam format.").type(File.class).hasArg().build());
		options.addOption(Option.builder("b").longOpt("non-bis").desc("Non-bisulfite alignment file in bam format.").type(File.class).hasArg().build());
		options.addOption(Option.builder("c").longOpt("pre-parsed").desc("Preparsed file generated from this program.").type(File.class).hasArg().build());
		options.addOption(Option.builder("d").longOpt("out-prefix").desc("Output file prefix.").type(File.class).hasArg().required().build());
		
		options.addOption(Option.builder("e").longOpt("ann-file").desc("Gene models in ucsc refflat format.").type(File.class).hasArg().required(true).build());
		options.addOption(Option.builder("f").longOpt("ref-file").desc("Reference file in fasta format.").type(File.class).hasArg().required(true).build());
		options.addOption(Option.builder("g").longOpt("repbase-file").desc("Repbaes entries in bed format.").type(File.class).hasArg().required(true).build());
		options.addOption(Option.builder("h").longOpt("biomart-file").desc("Biomart annotations.").type(File.class).hasArg().required(true).build());
		
		options.addOption(Option.builder("i").longOpt("min-del-obs").desc("Minimum deletions in bisulfite sample. Default 5.").type(Number.class).hasArg().build());
		options.addOption(Option.builder("j").longOpt("min-bis-cov").desc("Minimum coverage in bisulfite sample. Default 5.").type(Number.class).hasArg().build());
		options.addOption(Option.builder("k").longOpt("min-bs-frac").desc("Minimum fraction deletion in bisulfite sample. Default 0.005.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("l").longOpt("min-nbs-cov").desc("Minimum coverage in non-bisulfite sample. Default 10.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("m").longOpt("max-nbs-frac").desc("Max fraction deltion in non-bisulfite sample. Default 0.01.").type(Number.class).hasArg().build());	
		
		options.addOption(Option.builder("n").longOpt("homopolymer").desc("Deletions found at the end of homopolymers of length equal or greater to this number will be flagged. Default 6.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("o").longOpt("error-rate").desc("Error rate used in binomial test.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("p").longOpt("p-value").desc("p-value threshold. Collapsed positions greater than this threshold will not be annontated.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("q").longOpt("split-thresh").desc("Positions aren't collapsed if the fraction deletion is less than q * max fraction deletion of group. Default 0.5.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("r").longOpt("del-dist").desc("Positions within r bp of a NBS deletion are filtered out.").type(Number.class).hasArg().build());	
		options.addOption(Option.builder("s").longOpt("flank-dist").desc("Positions within s bp a exon boundary are filtered out.").type(Number.class).hasArg().build());	
		
		options.addOption("x","help",false,"Print help message and exit");
		
		try {
			if (args.length == 0) {
				printHelp(options);
				System.exit(0);
			}
			
			CommandLine line = parser.parse(options, args);
						
			if (line.hasOption("help")) {
				printHelp(options);
				System.exit(0);
			}
			
			//Parse input / output files
			if (line.hasOption("bis")) {
				bisulfiteAlignment = (File)line.getParsedOptionValue("bis");
				if (!bisulfiteAlignment.exists()) {
					System.out.println("Bisulfite alignment file specfied does not exist, exiting.");
					System.exit(1);
				}
			} 
			
			if (line.hasOption("non-bis")) {
				nonBisulfiteAlignment = (File)line.getParsedOptionValue("non-bis");
				if (!nonBisulfiteAlignment.exists()) {
					System.out.println("Non-bisulfite alignment file specfied does not exist, exiting.");
					System.exit(1);
				}
			} 
			
			if (line.hasOption("pre-parsed")) {
				preParsedFile = (File)line.getParsedOptionValue("pre-parsed");
				if (!preParsedFile.exists()) {
					System.out.println("Pre-parsed file specfied does not exist, exiting.");
					System.exit(1);
				}
			} 
			
			if (bisulfiteAlignment == null && preParsedFile == null) {
				System.out.println("Neither the bisuflite alignment file or alignment file were specified, exiting");
				System.exit(1);
			}
			
			//Parse references
			outputPrefix = (File)line.getParsedOptionValue("out-prefix");
			referenceFile = (File)line.getParsedOptionValue("ref-file");
			biomartFile = (File)line.getParsedOptionValue("biomart-file");
			ucscFile = (File)line.getParsedOptionValue("ann-file");
			repbaseFile = (File)line.getParsedOptionValue("repbase-file");
			
			//Parse optional settings
			if (line.hasOption("min-del-obs")) {
				minBsDel = (Integer)line.getParsedOptionValue("min-del-obs");
			}
			if (line.hasOption("min-bis-cov")) {
				minBsCov = (Integer)line.getParsedOptionValue("min-bis-cov");
			}
			if (line.hasOption("min-bs-frac")) {
				minBsFrac = (Double)line.getParsedOptionValue("min-bs-frac");
			}
			if (line.hasOption("min-nbs-cov")) {
				minNbsCov = (Integer)line.getParsedOptionValue("min-nbs-cov");
			}
			if (line.hasOption("max-nbs-frac")) {
				maxNbsFrac = (Double)line.getParsedOptionValue("max-nbs-frac");
			}
			if (line.hasOption("homopolymer")) {
				hpLength = (Integer)line.getParsedOptionValue("homopolymer");
			}
			if (line.hasOption("error-rate")) {
				errorRate = (Double)line.getParsedOptionValue("error-rate");
			}
			if (line.hasOption("p-value")) {
				pval = (Float)line.getParsedOptionValue("p-value");
			}
			if (line.hasOption("split-thresh")) {
				splitThresh = (Integer)line.getParsedOptionValue("split-thresh");
			}
			
			
				
		} catch (ParseException exp) {
			System.out.println("Error parsing command line arguments: " + exp.getMessage());
			printHelp(options);
			System.exit(1);
		}
	
	}
	
	private void printHelp(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.setWidth(500);
	    formatter.printHelp("This application scans alignments or pre-processed files for deletions.  Deletions that pass basic thresholds are collapsed and tested for significance "
	    		+ "using a binomial test.  Collapsed positions are also checked to make sure they don't meet criteria for potenial false positives, such as proximity to exons, homopolymers "
	    		+ "or deletions in nbs regions.\n\n", options );
	}
	
	/**Code borrowed from USeq**/
	/**Assumes sorted pvalues, descending. Alters the input array.*/
	public static void benjaminiHochbergCorrect(double[] sortedPValDecending){
		double num = sortedPValDecending.length;
		double prior = 1;
		for (int i=1; i< sortedPValDecending.length; i++){
			sortedPValDecending[i] = sortedPValDecending[i] * num / (num-i);
			if(sortedPValDecending[i] < prior) prior = sortedPValDecending[i]; 
			else sortedPValDecending[i] = prior;
		}
	}

	
	
	
	
	
}



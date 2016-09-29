package rbsseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.SamLocusIterator; 
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CreateMethTable {

	//File settings
	private File alignmentFile = null;
	private File forwardOut = null;
	private File reverseOut = null;
	private File referenceFile = null;
	private String mode = "nope";
	
	private HashMap<String,String> refSeq = new HashMap<String,String>();
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CreateMethTable(args);

	}
	
	public CreateMethTable(String[] args) {
		this.processArgs(args);
		if (mode.equals("meth")) {
			methTable(args);
		} else if (mode.equals("allCG")) {
			allCG(args);
		} else if (mode.equals("allAT")) {
			allAT(args);
		} else if (mode.equals("highAT")) {
			highAT(args);
		} else {
			System.out.println("Don't recognize mode: " + this.mode + ", exiting!");
			System.exit(1);
		}
	}
	
	private void methTable(String args[]) {
		
		this.readReferenceSequence();
		try {
			BufferedWriter bwF = new BufferedWriter(new FileWriter(this.forwardOut));
			BufferedWriter bwR = new BufferedWriter(new FileWriter(this.reverseOut));
			String header = "Chrom\tCoord\tStrand\tNuc\tDepth\t#A\t#C\t#G\t#T\t#N\t#-\t%Methylated\n";
			bwF.write(header);
			bwR.write(header);
			SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.alignmentFile);
			SamLocusIterator sli = new SamLocusIterator(sr);
			int counter = 0;
			for (LocusInfo li: sli) {
				String chrom = li.getSequenceName();
				int pos = li.getPosition();
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(chrom + " " + pos );
				}
				counter += 1;
				int[] forward = li.checkForward();
				int[] reverse = li.checkReverse();
				String base = refSeq.get(chrom).substring(pos-1, pos);
				if (forward[1] > 0) {
					int del = li.getDeletionCount();
					int cov = forward[0] + forward[1] + forward[2] + forward[3] + forward[4];
				   
			    	double methylation = (double)forward[1] / (cov);
					String result = String.format("%s\t%d\tF\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, forward[0], 
							forward[1], forward[2], forward[3], forward[4], del, methylation);
					bwF.write(result);
				    
					
				} else if (reverse[2] > 0) {
					int del = li.getDeletionCount();
					int cov = reverse[0] + reverse[1] + reverse[2] + reverse[3] + reverse[4];
					double methylation = (double)reverse[2] / (cov);
					String result = String.format("%s\t%d\tR\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, reverse[0], 
							reverse[1], reverse[2], reverse[3], reverse[4], del, methylation);
					bwR.write(result);
					
					
				}
			}
			bwF.close();
			bwR.close();
			sli.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
		
	}
	
	private void allCG(String args[]) {
		this.readReferenceSequence();
		try {
			BufferedWriter bwF = new BufferedWriter(new FileWriter(this.forwardOut));
			BufferedWriter bwR = new BufferedWriter(new FileWriter(this.reverseOut));
			String header = "Chrom\tCoord\tStrand\tNuc\tDepth\t#A\t#C\t#G\t#T\t#N\t#-\t%Methylated\n";
			bwF.write(header);
			bwR.write(header);
			SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.alignmentFile);
			SamLocusIterator sli = new SamLocusIterator(sr);
			int counter = 0;
			for (LocusInfo li: sli) {
				String chrom = li.getSequenceName();
				int pos = li.getPosition();
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(chrom + " " + pos );
				}
				counter += 1;
				int[] forward = li.checkForward();
				int[] reverse = li.checkReverse();
				String base = refSeq.get(chrom).substring(pos-1, pos);
				if (base.equals("C")) {
					int del = li.getDeletionCount();
					int cov = forward[0] + forward[1] + forward[2] + forward[3] + forward[4];
				    
				    double methylation = (double)forward[1] / (cov);
					String result = String.format("%s\t%d\tF\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, forward[0], 
							forward[1], forward[2], forward[3], forward[4], del, methylation);
					bwF.write(result);
				    
					
				} else if (base.equals("G")) { 
					
					int del = li.getDeletionCount();
					int cov = reverse[0] + reverse[1] + reverse[2] + reverse[3] + reverse[4];
					
					double methylation = (double)reverse[2] / (cov);
					String result = String.format("%s\t%d\tR\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, reverse[0], 
							reverse[1], reverse[2], reverse[3], reverse[4], del, methylation);
					bwR.write(result);
					
					
				}
			}
			bwF.close();
			bwR.close();
			sli.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
	}
	
	private void highAT(String args[]) {
		this.readReferenceSequence();
		try {
			BufferedWriter bwF = new BufferedWriter(new FileWriter(this.forwardOut));
			BufferedWriter bwR = new BufferedWriter(new FileWriter(this.reverseOut));
			String header = "Chrom\tCoord\tStrand\tNuc\tDepth\t#A\t#C\t#G\t#T\t#N\t#-\t%Methylated\n";
			bwF.write(header);
			bwR.write(header);
			SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.alignmentFile);
			SamLocusIterator sli = new SamLocusIterator(sr);
			int counter = 0;
			for (LocusInfo li: sli) {
				String chrom = li.getSequenceName();
				int pos = li.getPosition();
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(chrom + " " + pos );
				}
				counter += 1;
				int[] forward = li.checkForward();
				int[] reverse = li.checkReverse();
				String base = refSeq.get(chrom).substring(pos-1, pos);
				if (base.equals("A")) {
					int del = li.getDeletionCount();
					int cov = forward[0] + forward[1] + forward[2] + forward[3] + forward[4];
				    if (cov >= 100){
				    	double methylation = (double)forward[1] / (cov);
						String result = String.format("%s\t%d\tF\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, forward[0], 
								forward[1], forward[2], forward[3], forward[4], del, methylation);
						bwF.write(result);
				    }
					
				} else if (base.equals("T")) { 
					int del = li.getDeletionCount();
					int cov = reverse[0] + reverse[1] + reverse[2] + reverse[3] + reverse[4];
					if (cov >= 100) {
						double methylation = (double)reverse[2] / (cov);
						String result = String.format("%s\t%d\tR\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, reverse[0], 
								reverse[1], reverse[2], reverse[3], reverse[4], del, methylation);
						bwR.write(result);
					}
					
				}
			}
			bwF.close();
			bwR.close();
			sli.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
	}
	
	private void allAT(String args[]) {
		this.readReferenceSequence();
		try {
			BufferedWriter bwF = new BufferedWriter(new FileWriter(this.forwardOut));
			BufferedWriter bwR = new BufferedWriter(new FileWriter(this.reverseOut));
			String header = "Chrom\tCoord\tStrand\tNuc\tDepth\t#A\t#C\t#G\t#T\t#N\t#-\t%Methylated\n";
			bwF.write(header);
			bwR.write(header);
			SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.alignmentFile);
			SamLocusIterator sli = new SamLocusIterator(sr);
			int counter = 0;
			for (LocusInfo li: sli) {
				String chrom = li.getSequenceName();
				int pos = li.getPosition();
				if (counter % 5000000 == 0 && counter != 0) {
					System.out.println(chrom + " " + pos );
				}
				counter += 1;
				int[] forward = li.checkForward();
				int[] reverse = li.checkReverse();
				String base = refSeq.get(chrom).substring(pos-1, pos);
				if (base.equals("A")) {
					int del = li.getDeletionCount();
					int cov = forward[0] + forward[1] + forward[2] + forward[3] + forward[4];
			    	double methylation = (double)forward[1] / (cov);
					String result = String.format("%s\t%d\tF\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, forward[0], 
							forward[1], forward[2], forward[3], forward[4], del, methylation);
					bwF.write(result);
				    
					
				} else if (base.equals("T")) { 
					int del = li.getDeletionCount();
					int cov = reverse[0] + reverse[1] + reverse[2] + reverse[3] + reverse[4];
					double methylation = (double)reverse[2] / (cov);
					String result = String.format("%s\t%d\tR\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", chrom, pos, base, cov, reverse[0], 
							reverse[1], reverse[2], reverse[3], reverse[4], del, methylation);
					bwR.write(result);
				}
			}
			bwF.close();
			bwR.close();
			sli.close();
		} catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
	}
	
	
	
	private void readReferenceSequence() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.referenceFile));
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
					seq.append(temp.toUpperCase());
				}
			}
			refSeq.put(chrom,seq.toString());
			br.close();
		} catch (IOException ioex) {
			System.out.println("Error reading reference fasta file: " + ioex.getMessage());
		}
	}
	

	

	
	
	
	private void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': alignmentFile = new File(args[++i]); break;
					case 'g': referenceFile = new File(args[++i]); break;
					case 'f': forwardOut = new File(args[++i]); break;
					case 'r': reverseOut = new File(args[++i]); break;
					case 'm': mode = args[++i]; break;
					default: printErrorAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					printErrorAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (this.alignmentFile == null) {
			System.out.println("Alignment file not specified file not specfied: -a");
			System.exit(1);
		}
		if (!this.alignmentFile.exists()) {
			System.out.println("Specified alignment file does not exist: " + this.alignmentFile.getAbsolutePath());
			System.exit(1);
		}
		
		if (this.referenceFile == null) {
			System.out.println("Reference file not specified file not specfied: -g");
			System.exit(1);
		}
		if (!this.referenceFile.exists()) {
			System.out.println("Specified reference file does not exist: " + this.referenceFile.getAbsolutePath());
			System.exit(1);
		}
		
		if (this.forwardOut == null) {
			System.out.println("Forward output file not specified: -f");
			System.exit(1);
		}
		if (this.reverseOut == null) {
			System.out.println("Reverse output file not specified: -r");
			System.exit(1);
		}
		

	
	}
	
	public static void printDocs() { 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                       Make Meth Table                           **\n" +
				"**************************************************************************************\n" +
				"This replaces the python script that parses pairwise alignments\n" +
				

				"\nRequired Options:\n"+
				"-a alignment file. Path to bisulfite alignment file.\n" +
				"   AND \n" +
				"-f full path to the forward strand output file\n" + 
				"-r full path to the reverse strand output file\n" +
				"-g full path to the reference file\n" +
				"-m mode: \n" +
				"         meth = any observed C on forward or G on reverse.\n" +
				"         allCG = all C positions on forward or G on reverse.\n" +
				"         highAT = all A positions on forward or G on reverse, coverage 100x \n" +
				"         allAT = all A positions on forward or G on reverse.\n" +

				"\n"+

		        "**************************************************************************************\n");
	}
	
	private void printErrorAndExit(String message) {
		System.out.println(message);
		System.exit(1);
	}
}
	
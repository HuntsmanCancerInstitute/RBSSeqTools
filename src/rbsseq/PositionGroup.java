package rbsseq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;


public class PositionGroup {
	private ArrayList<Position> posList;
	private Position maxPos;
	private int lastPos;
	
	//Group stats
	private double bsPval;
	private double nbsPval;
	private double qvalue;
	
	//Flags
	private String filterFlag = "NA";
	private String qvalueFlag = "NA";
	
	//Annotations
	private String geneName = "NA";
	private String biotype = "NA";
	
	//Pattern matchers
	public final Pattern pF = Pattern.compile("^[CT]+");
	public final Pattern pR = Pattern.compile("^[GA]+");
	
	//Settings
	double splitThresh;
	int hpLength;
	
	public PositionGroup(double splitThresh, int hpLength) {
		posList = new ArrayList<Position>();
		this.splitThresh = splitThresh;
		this.hpLength = hpLength;
	}
	
	public void addPosition(Position pos) {
		posList.add(pos);
		double maxVal = 0;
		
		for (Position p: posList) {
			if (p.getBsPer() > maxVal) {
				maxPos = p;
				maxVal = p.getBsPer();
			}
		}
		lastPos = pos.getPos();
	}
	
	public void addAnnotation(String geneName, String filterFlag, String biotype) {
		this.geneName = geneName;
		this.filterFlag = filterFlag;
		this.biotype = biotype;
	}
	
	public String outputString(HashMap<String,String> refSeq) {
		/********** 
		 * This method finalizes the position information and returns a formatted string
		 **********/
		
		 //Generate position list
		 String allPosString = "";
		 
		 ArrayList<Integer> posLocList = new ArrayList<Integer>();
		 for (Position p: posList) {
			 posLocList.add(p.getPos());
		 }
		 
		 Collections.sort(posLocList);
		 for (Integer ap: posLocList) {
			 allPosString += ";" + String.valueOf(ap+1);
		 }
		 allPosString = allPosString.substring(1);
		 
		 String outputString = String.format("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%s\t%s\t%s\t%s\t%s\n",
				 maxPos.getChrom(),
				 maxPos.getPos()+1,
				 maxPos.getPredPos()+1,
				 allPosString,
				 posList.size(),
				 maxPos.getDirection(),
				 maxPos.getBase(),
				 maxPos.getBsDepth(),
				 maxPos.getBsDel(),
				 maxPos.getBsPer(),
				 bsPval,
				 qvalue,
				 maxPos.getNbsDepth(),
				 maxPos.getNbsDel(),
				 maxPos.getNbsPer(),
				 nbsPval,
				 filterFlag,
				 qvalueFlag,
				 maxPos.getBaseFlag(),
				 geneName,
				 biotype
				 );
		 
		 return outputString;	
	}
	
	public  String getHeader() {
		String header = "Chrom\tPositionCall\tPositionT\tPositionAll\tPositionCount\tStrand\tBase\tBSDepth\tBSDel\tBSFraction\tBSPval"
				+ "\tBSQval\tNBSDepth\tNBSDel\tNBSFraction\tNBSPval\tFilterFlag\tQvalueFlag\tBaseFlag\tGene\tBiotype\n";
		return header;
	}
	
	
	public void addQvalue(double qvalue, boolean flag) {
		this.qvalue = qvalue;
		if (flag) {
			this.qvalueFlag = "qvalue";
		} else {
			this.qvalueFlag = "OK";
		}
	}
	
	public void calculateBinomal(double errorRate, BinomialTest bt) {
		bsPval = bt.binomialTest(maxPos.getBsDepth(), maxPos.getBsDel(), errorRate, AlternativeHypothesis.GREATER_THAN);
		nbsPval = bt.binomialTest(maxPos.getNbsDepth(), maxPos.getNbsDel(), errorRate, AlternativeHypothesis.GREATER_THAN);
	}
	
	public ArrayList<PositionGroup> splitGroup() {
		//Prune!
		ArrayList<PositionGroup> newGroups = new ArrayList<PositionGroup>();
		double thresh = maxPos.getBsPer() * splitThresh;
		PositionGroup newGroup = new PositionGroup(splitThresh, hpLength);
		

		for(Position p: posList) {
			if (p.getBsPer() < thresh) {
				if (newGroup.posList.size() > 0) {
					newGroups.add(newGroup);
				}
				newGroup = new PositionGroup(splitThresh, hpLength);
			} else {
				newGroup.addPosition(p);
			}
		}
		if (newGroup.posList.size() > 0) {
			newGroups.add(newGroup);
		}
		
		return newGroups;
	}
	
	public void determineBaseToReport(HashMap<String,String> refSeq) {
		/*This method figures out if the deletion originated from a 'T'. If it does, it sets
		The position to the first 'T' and sets the baseFlag. The baseFlag lets the user know
		how the software determined the deletion originated from a 'T'*/
		String sequence = refSeq.get(maxPos.getChrom()).substring(maxPos.getPos(),maxPos.getPos() + hpLength-1).toUpperCase();
		
		int beginRange = maxPos.getPos()-25;
		int endRange = maxPos.getPos()+25;
		if (beginRange < 0) {
			beginRange = 0;
		}
		if (endRange > refSeq.get(maxPos.getChrom()).length()) {
			endRange = refSeq.get(maxPos.getChrom()).length()-1;
		}
		
		String sequence2 = refSeq.get(maxPos.getChrom()).substring(beginRange,endRange).toUpperCase();
		
		
		maxPos.setPredPos(-1);
		String direction = maxPos.getDirection().substring(0,1);
		Matcher mF = pF.matcher(sequence);
		Matcher mR = pR.matcher(sequence);
		
		int altBase = 0;
		for (Position p: posList) {
			String indDir = p.getDirection().substring(0,1);
			if (indDir.equals("+") && p.getBase() == 'T') {
				altBase = p.getPos();
			} else if (indDir.equals("-") && p.getBase() == 'A') {
				altBase = p.getPos();
			}
		}
		
		if (direction.equals("+")) {
			if (maxPos.getBase() == 'T') {
				maxPos.setBaseFlag("FirstDelBase");
				maxPos.setPredPos(maxPos.getPos());
			}  else if (altBase != 0) {
				maxPos.setBaseFlag("AltPositionBase");
				maxPos.setPredPos(altBase);
				maxPos.setBase('T');
			} else if (mF.find() && mF.group(0).indexOf('T') != -1) {
				maxPos.setPredPos(maxPos.getPos() + mF.group(0).indexOf("T"));
				maxPos.setBaseFlag("DwnStrmDelBase");
				maxPos.setBase('T');
			} else {
				int bestDistance = 50;
				int bestPos = -1;
				int currPos = maxPos.getPos() - 25;
				for (char c: sequence2.toCharArray()) {
					if (c == 'T' ) {
						int currDistance = Math.abs(currPos - maxPos.getPos());
						if (currDistance < bestDistance) {
							bestDistance = currDistance;
							bestPos = currPos;
						}
					}
					currPos++;
				}
				maxPos.setBase('T');
				maxPos.setPredPos(bestPos);
				maxPos.setBaseFlag("Neighborhood");
			}
		} else {
			if (maxPos.getBase() == 'A') {
				maxPos.setBaseFlag("FirstDelBase");
				maxPos.setPredPos(maxPos.getPos());
				maxPos.setBase('T');
			} else if (altBase != 0) {
				maxPos.setBaseFlag("AltPositionBase");
				maxPos.setPredPos(altBase);
				maxPos.setBase('T');
			} else if (mR.find() && mR.group(0).indexOf("A") != -1) {
				maxPos.setPredPos(maxPos.getPos() + mR.group(0).indexOf("A"));
				maxPos.setBaseFlag("DwnStrmDelBase");
				maxPos.setBase('T');
			} else  {
				int bestDistance = 50;
				int bestPos = -1;
				int currPos = maxPos.getPos() - 25;
				for (char c: sequence2.toCharArray()) {
					if (c == 'A' ) {
						int currDistance = Math.abs(currPos - maxPos.getPos());
						if (currDistance < bestDistance) {
							bestDistance = currDistance;
							bestPos = currPos;
						}
					}
					currPos++;
				}
				maxPos.setBase('T');
				maxPos.setPredPos(bestPos);
				maxPos.setBaseFlag("Neighborhood");
			}
		}
	}

	public int getLastPos() {
		return lastPos;
	}

	public double getBsPval() {
		return bsPval;
	}

	public Position getMaxPos() {
		return maxPos;
	}

	public double getNbsPval() {
		return nbsPval;
	}

	public String getFilterFlag() {
		return filterFlag;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getBiotype() {
		return biotype;
	}

	public ArrayList<Position> getPosList() {
		return posList;
	}

	
	
	


	

	
	
}

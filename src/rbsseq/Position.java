package rbsseq;

public class Position {		
	//Basic information
	private String chrom;
	private String direction;
	private int pos;
	private int predPos;
	private char base;
	private int bsDepth;
	private int bsDel;
	private int nbsDepth;
	private int nbsDel;
	private double bsPer = -1;
	private double nbsPer = -1;
	private String baseFlag;
	
	public Position(String chrom, int pos, String direction, String base, int bsDepth, int bsDel, int nbsDepth, int nbsDel) {
		this.chrom = chrom;
		this.pos = pos;
		this.direction = direction;
		this.base = base.toUpperCase().charAt(0);
		this.bsDepth = bsDepth;
		this.bsDel = bsDel;
		this.nbsDepth = nbsDepth;
		this.nbsDel = nbsDel;
		
		if (this.nbsDepth == 0) {
			this.nbsPer = -1;
		} else {
			this.nbsPer = (float)this.nbsDel / this.nbsDepth;
		}
		this.bsPer = (float)this.bsDel / this.bsDepth;
	}

	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public String getDirection() {
		return direction;
	}

	public void setDirection(String direction) {
		this.direction = direction;
	}

	public int getPos() {
		return pos;
	}

	public void setPos(int pos) {
		this.pos = pos;
	}

	public int getPredPos() {
		return predPos;
	}

	public void setPredPos(int predPos) {
		this.predPos = predPos;
	}

	public char getBase() {
		return base;
	}

	public void setBase(char base) {
		this.base = base;
	}

	public int getBsDepth() {
		return bsDepth;
	}

	public void setBsDepth(int bsDepth) {
		this.bsDepth = bsDepth;
	}

	public int getBsDel() {
		return bsDel;
	}

	public void setBsDel(int bsDel) {
		this.bsDel = bsDel;
	}

	public int getNbsDepth() {
		return nbsDepth;
	}

	public void setNbsDepth(int nbsDepth) {
		this.nbsDepth = nbsDepth;
	}

	public int getNbsDel() {
		return nbsDel;
	}

	public void setNbsDel(int nbsDel) {
		this.nbsDel = nbsDel;
	}

	public double getBsPer() {
		return bsPer;
	}

	public void setBsPer(double bsPer) {
		this.bsPer = bsPer;
	}

	public double getNbsPer() {
		return nbsPer;
	}

	public void setNbsPer(double nbsPer) {
		this.nbsPer = nbsPer;
	}

	public String getBaseFlag() {
		return baseFlag;
	}

	public void setBaseFlag(String baseFlag) {
		this.baseFlag = baseFlag;
	}
	
	
}

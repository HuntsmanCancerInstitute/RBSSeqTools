package rbsseq;

class Feature  {
	private int start;
	private int end;
	private String biotype;
	private String name;
	
	public Feature(int start, int end, String name, String biotype) {
		this.start = start;
		this.end = end;
		this.name = name;
		this.biotype = biotype;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public String getBiotype() {
		return biotype;
	}

	public void setBiotype(String biotype) {
		this.biotype = biotype;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	

}

package rbsseq;

import java.util.Comparator;

public class PositionComparator implements Comparator<PositionGroup> {
	@Override
	public int compare(PositionGroup pg1, PositionGroup pg2) {
		Double pval1 = new Double(pg1.getBsPval());
		Double pval2 = new Double(pg2.getBsPval());
		return pval2.compareTo(pval1);
		
	}
}

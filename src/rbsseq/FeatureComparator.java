package rbsseq;

import java.util.Comparator;

public class FeatureComparator implements Comparator<Feature> {
    @Override
    public int compare(Feature o1, Feature o2) {
    	Integer start1 = new Integer(o1.getStart());
    	Integer start2 = new Integer(o2.getStart());
    	return start1.compareTo(start2);
    }
}
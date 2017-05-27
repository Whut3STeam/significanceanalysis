package edu.whut.significance.dataset;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by SunMing on 2017/5/23.
 */
public class ResultData {
    private Set<Region> regionSet;

    public ResultData(){
        regionSet=new HashSet<>();
    }

    public Set<Region> getRegionSet() {
        return regionSet;
    }

    public void setRegionSet(Set<Region> regionSet) {
        this.regionSet = regionSet;
    }

    public void addRegion(int start, int end ){
        regionSet.add(new Region(start, end));
    }
}

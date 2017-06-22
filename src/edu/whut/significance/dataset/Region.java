package edu.whut.significance.dataset;

import org.w3c.dom.ranges.Range;

/**
 * Created by SunMing on 2017/5/23.
 */
public class Region implements Comparable<Region>{
    private int startId;
    private int endId;
//    private int length;


    public Region() {
    }

    public Region(int startId, int endId) {
        this.startId = startId;
        this.endId = endId;
    }

    public Region(){}
    public Region(int startId,int endId){
        this.startId=startId;
        this.endId=endId;
        this.length=endId-startId+1;
    }
    public int getStartId(){
        return startId;
    }
    public int getEndId(){
        return endId;
    }
    public int getLength(){
        return endId - startId +1 ;
    }

    public void setStartId(int startId){
        this.startId=startId;
    }
    public void setEndId(int endId){
        this.endId=endId;
    }

    @Override
    public int compareTo(Region o) {
        return Integer.compare(this.getStartId(), o.getStartId());
    }

//    public void setLength(int length){
//        this.length=length;
//    }
}

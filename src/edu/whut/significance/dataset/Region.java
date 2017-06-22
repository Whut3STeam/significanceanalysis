package edu.whut.significance.dataset;

/**
 * Created by SunMing on 2017/5/23.
 */
public class Region {
    private int startId;
    private int endId;
    private int length;

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
        return length;
    }

    public void setStartId(int startId){
        this.startId=startId;
    }
    public void setEndId(int endId){
        this.endId=endId;
    }
    public void setLength(int length){
        this.length=length;
    }
}

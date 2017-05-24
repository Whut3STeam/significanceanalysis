package edu.whut.significance.dataset;

import java.util.Map;

/**
 * Created by SunMing on 2017/5/23.
 */
public class Index {
    private int id;
    private String chrName;
    private int pos;
    private String snpName;

    public int getId(){
        return id;
    }
    public String getChrName(){
        return chrName;
    }
    public int getPos(){
        return pos;
    }
    public String getSNPName(){
        return snpName;
    }

    public void setId(int id){
        this.id=id;
    }
    public void setChrName(String chrName){
        this.chrName=chrName;
    }
    public void setPos(int pos){
        this.pos=pos;
    }
    public void setSNPName(String snpName){
        this.snpName=snpName;
    }
}

package edu.whut.significance.methods;

/**
 * Created by SunMing on 2017/5/24.
 */
public class GlobalParameters {
    private double ampThreshold = 0.322;
    //private double ampThreshold=0.1;
    private double delThreshold = -0.415;

    public double getAmpThreshold() {
        return ampThreshold;
    }

    public void setAmpThreshold(double ampThreshold) {
        this.ampThreshold = ampThreshold;
    }

    public double getDelThreshold() {
        return delThreshold;
    }

    public void setDelThreshold(double delThreshold) {
        this.delThreshold = delThreshold;
    }
}

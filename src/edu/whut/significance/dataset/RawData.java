package edu.whut.significance.dataset;

import org.apache.commons.math3.linear.RealMatrix;

import java.util.Map;

/**
 * Created by SunMing on 2017/5/23.
 */
public class RawData {
    private RealMatrix dataMatrix;

    public RealMatrix getDataMatrix() {
        return dataMatrix;
    }

    public void setDataMatrix(RealMatrix dataMatrix) {
        this.dataMatrix = dataMatrix;
    }
}

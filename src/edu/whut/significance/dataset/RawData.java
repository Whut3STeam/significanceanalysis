package edu.whut.significance.dataset;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Map;

/**
 * Created by SunMing on 2017/5/23.
 */
public class RawData {
    private RealMatrix m_data;

    public RealMatrix getDataMatrix() {
        return m_data;
    }

    public void setData(RealMatrix data) {
        this.m_data = data;
    }

    public void setData(double[][] data){
        m_data = new BlockRealMatrix(data);
    }

    public void setRow(int i, double[] aRow){
        m_data.setRow(i, aRow);
    }

}

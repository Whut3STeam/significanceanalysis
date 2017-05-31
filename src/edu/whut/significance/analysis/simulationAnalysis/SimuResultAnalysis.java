package edu.whut.significance.analysis.simulationAnalysis;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.GlobalParameters;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/24.
 */
public class SimuResultAnalysis {
    private Logger m_log;

    private RawData rawData;
    private RealMatrix rawMatrix;
    private ResultData resultData;
    private int probeNum;
    private int sampleNum;
    private double ampThreshold;
    private double delThreshold;
    private boolean[] PArray;//true则实际为正，false则实际为负
    private boolean[] TArray;//true则预测为正，false则预测为负
    private int P;
    private int T;
    private int TP;//正确找到包含在模拟CNA区域的一条CNA探针
    private int FP;//找到的CNA探针不包含在模拟CNA区域
    private int FN;//没有被找到的包含在模拟CNA区域的探针
    private double precision;//准确率
    private double recall;//召回率
    private double FMeasure;//F值

    public SimuResultAnalysis(RawData rawData, ResultData resultData) {
        this.rawData = rawData;
        this.resultData = resultData;
        rawMatrix = rawData.getDataMatrix();
        GlobalParameters globalParameters = new GlobalParameters();
        ampThreshold = globalParameters.getAmpThreshold();
        delThreshold = globalParameters.getDelThreshold();
        probeNum = rawMatrix.getColumnDimension();
        sampleNum = rawMatrix.getRowDimension();
        PArray = new boolean[probeNum];
        TArray = new boolean[probeNum];
        P = 0;
        T = 0;
        analysis();
    }

    public double getFMeasure() {
        return FMeasure;
    }

    public void analysis() {
        //根据原始数据初始化P
//        for (int i = 0; i < probeNum; i++) {
//            for (int j = 0; j < sampleNum; j++) {
//                if (rawMatrix.getEntry(j, i) > 0.55 || rawMatrix.getEntry(j, i) < delThreshold) {
//                    PArray[i] = true;
//                    P++;
//                    break;
//                }
//            }
//        }
        for (int i = 401; i <= 600; i++){
            PArray[i] = true;
            P++;
        }

        //根据结果数据初始化T
        for (Region resultRegion : resultData.getRegionSet()) {
            for (int i = resultRegion.getStartId(); i <= resultRegion.getEndId(); i++) {
                TArray[i] = true;
                T++;
            }
        }

        //给TP，FP，FN赋值
        for (int i = 0; i < probeNum; i++) {
            if (PArray[i] && TArray[i]) {
                TP++;
            } else if ((!PArray[i]) && TArray[i]) {
                FP++;
            } else if (PArray[i] && (!TArray[i])) {
                FN++;
            }
        }

        //计算准确率，召回率
        precision = (double) TP / (TP + FP);
        recall = (double) TP / (TP + FN);

        //计算F值
        FMeasure = 2 * precision * recall / (precision + recall);
    }
}

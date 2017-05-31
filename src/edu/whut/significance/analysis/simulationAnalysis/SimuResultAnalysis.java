package edu.whut.significance.analysis.simulationAnalysis;

import com.alibaba.fastjson.JSON;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import edu.whut.significance.dataset.ExampleJ;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.GlobalParameters;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
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
    private ExampleJ exampleJ;
    private boolean[] PArray;//true��ʵ��Ϊ����false��ʵ��Ϊ��
    private boolean[] TArray;//true��Ԥ��Ϊ����false��Ԥ��Ϊ��
    private int P;
    private int T;
    private int TP;//��ȷ�ҵ�������ģ��CNA�����һ��CNA̽��
    private int FP;//�ҵ���CNA̽�벻������ģ��CNA����
    private int FN;//û�б��ҵ��İ�����ģ��CNA�����̽��
    private double precision;//׼ȷ��
    private double recall;//�ٻ���
    private double FMeasure;//Fֵ

    public SimuResultAnalysis(RawData rawData, ResultData resultData, String filePath) {
        m_log = Logger.getLogger("significanceAnalysis");

        try {
            File infile = new File(filePath);
            InputStreamReader isr = new InputStreamReader(new FileInputStream(infile));
            BufferedReader br = new BufferedReader(isr);
            String jsonString = "";
            String line;
            while ((line = br.readLine()) != null) {
                jsonString += line;
            }
            exampleJ = JSON.parseObject(jsonString, ExampleJ.class);
        } catch (Exception e) {
            e.printStackTrace();
        }

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
        //����json�ļ��ó���exampleJ��ʼ����ʵ����
        int midPos, width;
        RangeSet<Integer> range = TreeRangeSet.create();
        for (ExampleJ.Sample sample : exampleJ.getSamples()) {
            midPos = sample.getWindows().getMidPos();
            width = sample.getWindows().getWidth();
            range.add(Range.closed(midPos - width / 2, midPos + width / 2));
        }

        //����ԭʼ���ݳ�ʼ��P
        for (int i = 0; i < probeNum; i++) {
            if (range.contains(i)) {
                PArray[i] = true;
                P++;
            }
        }

        //���ݽ�����ݳ�ʼ��T
        for (Region resultRegion : resultData.getRegionSet()) {
            for (int i = resultRegion.getStartId(); i <= resultRegion.getEndId(); i++) {
                TArray[i] = true;
                T++;
            }
        }

        //��TP��FP��FN��ֵ
        for (int i = 0; i < probeNum; i++) {
            if (PArray[i] && TArray[i]) {
                TP++;
            } else if ((!PArray[i]) && TArray[i]) {
                FP++;
            } else if (PArray[i] && (!TArray[i])) {
                FN++;
            }
        }

        //����׼ȷ�ʣ��ٻ���
        precision = (double) TP / (TP + FP);
        recall = (double) TP / (TP + FN);
        m_log.info(String.format("׼ȷ�ʣ�%.6f\t �ٻ��ʣ�%.6f", precision, recall));

        //����Fֵ
        FMeasure = 2.0 * precision * recall / (precision + recall);
    }
}

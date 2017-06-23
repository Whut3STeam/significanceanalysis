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
    private boolean[] PArray;//true则实际为正，false则实际为负
    private boolean[] TArray;//true则预测为正，false则预测为负
    private int P;
    private int T;
    private int TP;//正确找到包含在模拟CNA区域的一条CNA探针
    private int FP;//找到的CNA探针不包含在模拟CNA区域
    private int FN;//没有被找到的包含在模拟CNA区域的探针
    private int TN;
    private double precision;//准确率
    private double recall;//召回率(TPR)
    private double FPR;
    private double FMeasure;//F值

    public SimuResultAnalysis(RawData rawData,ResultData resultData,String filePath){
        m_log = Logger.getLogger("significanceAnalysis");

        try {
            File infile = new File(filePath);
            InputStreamReader isr = new InputStreamReader(new FileInputStream(infile));
            BufferedReader br = new BufferedReader(isr);
            String jsonString="";
            String line;
            while ((line=br.readLine())!=null){
                jsonString+=line;
            }
            exampleJ= JSON.parseObject(jsonString,ExampleJ.class);
        }
        catch (Exception e){
            e.printStackTrace();
        }

        this.rawData=rawData;
        this.resultData=resultData;
        rawMatrix=rawData.getDataMatrix();
        GlobalParameters globalParameters=new GlobalParameters();
        ampThreshold=globalParameters.getAmpThreshold();
        delThreshold=globalParameters.getDelThreshold();
        probeNum=rawMatrix.getColumnDimension();
        sampleNum=rawMatrix.getRowDimension();
        PArray =new boolean[probeNum];
        TArray =new boolean[probeNum];
        P=0;
        T=0;
        analysis();
    }

    public double getFMeasure() {
        return FMeasure;
    }

    public void analysis(){
        //根据json文件得出的exampleJ初始化真实区域range
        int midPos,width;
        RangeSet<Integer> range= TreeRangeSet.create();
        for(ExampleJ.Sample sample:exampleJ.getSamples()){
            for (ExampleJ.Sample.Window window:sample.getWindows()){
                midPos=window.getMidPos();
                width=window.getWidth();
                range.add(Range.closed(midPos-width/2,midPos+width/2));
            }
        }

        //根据原始数据初始化P
        for(int i=0;i<probeNum;i++){
            if(range.contains(i)) {
                PArray[i] = true;
                P++;
            }
        }

        //根据结果数据初始化T
        for(Region resultRegion:resultData.getRegionSet()){
            for(int i=resultRegion.getStartId();i<=resultRegion.getEndId();i++){
                TArray[i]=true;
                T++;
            }
        }

        //给TP，FP，FN赋值
        for(int i=0;i<probeNum;i++){
            if(PArray[i]&&TArray[i]){
                TP++;
            }
            else if((!PArray[i])&&TArray[i]){
                FP++;
            }
            else if(PArray[i]&&(!TArray[i])){
                FN++;
            }
            else if(!PArray[i]&&(!TArray[i])){
                TN++;
            }
        }

        //计算准确率，召回率
        precision=(double)TP/(TP+FP);
        recall=(double)TP/(TP+FN);
        FPR=(double)FP /(FP + TN);

        m_log.info("\n");
        m_log.info("准确率："+precision);
        m_log.info("召回率(TPR)："+recall);
        m_log.info("FPR:"+FPR);

        //计算F值
        FMeasure=2*precision*recall/(precision+recall);
        m_log.info("F值："+FMeasure);
    }
}

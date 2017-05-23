package edu.whut.significance.dataset;

import org.apache.commons.math3.random.RandomDataGenerator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/16.
 */
public class SampleGEN {
    private  String filePath;
    private RandomDataGenerator rdg;
    private Logger m_log;

    public SampleGEN(){
        m_log = Logger.getLogger("significanceDetection");
        rdg=new RandomDataGenerator();
        this.filePath="data//simulatedData//sampleGenerated.txt";
    }

    public SampleGEN(String filePath){
        m_log = Logger.getLogger("significanceDetection");
        rdg=new RandomDataGenerator();
        this.filePath=filePath;
    }

    public void generate(int number,int midPos,int width,double LRRVal,double sigma,double alpha,double beta){
        double[] singleSample=new double[1000];
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
            for (int i = 0; i < number; i++) {
                singleSample = singleSamplegenerate((int) (midPos + width*rdg.nextGaussian(0, beta)),
                        (int) (width + rdg.nextGaussian(1, alpha * width)),
                        LRRVal, sigma);
                String newLine="";
                for(int j=0;j<1000-1;j++){
                    newLine+=singleSample[j]+"\t";
                }
                newLine+=singleSample[1000-1];
                bw.write(newLine);
                bw.newLine();
            }
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public double[] singleSamplegenerate(int midPos,int width,double LRRVal,double sigma){
        double[] singleSample=new double[1000];
        for(int i=0;i<midPos-width/2;i++){
            singleSample[i]=rdg.nextGaussian(0,sigma);
        }
        for(int i=midPos-width/2;i<midPos+width/2;i++){
            singleSample[i]=LRRVal+rdg.nextGaussian(0,sigma);
        }
        for(int i=midPos+width/2;i<1000;i++){
            singleSample[i]=rdg.nextGaussian(0,sigma);
        }
        return singleSample;
    }
}

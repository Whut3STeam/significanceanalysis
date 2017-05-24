package edu.whut.significance.dataset;

import com.alibaba.fastjson.JSON;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/16.
 */
public class SampleGEN {
    private  String filePath;
    private RandomDataGenerator rdg;
    private Logger m_log;
    private ExampleJ exampleJ;

    public SampleGEN(){
        m_log = Logger.getLogger("significanceAnalysis");
        rdg=new RandomDataGenerator();
        this.filePath="data//simulatedData//sampleGenerated.json";
        exampleJ=new ExampleJ();
        exampleJ.samples=new ArrayList<>();
    }

    public SampleGEN(String filePath){
        m_log = Logger.getLogger("significanceDetection");
        rdg=new RandomDataGenerator();
        this.filePath=filePath;
        exampleJ=new ExampleJ();
        exampleJ.samples=new ArrayList<>();
    }

    public void generate(int number,int midPos,int width,double LRRVal,double sigma,double alpha,double beta){
        ExampleJ.Sample sample=new ExampleJ.Sample();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMddHHmmss");
        sample.id=df.format(new Date());
        sample.length=1000;
        sample.value=0;
        sample.count=number;
        sample.sigma=sigma;
        sample.windows=new ExampleJ.Sample.Windows();
        sample.windows.midPos=midPos;
        sample.windows.width=width;
        sample.windows.value=LRRVal;
        sample.windows.alpha=alpha;
        sample.windows.beta=beta;
        sample.data=new double[number][1000];

        double[] singleSample=new double[1000];
        for(int i=0;i<number;i++) {
            sample.data[i] = singleSamplegenerate((int) (midPos + width * rdg.nextGaussian(0, beta)),
                    (int) (width + rdg.nextGaussian(1, alpha * width)),
                    LRRVal, sigma);
        }
        exampleJ.samples.add(sample);
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

    public void saveSamples(){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
            String jsonString= JSON.toJSONString(exampleJ);

            //ExampleJ exampleJ1=JSON.parseObject(jsonString,ExampleJ.class);

            bw.write(jsonString);
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}

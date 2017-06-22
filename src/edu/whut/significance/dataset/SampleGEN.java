package edu.whut.significance.dataset;

import com.alibaba.fastjson.JSON;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Random;
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

    public void generate(int length,int sCount,double baselineVal,double sigma,
                         int breakPoint, double baselineInrcement, double baselineFluctuation,
                         int passengerCount,double passengerVal,
                         int winCount,int[] midPos,int[] width,double[] LLRVal,double[] alpha,double[] beta){
        //��ʼ��json����
        ExampleJ.Sample sample=new ExampleJ.Sample();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMddHHmmss");
        sample.id=df.format(new Date());
        sample.length=length;//��������
        sample.sCount =sCount;//��������
        sample.baselineValue =baselineVal;//����ֵ
        sample.sigma=sigma;//ǿ������
        sample.breakPoint=breakPoint;//���߷ָ��
        sample.baselineInrcement=baselineInrcement;//�ڶ�����������ڵ�һ�����ߵı仯ֵ
        sample.baselineFluctuation=baselineFluctuation;//����Ĳ����ٷֱ�
        sample.passengerCount=passengerCount;//�˿ͻ�������
        sample.passengerVal=passengerVal;//�˿͵�ȡֵ
        sample.winCount=winCount;//��������
        sample.data=new double[sCount][length];

        ExampleJ.Sample.Window[] windows=new ExampleJ.Sample.Window[winCount];
        for(int i=0;i<winCount;i++){
            ExampleJ.Sample.Window window=new ExampleJ.Sample.Window();
            window.midPos=midPos[i];
            window.width=width[i];
            window.LRRVal=LLRVal[i];
            window.alpha=alpha[i];
            window.beta=beta[i];

            windows[i]=window;
        }
        sample.windows=windows;

        double[] singleSample=new double[sample.length];

        //��������
        for(int i=0;i<sCount;i++){
            //���ɴ��ڵ��е�Ϳ��ֵ����
            int[] tempMidPos=new int[winCount];
            int[] tempWidth=new int[winCount];
            for(int j=0;j<winCount;j++){
                tempMidPos[j]=(int) (midPos[j] + width[j] * rdg.nextGaussian(0, beta[j]));
                tempWidth[j]=(int) (width[j] + rdg.nextGaussian(1, alpha[j] * width[j]));
            }

            //���մ��ڲ�������һ������
            sample.data[i] = singleSamplegenerate(length,baselineVal,passengerCount,passengerVal,
                    breakPoint,baselineInrcement,baselineFluctuation,
                    winCount,tempMidPos,tempWidth,LLRVal,sigma);


        }
        exampleJ.samples.add(sample);
    }


    public double[] singleSamplegenerate(int length,double baselineVal,int passengerCount,double passengerVal,
                                         int breakPoint,double baselineInrcement,double baselineFluction,
                                         int winCount,int[] midPos,int[] width,double[] LRRVal,double sigma){
        double[] singleSample=new double[length];

        //��ӻ���ֵ
        for(int i=0;i<length;i++){
            singleSample[i]=baselineVal;
        }

        //���passengers
        int passengerStart,passengerEnd;
        for(int i=0;i<passengerCount;i++){
            passengerStart=rdg.nextInt(0,length);
            passengerEnd=passengerStart+50;
            passengerEnd=passengerEnd>length?length:passengerEnd;

            for(int j=passengerStart;j<passengerEnd;j++){
                singleSample[j]=passengerVal;
            }
        }

        //���drivers
        for(int i=0;i<winCount;i++){
            for(int j=midPos[i]-width[i]/2;j<midPos[i]+width[i]/2;j++){
                singleSample[j]=LRRVal[i];
            }
        }

        //��ӵڶ�������
        int range = (int)(breakPoint*baselineFluction);//������Χ
        breakPoint+=rdg.nextInt(-range,range);
        breakPoint=breakPoint<0?0:(breakPoint>(length-1)?(length-1):breakPoint);
        for(int i=breakPoint;i<length;i++){
            singleSample[i]+=baselineInrcement;
        }

        //���ǿ������
        if(sigma>0){
            for(int i=0;i<length;i++){
                singleSample[i]+=rdg.nextGaussian(0,sigma);
            }
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

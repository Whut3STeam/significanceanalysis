package edu.whut.significance.dataset;

import com.alibaba.fastjson.JSON;
import org.apache.commons.math3.linear.BlockRealMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

/**
 * Created by SunMing on 2017/5/23.
 */
public class Reader {
    public static void readSimulationData(RawData rawData,String filePath){
        ExampleJ exampleJ=new ExampleJ();
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

        int num=0;
        for(ExampleJ.Sample sample:exampleJ.samples){
            num+=sample.sCount;
        }

        rawData.setDataMatrix(new BlockRealMatrix(num+1,exampleJ.samples.get(0).length));
        //写第一行
        double[] ids=new double[exampleJ.samples.get(0).length];
        for(int i=0;i<ids.length;i++){
            ids[i]=i;
        }
        rawData.getDataMatrix().setRow(0,ids);

        //写数据矩阵
        int i=1;
        for(ExampleJ.Sample sample:exampleJ.samples){
            for(int j=0;j<sample.data.length;j++){
                rawData.getDataMatrix().setRow(i,sample.data[j]);
                i++;
            }
        }
    }

    public static void readTrueData(RawData rawData,String filePath){

    }
}

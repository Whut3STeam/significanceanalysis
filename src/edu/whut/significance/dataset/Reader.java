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
    public static void readSimulationData(RawData rawData, String filePath) {
        ExampleJ exampleJ = new ExampleJ();
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

        int num = 0;
        for (ExampleJ.Sample sample : exampleJ.samples) {
            num += sample.sCount;
        }

        rawData.setData(new BlockRealMatrix(num, exampleJ.samples.get(0).length));
        //Ð´Êý¾Ý¾ØÕó
        int i = 0;
        for (ExampleJ.Sample sample : exampleJ.samples) {
            for (int j = 0; j < sample.data.length; j++) {
                rawData.setRow(i, sample.data[j]);
                i++;
            }
        }
    }

    public static void readTrueData(RawData rawData, String filePath) {

    }
}

package edu.whut.significance.methods;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.ResultData;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by SunMing on 2017/6/1.
 */
public class RandomForest {
    private RawData rawData;
    private ResultData resultData;
    private AbstractSig saic;

    public RandomForest(RawData rawData,ResultData resultData){
        this.rawData=rawData;
        this.resultData=resultData;

        process();
    }

    public void process(){
        List<RawData> rawDataList=new ArrayList<>();
        List<ResultData> resultDataList=new ArrayList<>();
        sample(rawDataList);

        saic=new SAIC();
        for(RawData tempRawData:rawDataList){
            ResultData tempResultData=new ResultData();
            saic.preprocess(tempRawData);
            saic.process(tempResultData);
            resultDataList.add(tempResultData);
        }

        vote(resultDataList);
    }

    public void sample(List<RawData> rawDataList){
        int sampleNum=rawData.getDataMatrix().getRowDimension()-1;
        
    }

    public void vote(List<ResultData> resultDataList){

    }

    static class Parameters{
        static int sampleFrequency=100;
        static int sampleSize=10;
    }
}

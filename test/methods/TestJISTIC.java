package methods;

import edu.whut.significance.analysis.simulationAnalysis.SimuResultAnalysis;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Reader;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.AbstractSig;
import edu.whut.significance.methods.JISTIC;
import org.junit.Test;

/**
 * Created by SunMing on 2017/6/20.
 */
public class TestJISTIC {
    @Test
    public void test(){
        RawData rawData=new RawData();
        ResultData resultData=new ResultData();
        String filePath="data//simulatedData//20170622113917.json";
        Reader.readSimulationData(rawData,filePath);
        AbstractSig jistic=new JISTIC();
        jistic.preprocess(rawData);
        jistic.process(resultData);
        System.out.println("FÖµ£º"+new SimuResultAnalysis(rawData,resultData,filePath).getFMeasure());
        System.out.println("OK");
    }
}

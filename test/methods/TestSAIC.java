package methods;

import edu.whut.significance.analysis.simulationAnalysis.SimuResultAnalysis;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Reader;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.AbstractSig;
import edu.whut.significance.methods.SAIC;
import org.junit.Test;

/**
 * Created by SunMing on 2017/5/24.
 */
public class TestSAIC {
    @Test
    public void test(){
        RawData rawData=new RawData();
        ResultData resultData=new ResultData();
        String filePath="data//simulatedData//20170531155103.json";
        Reader.readSimulationData(rawData,filePath);
        AbstractSig saic=new SAIC();
        saic.preprocess(rawData);
        saic.process(resultData);
        System.out.println("FÖµ£º"+new SimuResultAnalysis(rawData,resultData,filePath).getFMeasure());
    }
}

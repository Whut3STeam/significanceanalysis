package methods;

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
        Reader.readSimulationData(rawData,"data//simulatedData//20170524111638.json");
        AbstractSig saic=new SAIC();
        saic.preprocess(rawData);
        saic.process(resultData);
        System.out.println("OK");
    }
}

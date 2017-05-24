package dataset;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Reader;
import org.junit.Test;

/**
 * Created by SunMing on 2017/5/24.
 */
public class TestReader {
    @Test
    public void test(){
        RawData rawData=new RawData();
        Reader.readSimulationData(rawData,"data//simulatedData//20170524111638.json");
        System.out.println("OK");
    }
}

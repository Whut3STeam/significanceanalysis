package dataset;

import edu.whut.significance.dataset.SampleGEN;
import org.junit.Test;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Created by SunMing on 2017/5/16.
 */
public class TestSampleGenerator {
    @Test
    public void test(){
        String filePath=String.format("data//simulatedData//%s.json",
                new SimpleDateFormat("yyyyMMddHHmmss").format(new Date()));
        SampleGEN sg=new SampleGEN(filePath);
        //SampleGEN sg=new SampleGEN();
        int winCount=3;
        //以下数组的长度必须为winCount的值
        int[] midPos={400,1600,2500};
        int[] width={200,150,100};
        double[] LRRVal={0.58,0.6,0.62};
        double[] alpha={0.1,0.1,0.1};
        double[] beta={0.1,0.1,0.1};
        sg.generate(3000,20,0,0.1,
                0,0,0,
                10,0.35,
                3,midPos,width,LRRVal,alpha,beta);
        sg.saveSamples();
    }
}

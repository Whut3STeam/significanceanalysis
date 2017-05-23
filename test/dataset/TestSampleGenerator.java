package dataset;

import edu.whut.significance.dataset.SampleGEN;
import org.junit.Test;

/**
 * Created by SunMing on 2017/5/16.
 */
public class TestSampleGenerator {
    @Test
    public void test(){
        SampleGEN sg=new SampleGEN();
        /*double[] singleSample=sg.singleSamplegenerate(200,200,0.58,0.25);
        for(int i=0;i<1000;i++){
            if(i%100==0)
                System.out.println();
            System.out.print(singleSample[i]+"\t");
        }*/
        sg.generate(4,500,200,0.58,0.25,0.25,0.25);
    }
}

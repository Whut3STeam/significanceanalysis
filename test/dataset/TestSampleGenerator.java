package dataset;

import edu.whut.significance.dataset.SampleGEN;
import org.junit.Test;

import java.text.DateFormat;
import java.text.FieldPosition;
import java.text.ParsePosition;
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
        sg.generate(4,500,200,0.58,0.25,0.25,0.25);
        sg.generate(3,400,200,-1,0.25,0.25,0.25);
        sg.saveSamples();
    }
}

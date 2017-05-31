package methods;

import edu.whut.significance.analysis.simulationAnalysis.SimuResultAnalysis;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Reader;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.AbstractSig;
import edu.whut.significance.methods.SAIC;
import edu.whut.significance.util.BioLogger;
import org.junit.Test;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/24.
 */
    public class TestSAIC {
    @Test
    public void test(){
        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String filename = String.format("Result_%s.log", df.format(date));
        new BioLogger("data",filename);
        date = null;
        df = null;

        RawData rawData = new RawData();
        ResultData resultData = new ResultData();
        Reader.readSimulationData(rawData, "data//simulatedData//20170531160338.json");
        Logger.getLogger("significanceAnalysis").info(String.format("Input Data >>> Row = %d, Col = %d",
                rawData.getDataMatrix().getRowDimension(), rawData.getDataMatrix().getColumnDimension()));
        ;
        AbstractSig saic = new SAIC();
        saic.preprocess(rawData);
        saic.process(resultData);
        Logger.getLogger("significanceAnalysis").info("FÖµ£º" + new SimuResultAnalysis(rawData, resultData).getFMeasure());
    }
}

package methods;

import edu.whut.significance.analysis.simulationAnalysis.SimuResultAnalysis;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Reader;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.RandomForest;
import edu.whut.significance.util.BioLogger;
import org.junit.Test;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/6/1.
 */
public class TestRandomForest {
    @Test
    public void test(){
        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String filename = String.format("Result_%s.log", df.format(date));
        new BioLogger("data", filename);
        date = null;
        df = null;

        String filePath = "data/simulatedData/20170622183459.json";
        RawData rawData = new RawData();
        ResultData resultData = new ResultData();
        Reader.readSimulationData(rawData, filePath);

        Logger m_log = Logger.getLogger("significanceAnalysis");
        m_log.info(String.format("Input Data >>> Row = %d, Col = %d",
                rawData.getDataMatrix().getRowDimension(), rawData.getDataMatrix().getColumnDimension()));

        RandomForest randomForest=new RandomForest(rawData,resultData);
        new SimuResultAnalysis(rawData, resultData, filePath).getFMeasure();
    }
}

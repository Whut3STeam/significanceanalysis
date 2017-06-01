package edu.whut.significance.methods;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/6/1.
 */
public class RandomForest {
    int sampleNum;
    int probeNum;
    private RawData rawData;
    private ResultData resultData;
    private AbstractSig saic;
    private boolean enableDedugeInfo = false;
    private Logger m_log;

    public RandomForest(RawData rawData, ResultData resultData) {
        this.rawData = rawData;
        this.resultData = resultData;
        sampleNum = rawData.getDataMatrix().getRowDimension() - 1;
        probeNum = rawData.getDataMatrix().getColumnDimension();
        m_log = Logger.getLogger("significanceAnalysis");

        process();
    }

    public void process() {
        List<RawData> rawDataList = new ArrayList<>();
        List<ResultData> resultDataList = new ArrayList<>();
        sample(rawDataList);

        for (RawData tempRawData : rawDataList) {
            saic = new SAIC();
            ResultData tempResultData = new ResultData();
            saic.preprocess(tempRawData);
            saic.process(tempResultData);
            if (tempResultData.getRegionSet().size() > 0)
                resultDataList.add(tempResultData);
        }

        vote(resultDataList);
    }

    public void sample(List<RawData> rawDataList) {
        ArrayList<Integer> randomList = new ArrayList<>();
        for (int i = 0; i < sampleNum; i++) {
            randomList.add(i + 1);
        }

        for (int i = 0; i < Parameters.sampleFrequency; i++) {
            Collections.shuffle(randomList);
            RawData tempRawData = new RawData();
            RealMatrix tempRawMatrix = new BlockRealMatrix(Parameters.sampleSize, probeNum);

            for (int j = 0; j < Parameters.sampleSize; j++) {
                tempRawMatrix.setRow(j, rawData.getDataMatrix().getRow(randomList.get(j)));
            }

            tempRawData.setData(tempRawMatrix);
            rawDataList.add(tempRawData);
        }
    }

    public void vote(List<ResultData> resultDataList) {
        int[] voteNum = new int[probeNum];
        for (ResultData tempResultData : resultDataList) {
            RangeSet<Integer> rangeSet = TreeRangeSet.create();
            for (Region region : tempResultData.getRegionSet()) {
                rangeSet.add(Range.closed(region.getStartId(), region.getEndId()));
            }

            for (int i = 0; i < probeNum; i++) {
                if (rangeSet.contains(i)) {
                    voteNum[i]++;
                }
            }
        }

        Set<Integer> idSet = new HashSet<>();
        for (int i = 0; i < probeNum; i++) {
            if (enableDedugeInfo) System.out.println(i + "\t" + voteNum[i]);
            if (voteNum[i] >= Parameters.voteThreshold) {
                idSet.add(i);
            }
        }

        List<Integer> idList = new ArrayList<>(idSet);
        Collections.sort(idList);

        int tempStart = idList.get(0), tempId = idList.get(0), id;
        for (int i = 1; i < idList.size(); i++) {
            id = idList.get(i);
            if (id - tempId > 1) {
                Region tempRegion = new Region();
                tempRegion.setStartId(tempStart);
                tempRegion.setEndId(tempId);
                resultData.getRegionSet().add(tempRegion);

                tempStart = id;
            }
            tempId = id;
        }
        if (tempId == idList.get(idList.size() - 1)) {
            Region tempRegion = new Region();
            tempRegion.setStartId(tempStart);
            tempRegion.setEndId(tempId);
            resultData.getRegionSet().add(tempRegion);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("\n\nthe fianl regions: ");
        for (Region region : resultData.getRegionSet()) {
            int start = region.getStartId();
            int end = region.getEndId();
            m_log.info(String.format("Region [%d : %d : %d] Length = %d", start, (start + end) >> 1, end, region.getLength()));
            sb.append(String.format("[%d, %d], ",start,end));
        }
        m_log.info(sb.substring(0, sb.length() - 2));
    }

    static class Parameters {
        static int sampleSize = 6;
        static int sampleFrequency = 100;
        static int voteThreshold = 50;
    }
}

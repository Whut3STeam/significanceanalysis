package edu.whut.significance.methods;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.util.BioToolbox;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.StatUtils;

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
        sampleNum = rawData.getDataMatrix().getRowDimension();
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

        //vote(resultDataList);
        vote2(resultDataList);
    }

    public void sample(List<RawData> rawDataList) {
        ArrayList<Integer> randomList = new ArrayList<>();
        for (int i = 0; i < sampleNum; i++) {
            randomList.add(i);
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

        Set<Integer> idSet = new TreeSet<>();
        for (int i = 0; i < probeNum; i++) {
            if (enableDedugeInfo) System.out.println(i + "\t" + voteNum[i]);
            if (voteNum[i] >= Parameters.voteThreshold) {
                idSet.add(i);
            }
        }

        List<Integer> idList = new ArrayList<>(idSet);

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

    public void vote2(List<ResultData> resultDataList) {
        double[] voteNum = new double[probeNum];
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

        resultData.setRegionSet(getRegions(voteNum));

        StringBuilder sb = new StringBuilder();
        sb.append("\n\t\tthe fianl regions: ");
        for (Region region : resultData.getRegionSet()) {
            int start = region.getStartId();
            int end = region.getEndId();
            m_log.info(String.format("\n\t\t RF - Region [%d : %d : %d] Length = %d",
                    start, (start + end) >> 1, end, region.getLength()));
            sb.append(String.format("[%d, %d], ",start,end));
        }
        m_log.info(sb.substring(0, sb.length() - 2));

    }

    private Set<Region> getRegions(double[] voteNum){
        double[] data = BioToolbox.GaussianBlur(voteNum,5,1);

        double[] diff = new double[data.length];
        for (int i = 1; i < data.length; i++){
            diff[i] = data[i] - data[i - 1];
        }

        double max = StatUtils.max(diff);
        double min = StatUtils.min(diff);

        int k1 = 0, k2 = 0;
        for (int i = 0; i < diff.length; i++) {
            if (diff[i] == max) k1 = i;
            if (diff[i] == min) k2 = i;
        }

        //if (enableDedugeInfo){
            for (int i = -10; i < 11; i++){
                m_log.info(String.format("left = [%d:%.2f], right = [%d:%.2f]",
                        k1 + i,data[k1 + i], k2 + i,data[k2 + i]));
            }
        //}
        int posLeft = k1 < k2 ? k1:k2;
        int posRight = k1 > k2 ? k1:k2;

        Set<Region> result = new TreeSet<>();
        double thresh = Math.min(data[k1],data[k2]);
        for (int i = posLeft; i <= posRight; i++){
            if (data[i] < thresh){ //应该分成不同的区域, 通过寻找上升下降沿最陡峭之处
                return null;
            }
        }
        result.add(new Region(posLeft,posRight));
        return result;//如果有多个区域还有问题的
    }


    static class Parameters {
        static int sampleSize = 3;
        static int sampleFrequency = 100;
        static double voteThreshold = sampleFrequency / 2;
    }
}

package edu.whut.significance.methods;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/24.
 */
public class SAICEV extends AbstractSig {
    public ResultData resultData;
    public RealMatrix rawMatrix;
    public int rowNum;
    public int colNum;
    private Logger m_log;
    private boolean enableDedugeInfo = false;

    public void preprocess(RawData rawData) {
        rawMatrix = rawData.getDataMatrix().transpose();
        rowNum = rawMatrix.getRowDimension();
        colNum = rawMatrix.getColumnDimension();
        m_log = Logger.getLogger("significanceAnalysis");

        m_log.info(String.format("Current Row = %d, Col = %d", rowNum, colNum));
    }

    public void process(ResultData resultData) {
        RealMatrix[] rawMatrixs = classify(rawMatrix);
        RealMatrix ampRawMatrix = rawMatrixs[0];
        RealMatrix delRawMatrix = rawMatrixs[1];

        ResultData ampResultData = new ResultData();
        ResultData delResultData = new ResultData();

        Permute ampPermute = new Permute(ampRawMatrix, ampResultData);
        m_log.info(String.format("Staring to process amp : --------------->"));
        ampPermute.processing();

        Permute delPermute = new Permute(delRawMatrix, delResultData);
        m_log.info(String.format("Staring to process del : --------------->"));
        delPermute.processing();

        mergeResult(resultData, ampResultData, delResultData);
    }

    //��������ȱʧ������������ϲ�
    public void mergeResult(ResultData resultData, ResultData ampResultData, ResultData delResultData) {
        Set<Integer> idSet = new HashSet<>();
        for (Region ampRegion : ampResultData.getRegionSet()) {
            for (int i = ampRegion.getStartId(); i < ampRegion.getEndId() + 1; i++) {
                idSet.add(i);
            }
        }
        for (Region delRegion : delResultData.getRegionSet()) {
            for (int i = delRegion.getStartId(); i < delRegion.getEndId() + 1; i++) {
                idSet.add(i);
            }
        }
        List<Integer> idList = new ArrayList<>(idSet);
        Collections.sort(idList);

        if (idList.size() == 0) return;

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

        //if (enableDedugeInfo){
            StringBuilder sb = new StringBuilder();
            sb.append("the regions of one sample set: ");
            for (Region region : resultData.getRegionSet()) {
                int start = region.getStartId();
                int end = region.getEndId();
                m_log.info(String.format("Region [%d : %d : %d] Length = %d", start, (start + end) >> 1, end, region.getLength()));
                sb.append(String.format("[%d, %d], ",start,end));
            }
            m_log.info(sb.substring(0, sb.length() - 2));
        //}

    }

    public RealMatrix[] classify(RealMatrix rawMatrix) {
        GlobalParameters globalParameters = new GlobalParameters();
        double ampThreshold = globalParameters.getAmpThreshold();
        double delThreshold = globalParameters.getDelThreshold();

        RealMatrix ampRawMatrix = rawMatrix.copy();
//        for (int i = 0; i < rowNum; i++) {
//            for (int j = 1; j < colNum; j++) {
//                if (ampRawMatrix.getEntry(i, j) <= ampThreshold)
//                    ampRawMatrix.setEntry(i, j, 0);
//            }
//        }
        ampRawMatrix.walkInOptimizedOrder(new RealMatrixChangingVisitor() {

            @Override
            public void start(int rows, int columns, int startRow, int endRow, int startColumn, int endColumn) {

            }

            @Override
            public double end() {
                return 0;
            }

            @Override
            public double visit(int row, int column, double value) {
                return (value <= ampThreshold) ? 0 : value;
            }
        });

        RealMatrix delRawMatrix = rawMatrix.copy();
        for (int i = 0; i < rowNum; i++) {
            for (int j = 0; j < colNum; j++) {
                if (delRawMatrix.getEntry(i, j) >= delThreshold)
                    delRawMatrix.setEntry(i, j, 0);
            }
        }

        RealMatrix[] rawMatrixs = new RealMatrix[]{ampRawMatrix, delRawMatrix};
        return rawMatrixs;
    }

    static class Parameters {
        static double pccThreshold = 0.9;
        static int permuteNum = 50;
        static double sigValueThreshold = 0.0476;
        static int minCNALength = 6;
    }

    class Permute {
        //private RealMatrix oneRawMatrix;
        private RealMatrix oneRawDataMatrix;
        private ResultData oneResultData;
        private List<IdRegion> idRegionSet = new ArrayList<>();
        private List<CNARegion> CNARegionSet = new ArrayList<>();
        private Set<Integer> lengthSet = new TreeSet<>();
        private int permuteProbeSize;

        public Permute(RealMatrix oneRawMatrix, ResultData oneResultData) {
            //this.oneRawMatrix = oneRawMatrix;
            this.oneRawDataMatrix = oneRawMatrix;
            this.oneResultData = oneResultData;
        }

        public void processing() {
            getCNAs();
            if (enableDedugeInfo){
                m_log.info(String.format("candidate Region = %d", idRegionSet.size()));
                m_log.info(String.format("CNA region = %d", CNARegionSet.size()));
            }

            if (CNARegionSet.size() >= 2)
                permuteDetection();
            else
                m_log.info("There is not enough CNA units to be permuted");
            //m_log.info("There is not enough CNA units to be permuted");
        }

        //��ó�ʼCNA��Ԫ�ı��
        public void getCNAs() {
            Set<Integer> candidateIdSet = new TreeSet<>();
            ArrayList<IdRegion> tempIdRegionSet = new ArrayList<>();

            getCandidateId(candidateIdSet);

            if (candidateIdSet.size() > 0) { //Ӧ�ö�ȫ�ֲ�����
                getTempIdRegion(candidateIdSet, tempIdRegionSet);
                getIdRegion(tempIdRegionSet);
                regionMerge();
            } else {
                m_log.info("There are no enough CNA probes");
            }
        }

        //��ò�Ϊ0��Id�ŵļ���
        public void getCandidateId(Set<Integer> candidateIdSet) {
            for (int i = 0; i < rowNum; i++) {
                for (int j = 0; j < colNum; j++) {
                    if (oneRawDataMatrix.getEntry(i, j) != 0) {
                        candidateIdSet.add(i);
                        break;
                    }
                }
            }
        }

        //��ȡ���Ȳ�С��6�ĵ�Ԫ�ı�ż���
        public void getTempIdRegion(Set<Integer> candidateIdSet, List<IdRegion> tempIdRegionSet) {

            int count = candidateIdSet.size();
            m_log.info(String.format("the number of candidate probes = %d", count));

            int previous = -1;
            List<Integer> IDArray = new LinkedList<>();
            for (int current : candidateIdSet) {
                if ((previous < 0) || (previous >= 0 && current - previous == 1)) {
                    IDArray.add(current);
                    previous = current;
                } else {
                    //�Ͽ���
                    int begin = IDArray.get(0);
                    int end = IDArray.get(IDArray.size() - 1);
                    tempIdRegionSet.add(new IdRegion(begin, end));
                    IDArray.clear();

                    //�µĿ�ʼ
                    IDArray.add(current);
                    previous = -1;
                }
            }

            int begin = IDArray.get(0);
            int end = IDArray.get(IDArray.size() - 1);
            tempIdRegionSet.add(new IdRegion(begin, end));


        }

        //��û���Ƥ��ɭ���ϵ���ֶκ��CNA��Ԫ
        public void getIdRegion(List<IdRegion> tempIdRegionSet) {
            int tempStart = 0, tempEnd = 0, tempId = 0;
            double peaCorCoe = 0.0;
            PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();
            for (IdRegion tempIdRegion : tempIdRegionSet) {
                tempStart = tempIdRegion.getStart();
                tempEnd = tempIdRegion.getEnd();
                for (tempId = tempStart; tempId < tempEnd; tempId++) {
                    peaCorCoe = pearsonsCorrelation.correlation(oneRawDataMatrix.getRow(tempId),
                            oneRawDataMatrix.getRow(tempId + 1));
                    if (peaCorCoe < Parameters.pccThreshold) {
                        if ((tempId - tempStart + 1) >= Parameters.minCNALength) {
                            IdRegion validIdRegion = new IdRegion(tempStart, tempId);
                            idRegionSet.add(validIdRegion);
                        }
                        tempStart = tempId + 1;
                    }
                }
                //���һ��̽�뵥������
                if ((tempId == tempEnd) && ((tempId - tempStart + 1) >= Parameters.minCNALength)) {
                    IdRegion validIdRegion = new IdRegion(tempStart, tempId);
                    idRegionSet.add(validIdRegion);
                }
            }
        }


        //��֯CNA��Ԫ�������յ�CNA��Ԫ�ںϳ�CNARegionSet�����һ��lengthSet
        public void regionMerge() {
            int cnaId = 0;
            int cnaLength = 0;

            if (idRegionSet.size() > 0) {

                for (IdRegion idRegion : idRegionSet) {
                    cnaLength = idRegion.getEnd() - idRegion.getStart() + 1;

                    CNARegion tempCNARegion = new CNARegion(cnaId, idRegion);
                    CNARegionSet.add(tempCNARegion);
                    lengthSet.add(cnaLength);
                    permuteProbeSize += cnaLength;
                    cnaId++;
                }
            }
        }

        //�����Լ��
        public void permuteDetection() {

            double[][] maxUScore = new double[Parameters.permuteNum][lengthSet.size()];

            calUScore();

            boolean loopFlag = true;
            List<CNARegion> permuteCNARegions = new LinkedList<>(CNARegionSet);
            //int i = 0;
            while (loopFlag == true) {
                permute(permuteCNARegions, maxUScore);
                loopFlag = SCAExclude(maxUScore, permuteCNARegions);
                //System.out.println(i++);
            }
            calPScore(maxUScore);
            getResultData();
        }

        //��������CNA��Ԫ��Uֵ
        public void calUScore() {
            int start, end;
            int cnaLength;
            double tempUScore;

            for (CNARegion region : CNARegionSet) {
                start = region.getIdRegion().getStart();
                end = region.getIdRegion().getEnd();
                cnaLength = region.getLength();
                tempUScore = 0.0;
                for (int i = start; i <= end; i++) {
                    tempUScore += StatUtils.sum(oneRawDataMatrix.getRow(i));
                }
                tempUScore = Math.abs(tempUScore / (cnaLength * colNum));
                region.setuValue(tempUScore);
            }
        }

        public double getSum(double[] nums) {
            double sum = 0;
            for (int i = 0; i < nums.length; i++) {
                sum += nums[i];
            }
            return sum;
        }

        //CNA��Ԫ�������������Uֵ����
        public void permute(List<CNARegion> permuteCNARegions, double[][] maxUScore) {
            List<Integer> idArray = new LinkedList<>();
            for (int i = 0; i < permuteCNARegions.size(); i++) {
                idArray.add(i);
            }

            for (int i = 0; i < Parameters.permuteNum; i++) {
                //һ�ν�����ʼ
                RealMatrix randomPermuteMatrix = new BlockRealMatrix(permuteProbeSize, colNum);
                for (int j = 0; j < colNum; j++) {
                    Collections.shuffle(idArray);
                    List<Integer> selectedIDs = idArray.subList(0, idArray.size() - 3);

                    double[] newCol = new double[permuteProbeSize];
                    double[] orignalCol = oneRawDataMatrix.getColumn(j);

                    int k = 0;
                    for (int id : selectedIDs) {
                        CNARegion region = permuteCNARegions.get(id);
                        int begin = region.getIdRegion().getStart();
                        int len = region.getLength();
                        System.arraycopy(orignalCol, begin, newCol, k, len);
                        k += len;
                    }
                    randomPermuteMatrix.setColumn(j, newCol);
                }

                findMaxUScore(randomPermuteMatrix, maxUScore[i]);
            }
        }

        //�����Uֵ���ڵ� i ��ʵ��
        public void findMaxUScore(RealMatrix randomPermuteMatrix, double[] maxUScoreAti) {
            double tempUScore;
            int probeNum = randomPermuteMatrix.getRowDimension();
            double[] probeSum = new double[probeNum + 1];
            //����ÿ��̽��ĺ�
            probeSum[0] = 0;
            for (int i = 1; i < probeNum + 1; i++) {
                probeSum[i] = probeSum[i - 1] + StatUtils.sum(randomPermuteMatrix.getRow(i - 1));
            }

            List<Integer> uniqueLengthSet = new ArrayList<>(lengthSet);
            for (int length : uniqueLengthSet) {
                int index = uniqueLengthSet.indexOf(length);
                maxUScoreAti[index] = 0.0;

                int endPos = probeNum - length + 1;
                int front, back;

                for (front = 0, back = front + length; back < endPos; front++, back++) {
                    double regionSum = probeSum[back] - probeSum[front];
                    tempUScore = Math.abs(regionSum) / (colNum * length);
                    if (tempUScore > maxUScoreAti[index]) {
                        maxUScoreAti[index] = tempUScore;
                    }
                }
            }
        }

        //�ų�SCAs������permuteCNARegions
        public boolean SCAExclude(double[][] maxUScore, List<CNARegion> permuteCNARegions) {
            double sigValue;
            int index;
            boolean flag = false;
            int maxLen = Collections.max(lengthSet);//To make sure there are enough probes
            int permuteNum = Parameters.permuteNum;
            double sigValueThreshold = Parameters.sigValueThreshold;

            int permuteRegionCount = permuteCNARegions.size();
            ArrayList<Integer> tempLengthSet = new ArrayList<>(lengthSet);

            //���ݳ�������CNA��Ԫ��Uֵ
            Iterator<CNARegion> itr = permuteCNARegions.iterator();
            while (itr.hasNext()) {
                CNARegion region = itr.next();
                sigValue = 0.0;
                index = tempLengthSet.indexOf(region.getLength());

                //�ۻ�����Uֵ�����Uֵ�ø���
                for (int j = 0; j < permuteNum; j++) {
                    if (maxUScore[j][index] > region.getuValue()) {
                        sigValue++;
                    }
                }
                sigValue /= (permuteNum + 1);

                if (sigValue < sigValueThreshold) {
                    flag = true;
                    if (permuteRegionCount - region.getLength() > maxLen)
                        permuteProbeSize -= region.getLength();
                    //permuteCNARegions.remove(i);
                    itr.remove();
                    permuteRegionCount = permuteCNARegions.size();
                }
            }
            return flag;
        }

        //�ٴμ�������CNA��Ԫ��Pֵ
        public void calPScore(double[][] maxUScore) {
            double pValue;
            int index;
            ArrayList<Integer> uniqueLengthSet = new ArrayList<>(lengthSet);

            if (enableDedugeInfo){
                m_log.info("<<<< all CNARegion >>>>");
            }

            for (CNARegion region : CNARegionSet) {
                pValue = 0.0;
                index = uniqueLengthSet.indexOf(region.getLength());

                for (int j = 0; j < Parameters.permuteNum; j++) {
                    if (maxUScore[j][index] > region.getuValue()) {
                        pValue++;
                    }
                }

                pValue /= (Parameters.permuteNum + 1);
                region.setpValue(pValue);
                if (enableDedugeInfo) m_log.info("\t>>>>" + region.toString());
            }
        }

        //get ResultDatas and update the UScore,PScore and SCATag of the CNARegionSet
        //��ý��������CNARegionSet��Uֵ��Pֵ��SCATag
        public void getResultData() {
            double pScoreThreshold = Parameters
                    .sigValueThreshold;
            //Set<Region> tempResultRegions = new HashSet<>();

            for (CNARegion region : CNARegionSet) {
                region.setSCATag(0);
                if (region.getpValue() < pScoreThreshold) {
                    region.setSCATag(1);
                    oneResultData.addRegion(region.getIdRegion().getStart(), region.getIdRegion().getEnd());
                    if (enableDedugeInfo) m_log.info(region.toString());
                }
            }
        }
    }

    class IdRegion {
        private int start;
        private int end;

        public IdRegion() {
        }

        public IdRegion(int start, int end) {
            this.start = start;
            this.end = end;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public int getEnd() {
            return end;
        }

        public void setEnd(int end) {
            this.end = end;
        }

        @Override
        public String toString() {
            return "{" +
                    "start=" + start +
                    ", end=" + end +
                    '}';
        }
    }

    class CNARegion {
        private int cnaId;
        private IdRegion idRegion;
        //private int length;
        private double uValue;
        private double pValue;
        private int SCATag;

        public CNARegion(int id, IdRegion region) {
            cnaId = id;
            idRegion = region;
        }

        public int getCnaId() {
            return cnaId;
        }

        public void setCnaId(int cnaId) {
            this.cnaId = cnaId;
        }

        public IdRegion getIdRegion() {
            return idRegion;
        }

        public void setIdRegion(IdRegion idRegion) {
            this.idRegion = idRegion;
        }

        public int getLength() {
            return idRegion.getEnd() - idRegion.getStart() + 1;
        }

        public double getuValue() {
            return uValue;
        }

        public void setuValue(double uValue) {
            this.uValue = uValue;
        }

        public double getpValue() {
            return pValue;
        }

        public void setpValue(double pValue) {
            this.pValue = pValue;
        }

        public int getSCATag() {
            return SCATag;
        }

        public void setSCATag(int SCATag) {
            this.SCATag = SCATag;
        }

        @Override
        public String toString() {
            return String.format("CNARegion{ %s, uValue = %.4f, pValue = %.4e }", idRegion, uValue, pValue);
        }
    }
}
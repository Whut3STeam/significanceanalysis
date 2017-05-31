package edu.whut.significance.methods;

import com.sun.xml.internal.bind.v2.model.core.ID;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.util.BioLogger;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathUtils;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/24.
 */
public class SAIC extends AbstractSig {
    private Logger m_log;

    public ResultData resultData;
    public RealMatrix rawMatrix;
    public int rowNum;
    public int colNum;

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
        m_log.info(String.format("Staring to process amp :"));
        ampPermute.processing();

        Permute delPermute = new Permute(delRawMatrix, delResultData);
        m_log.info(String.format("Staring to process del :"));
        delPermute.processing();

        mergeResult(resultData, ampResultData, delResultData);
    }

    //将扩增和缺失的显著性区域合并
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

        for (Region region : resultData.getRegionSet()) {
            int start = region.getStartId();
            int end = region.getEndId();
            m_log.info(String.format("Region [%d : %d : %d] Length = %d", start, (start + end) >> 1, end, region.getLength()));
        }
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

//        RealMatrix[] rawMatrixs = new RealMatrix[2];
//        rawMatrixs[0] = ampRawMatrix;
//        rawMatrixs[1] = delRawMatrix;
        RealMatrix[] rawMatrixs = new RealMatrix[]{ampRawMatrix, delRawMatrix};
        return rawMatrixs;
    }

    class Permute {
        //private RealMatrix oneRawMatrix;
        private RealMatrix oneRawDataMatrix;
        private ResultData oneResultData;
        private List<IdRegion> idRegionSet = new ArrayList<>();
        private List<CNARegion> CNARegionSet = new ArrayList<>();
        private Set<Integer> lengthSet = new TreeSet<>();
        //private List<Integer> lengthSet_copy;
        //private List<Integer> uniqueLengthSet;
        //private List<Double> UScore = new ArrayList<>();
        //private List<Double> UScore_copy;
        //private List<Double> PScore = new ArrayList<>();
        private int permuteProbeSize;

        public Permute(RealMatrix oneRawMatrix, ResultData oneResultData) {
            //this.oneRawMatrix = oneRawMatrix;
            this.oneRawDataMatrix = oneRawMatrix;
            this.oneResultData = oneResultData;
        }

        public void processing() {
            getCNAs();
            m_log.info(String.format("candidate Region = %d", idRegionSet.size()));
            m_log.info(String.format("CNA region = %d", CNARegionSet.size()));

            if (CNARegionSet.size() >= 2)
                permuteDetection();
            else
                m_log.info("There is not enough CNA units to be permuted");
            //m_log.info("There is not enough CNA units to be permuted");
        }

        //获得初始CNA单元的编号
        public void getCNAs() {
            Set<Integer> candidateIdSet = new TreeSet<>();
            ArrayList<IdRegion> tempIdRegionSet = new ArrayList<>();

            getCandidateId(candidateIdSet);

            if (candidateIdSet.size() > 0) { //应该读全局参数？
                getTempIdRegion(candidateIdSet, tempIdRegionSet);
                getIdRegion(tempIdRegionSet);
                regionMerge();
            } else {
                m_log.info("There are no enough CNA probes");
            }
        }

        //获得不为0的Id号的集合
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

        //获取长度不小于6的单元的编号集合
        public void getTempIdRegion(Set<Integer> candidateIdSet, List<IdRegion> tempIdRegionSet) {
//            int tempStart = candidateIdSet.get(0), tempEnd = candidateIdSet.get(0), tempId;
//            for (int i = 1; i < candidateIdSet.size(); i++) {
//                tempId = candidateIdSet.get(i);
//                if (tempId - tempEnd != 1) {
//                    if ((tempEnd - tempStart + 1) >= Parameters.minCNALength) {
//                        IdRegion tempIdRegion = new IdRegion();
//                        tempIdRegion.setStart(tempStart);
//                        tempIdRegion.setEnd(tempEnd);
//                        tempIdRegionSet.add(tempIdRegion);
//                    }
//                    tempStart = tempId;
//                    tempEnd = tempId;
//                } else {
//                    tempEnd = tempId;
//                }
//                //最后一个探针单独处理
//                if ((i == candidateIdSet.size() - 1) && ((tempId - tempStart + 1) > Parameters.minCNALength)) {
//                    IdRegion tempIdRegion = new IdRegion();
//                    tempIdRegion.setStart(tempStart);
//                    tempIdRegion.setEnd(tempId);
//                    tempIdRegionSet.add(tempIdRegion);
//                }
//            }
            int count = candidateIdSet.size();
            m_log.info(String.format("the number of candidate probes = %d", count));

            int previous = -1;
            List<Integer> IDArray = new LinkedList<>();
            for (int current : candidateIdSet) {
                if ((previous < 0) || (previous >= 0 && current - previous == 1)) {
                    IDArray.add(current);
                    previous = current;
                } else {
                    //断开了
                    int begin = IDArray.get(0);
                    int end = IDArray.get(IDArray.size() - 1);
                    tempIdRegionSet.add(new IdRegion(begin, end));
                    IDArray.clear();

                    //新的开始
                    IDArray.add(current);
                    previous = -1;
                }
            }

            int begin = IDArray.get(0);
            int end = IDArray.get(IDArray.size() - 1);
            tempIdRegionSet.add(new IdRegion(begin, end));


        }

        //获得基于皮尔森相关系数分段后的CNA单元
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
                //最后一个探针单独处理
                if ((tempId == tempEnd) && ((tempId - tempStart + 1) >= Parameters.minCNALength)) {
                    IdRegion validIdRegion = new IdRegion(tempStart, tempId);
                    idRegionSet.add(validIdRegion);
                }
            }
        }

//        //calculate the Pearson Correlation Coefficient
//        public double getPeaCorCoe(int id1,int id2,ArrayList<RawData>  candidateRawDatas){
//            double peaCorCoe;
//            double mean1,mean2,std1,std2,cov;
//
//            mean1=0.0;mean2=0.0;
//            for(int i=0;i<sampleNum;i++){
//                mean1+=candidateRawDatas.get(id1).getData()[i];
//                mean2+=candidateRawDatas.get(id2).getData()[i];
//            }
//            mean1/=sampleNum;
//            mean2/=sampleNum;
//
//            std1=0.0;std2=0.0;cov=0.0;
//            for (int i=0;i<sampleNum;i++){
//                std1+=(candidateRawDatas.get(id1).getData()[i]-mean1)*(candidateRawDatas.get(id1).getData()[i]-mean1);
//                std2+=(candidateRawDatas.get(id2).getData()[i]-mean2)*(candidateRawDatas.get(id2).getData()[i]-mean2);
//                cov+=(candidateRawDatas.get(id1).getData()[i]-mean1)*(candidateRawDatas.get(id2).getData()[i]-mean2);
//            }
//            std1/=sampleNum;
//            std2/=sampleNum;
//
//            if (std1==0||std2==0)
//                return 0;
//
//            peaCorCoe=cov/(Math.sqrt(std1)*Math.sqrt(std2)*(sampleNum-1));
//
//            return peaCorCoe;
//        }

        //组织CNA单元并将最终的CNA单元融合成CNARegionSet，并且获得lengthSet
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
            //lengthSet_copy = new ArrayList<>(lengthSet);
        }

        //显著性检测
        public void permuteDetection() {

//            //获得CNA单元的唯一长度
//            Set uniqueLength = new HashSet(lengthSet);
//            uniqueLengthSet = new ArrayList(uniqueLength);
//            Collections.sort(uniqueLengthSet);
//            //ArrayList uniqueLengthSet_copy=new ArrayList(uniqueLengthSet);

            double[][] maxUScore = new double[Parameters.permuteNum][lengthSet.size()];

            calUScore();

            boolean loopFlag = true;
            List<CNARegion> permuteCNARegions = new LinkedList<>(CNARegionSet);
            while (loopFlag == true) {
                permute(permuteCNARegions, maxUScore);
                loopFlag = SCAExclude(maxUScore, permuteCNARegions);
            }
            calPScore(maxUScore);
            getResultData();
        }

        //计算所有CNA单元的U值
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
                    //tempUScore += getSum(oneRawDataMatrix.getRow(i));
                    tempUScore += StatUtils.sum(oneRawDataMatrix.getRow(i));
                }
                tempUScore = Math.abs(tempUScore / (cnaLength * colNum));

                //UScore.add(tempUScore);
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

        //CNA单元交换，生成最大U值矩阵
        public void permute(List<CNARegion> permuteCNARegions, double[][] maxUScore) {
            List<Integer> idArray = new LinkedList<>();
            for (int i = 0; i < permuteCNARegions.size(); i++) {
                idArray.add(i);
            }

            for (int i = 0; i < Parameters.permuteNum; i++) {
                //一次交换开始
                RealMatrix randomPermuteMatrix = new BlockRealMatrix(permuteProbeSize, colNum);
                for (int j = 0; j < colNum; j++) {
                    Collections.shuffle(idArray);
                    double[] newCol = new double[permuteProbeSize];
                    double[] orignalCol = oneRawDataMatrix.getColumn(j);

                    int k = 0;
                    for (int id : idArray) {
                        CNARegion region = permuteCNARegions.get(id);
                        int begin = region.getIdRegion().getStart();
                        int len = region.getLength();
                        System.arraycopy(orignalCol, begin, newCol, k, len);
                        k += len;
                    }
                    randomPermuteMatrix.setColumn(j, newCol);
                }

//                int maxIndex = permuteCNARegions.size();
//                int index, len, cnaStart;
//
//                //Generate sampleNum*permuteMatrix.length random numbers
//                //生成 样本数*交换探针数 个随机数
//                Random rn = new Random();
//                for (int j = 0; j < colNum - 1; j++) {
//                    maxIndex = permuteCNARegions.size();
//                    boolean[] bool = new boolean[maxIndex];
//                    int randomInt = 0;
//                    int pStart = 0;
//
//                    //在每个样本之间进行交换
//                    for (int k = 0; k < permuteCNARegions.size(); k++) {
//
//                        //生成不重复的随机数
//                        do {
//                            randomInt = rn.nextInt(maxIndex);
//                        } while (bool[randomInt]);
//
//                        bool[randomInt] = true;
//                        index = randomInt;
//                        len = permuteCNARegions.get(index).getLength();
//                        cnaStart = permuteCNARegions.get(index).getIdRegion().getStart();
//
//                        //根据序号和原有矩阵更新交换矩阵
////                        randomPermuteMatrix.setSubMatrix(
////                                oneRawDataMatrix.getSubMatrix(cnaStart-pStart,cnaStart+len-pStart-1,j,j).getData(),pStart,j);
//                        randomPermuteMatrix.setSubMatrix(
//                                oneRawDataMatrix.getSubMatrix(cnaStart, cnaStart + len - 1, j, j).getData(), pStart, j);
//                        /*for(int m=pStart;m<(pStart+len);m++){
//                            randomPermuteMatrix[m][j]=candidateRawDatas.get(cnaStart+m-pStart).getData()[j];
//                        }*/
//                        pStart += len;
//                    }
//                }
                findMaxUScore(randomPermuteMatrix, maxUScore[i]);
            }
        }

        //找最大U值，在第 i 次实验
        public void findMaxUScore(RealMatrix randomPermuteMatrix, double[] maxUScoreAti) {
            double tempUScore;
            int probeNum = randomPermuteMatrix.getRowDimension();
            double[] probeSum = new double[probeNum + 1];
            //计算每个探针的和
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
//            double tempUScore;
//            int probeNum=randomPermuteMatrix.getRowDimension();
//            double[] probeSum=new double[probeNum];
//            int i,j,k;
//            int len;
//            int endPos;
//
//            //计算每个探针的和
//            for(i=0;i<probeNum;i++){
//                probeSum[i]=0;
//                probeSum[i]+=getSum(randomPermuteMatrix.getRow(i));
//                /*for(j=0;j<colNum-1;j++){
//                    probeSum[i]+=randomPermuteMatrix[i][j];
//                }*/
//                probeSum[i]=Math.abs(probeSum[i]);
//            }
//
//            //计算第一个唯一长度的U值
//            List<Integer> uniqueLengthSet = new ArrayList<>(lengthSet);
//            len=uniqueLengthSet.get(0);
//            tempUScore=0.0;
//            maxUScoreAti[0]=0.0;
//            endPos=probeNum-len;
//            double[] regionSum=new double[endPos];
//            for (j=0;j<endPos;j++){
//                for (k=j;k<j+len;k++){
//                    tempUScore+=probeSum[k];
//                }
//                regionSum[j]=Math.abs(tempUScore);
//                tempUScore=Math.abs(tempUScore/((colNum)*len));
//                if(tempUScore>maxUScoreAti[0]){
//                    maxUScoreAti[0]=tempUScore;
//                }
//            }
//
//            //计算余下唯一长度的U值
//            for(i=1;i<uniqueLengthSet.size();i++){
//                len=uniqueLengthSet.get(i);
//                endPos=probeNum-len;
//                maxUScoreAti[i]=0.0;
//                for (j=0;j<endPos;j++){
//                    tempUScore=regionSum[j];
//                    for (k=j+uniqueLengthSet.get(i-1);k<j+len;k++){
//                        tempUScore+=probeSum[k];
//                    }
//                    regionSum[j]=tempUScore;
//                    tempUScore=tempUScore/((colNum)*len);
//                    if(tempUScore>maxUScoreAti[i]){
//                        maxUScoreAti[i]=tempUScore;
//                    }
//                }
//            }

        }

        //排除SCAs并更新permuteCNARegions
        public boolean SCAExclude(double[][] maxUScore, List<CNARegion> permuteCNARegions) {
            double sigValue;
            int index;
            boolean flag = false;
            int maxLen = Collections.max(lengthSet);//To make sure there are enough probes
            int permuteNum = Parameters.permuteNum;
            double sigValueThreshold = Parameters.sigValueThreshold;

            int permuteRegionCount = permuteCNARegions.size();
            ArrayList<Integer> tempLengthSet = new ArrayList<>(lengthSet);

            //根据长度搜索CNA单元的U值
            Iterator<CNARegion> itr = permuteCNARegions.iterator();
            while (itr.hasNext()) {
                CNARegion region = itr.next();
                sigValue = 0.0;
                index = tempLengthSet.indexOf(region.getLength());

                //累积大于U值得最大U值得个数
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

//            //根据长度搜索CNA单元的U值
//            for (int i = 0; i < permuteSize; i++) {
//                sigValue = 0.0;
//                index = uniqueLengthSet.indexOf(lengthSet_copy.get(i));

//                //累积大于U值得最大U值得个数
//                for (int j = 0; j < permuteNum; j++) {
//                    if (maxUScore[j][index] > UScore.get(i)) {
//                        sigValue++;
//                    }
//                }
//                sigValue /= (permuteNum + 1);
//
//                if (sigValue < sigValueThreshold) {
//                    flag = true;
//                    if (permuteSize - tempLengthSet.get(i) > maxLen)
//                        permuteProbeSize -= tempLengthSet.get(i);
//                    permuteCNARegions.remove(i);
//                    UScore_copy.remove(i);
//                    lengthSet_copy.remove(i);
//                    permuteSize = permuteCNARegions.size();
//                }
//            }
//            return flag;
        }

        //再次计算所有CNA单元的P值
        public void calPScore(double[][] maxUScore) {
            double pValue;
            int index;
            ArrayList<Integer> uniqueLengthSet = new ArrayList<>(lengthSet);

            m_log.info("<<<< all CNARegion >>>>");
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
                m_log.info("\t>>>>" + region.toString());
            }

//
//
//            for (int i = 0; i < lengthSet.size(); i++) {
//                pValue = 0.0;
//                index = uniqueLengthSet.indexOf(it.next());
//                for (int j = 0; j < Parameters.permuteNum; j++) {
//                    if (maxUScore[j][index] > UScore.get(i)) {
//                        pValue++;
//                    }
//                }
//                pValue /= (Parameters.permuteNum + 1);
//                PScore.add(pValue);
//            }
        }

        //get ResultDatas and update the UScore,PScore and SCATag of the CNARegionSet
        //获得结果并更新CNARegionSet的U值，P值和SCATag
        public void getResultData() {
            double pScoreThreshold = Parameters
                    .sigValueThreshold;
            //Set<Region> tempResultRegions = new HashSet<>();

            for (CNARegion region : CNARegionSet) {
                region.setSCATag(0);
                if (region.getpValue() < pScoreThreshold) {
                    region.setSCATag(1);
                    oneResultData.addRegion(region.getIdRegion().getStart(), region.getIdRegion().getEnd());
                    m_log.info(region.toString());
                }
            }
        }

//            for (int i = 0; i < CNARegionSet.size(); i++) {
//                CNARegion tempCNARegion = new CNARegion();
//                tempCNARegion = CNARegionSet.get(i);
//                tempCNARegion.setuValue(UScore.get(i));
//                tempCNARegion.setpValue(PScore.get(i));
//                tempCNARegion.setSCATag(0);
//                if (PScore.get(i) < pScoreThreshold) {
//                    tempCNARegion.setSCATag(1);
//                    Region tempRegion = new Region();
//                    tempRegion.setStartId(CNARegionSet.get(i).getIdRegion().getStart());
//                    tempRegion.setEndId(CNARegionSet.get(i).getIdRegion().getEnd());
//                    //tempRegion.setLength(CNARegionSet.get(i).getLength());
//                    tempResultRegions.add(tempRegion);
//                }
//                CNARegionSet.set(i, tempCNARegion);
//            }
//            oneResultData.setRegionSet(tempResultRegions);
//        }
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

        public int getEnd() {
            return end;
        }

        public void setStart(int start) {
            this.start = start;
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

        public IdRegion getIdRegion() {
            return idRegion;
        }

        public int getLength() {
            return idRegion.getEnd() - idRegion.getStart() + 1;
        }

        public double getuValue() {
            return uValue;
        }

        public double getpValue() {
            return pValue;
        }

        public int getSCATag() {
            return SCATag;
        }

        public void setCnaId(int cnaId) {
            this.cnaId = cnaId;
        }

        public void setIdRegion(IdRegion idRegion) {
            this.idRegion = idRegion;
        }

        public void setuValue(double uValue) {
            this.uValue = uValue;
        }

        public void setpValue(double pValue) {
            this.pValue = pValue;
        }

        public void setSCATag(int SCATag) {
            this.SCATag = SCATag;
        }

        @Override
        public String toString() {
            return String.format("CNARegion{ %s, uValue = %.4f, pValue = %.4e }", idRegion, uValue, pValue);
        }
    }

    static class Parameters {
        static double pccThreshold = 0.9;
        static int permuteNum = 1000;
        static double sigValueThreshold = 0.0476;
        static int minCNALength = 6;
    }
}

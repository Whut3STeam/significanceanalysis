package edu.whut.significance.methods;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by SunMing on 2017/5/24.
 */
public class SAIC extends AbstractSig {
    ResultData resultData;
    RealMatrix rawMatrix;
    int rowNum;
    int colNum;

    public void preprocess(RawData rawData){
        rawMatrix=rawData.getDataMatrix().transpose();
        rowNum=rawMatrix.getRowDimension();
        colNum=rawMatrix.getRow(0).length;
    }

    public void process(ResultData resultData){
        RealMatrix[] rawMatrixs=classify(rawMatrix);
        RealMatrix ampRawMatrix=rawMatrixs[0];
        RealMatrix delRawMatrix=rawMatrixs[1];

        ResultData ampResultData=new ResultData();
        ResultData delResultData=new ResultData();

        Permute ampPermute=new Permute(ampRawMatrix,ampResultData);
        //m_log.info(conf.getChrName()+" amp:");
        ampPermute.processing();
        Permute delPermute=new Permute(delRawMatrix,delResultData);
        //m_log.info("\n");
        //m_log.info(conf.getChrName()+" del:");
        delPermute.processing();

        mergeResult(resultData,ampResultData,delResultData);
    }

    //将扩增和缺失的显著性区域合并
    public void mergeResult(ResultData resultData,ResultData ampResultData,ResultData delResultData){
        Set<Integer> idSet=new HashSet<>();
        for(Region ampRegion:ampResultData.getRegionSet()){
            for(int i=ampRegion.getStartId();i<ampRegion.getEndId()+1;i++){
                idSet.add(i);
            }
        }
        for(Region delRegion:delResultData.getRegionSet()){
            for(int i=delRegion.getStartId();i<delRegion.getEndId()+1;i++){
                idSet.add(i);
            }
        }
        List<Integer> idList=new ArrayList<>(idSet);
        Collections.sort(idList);

        int tempStart=idList.get(0),tempId=idList.get(0),id;
        for(int i=1;i<idList.size();i++){
            id=idList.get(i);
            if(id-tempId>1){
                Region tempRegion=new Region();
                tempRegion.setStartId(tempStart);
                tempRegion.setEndId(tempId);
                tempRegion.setLength(tempId-tempStart+1);
                resultData.getRegionSet().add(tempRegion);

                tempStart=id;
            }
            tempId=id;
        }
    }

    public RealMatrix[] classify(RealMatrix rawMatrix){
        GlobalParameters globalParameters=new GlobalParameters();
        double ampThreshold=globalParameters.getAmpThreshold();
        double delThreshold=globalParameters.getDelThreshold();

        RealMatrix ampRawMatrix=rawMatrix.copy();
        for(int i=0;i<rowNum;i++){
            for(int j=1;j<colNum;j++){
                if(ampRawMatrix.getEntry(i,j)<=ampThreshold)
                    ampRawMatrix.setEntry(i,j,0);
            }
        }

        RealMatrix delRawMatrix=rawMatrix.copy();
        for(int i=0;i<rowNum;i++){
            for(int j=1;j<colNum;j++){
                if(delRawMatrix.getEntry(i,j)>=delThreshold)
                    delRawMatrix.setEntry(i,j,0);
            }
        }

        RealMatrix[] rawMatrixs=new RealMatrix[2];
        rawMatrixs[0]=ampRawMatrix;
        rawMatrixs[1]=delRawMatrix;

        return rawMatrixs;
    }

    class Permute {
        private RealMatrix oneRawMatrix;
        private RealMatrix oneRawDataMatrix;
        private ResultData oneResultData;
        private List<IdRegion> idRegionSet=new ArrayList<>();
        private List<CNARegion> CNARegionSet=new ArrayList<>();
        private List<Integer> lengthSet=new ArrayList<>();
        private List<Integer> lengthSet_copy;
        private List<Integer> uniqueLengthSet;
        private List<Double> UScore=new ArrayList<>();
        private List<Double> UScore_copy;
        private List<Double> PScore=new ArrayList<>();
        private int permuteProbeSize;
        //private Logger m_log;

        public Permute(RealMatrix oneRawMatrix,ResultData oneResultData){
            this.oneRawMatrix=oneRawMatrix;
            this.oneRawDataMatrix=oneRawMatrix.getSubMatrix(0,rowNum-1,1,colNum-1);
            this.oneResultData=oneResultData;
            //m_log=Logger.getLogger("edu.vt.cbil.permute");
        }

        public void processing(){
            getCNAs();
            if(CNARegionSet.size()>=2)
                permuteDetection();
            else
                System.out.println("There is not enough CNA units to be permuted");
                //m_log.info("There is not enough CNA units to be permuted");
        }

        //获得初始CNA单元的编号
        public void getCNAs(){
            ArrayList<Integer> candidateIdSet=new ArrayList<>();
            ArrayList<IdRegion> tempIdRegionSet=new ArrayList<>();

            getCandidateId(candidateIdSet);
            if(candidateIdSet.size()>0){
                getTempIdRegion(candidateIdSet, tempIdRegionSet);
                getIdRegion(tempIdRegionSet);
                regionMerge();
            }
            else {
                System.out.println("There are no enough CNA probes");
            }
        }

        //获得不为0的Id号的集合
        public void getCandidateId(ArrayList<Integer> candidateIdSet){
            for(int i=0;i<rowNum;i++){
                for(int j=1;j<colNum;j++){
                    if(oneRawMatrix.getEntry(i,j)!=0){
                        candidateIdSet.add(i);
                        break;
                    }
                }
            }
        }

        //获取长度不小于6的单元的编号集合
        public void getTempIdRegion(ArrayList<Integer> candidateIdSet,ArrayList<IdRegion> tempIdRegionSet){
            int tempStart=candidateIdSet.get(0),tempEnd=candidateIdSet.get(0),tempId;
            for(int i=1;i<candidateIdSet.size();i++){
                tempId=candidateIdSet.get(i);
                if(tempId-tempEnd!=1){
                    if((tempEnd-tempStart+1)>=Parameters.minCNALength){
                        IdRegion tempIdRegion=new IdRegion();
                        tempIdRegion.setStart(tempStart);
                        tempIdRegion.setEnd(tempEnd);
                        tempIdRegionSet.add(tempIdRegion);
                    }
                    tempStart=tempId;
                    tempEnd=tempId;
                }
                else{
                    tempEnd=tempId;
                }
                //最后一个探针单独处理
                if((i==candidateIdSet.size()-1)&&((tempId-tempStart+1)> Parameters.minCNALength)){
                    IdRegion tempIdRegion=new IdRegion();
                    tempIdRegion.setStart(tempStart);
                    tempIdRegion.setEnd(tempId);
                    tempIdRegionSet.add(tempIdRegion);
                }
            }
        }

        //获得基于皮尔森相关系数分段后的CNA单元
        public void getIdRegion(ArrayList<IdRegion> tempIdRegionSet){
            int tempStart=0,tempEnd=0,tempId=0;
            double peaCorCoe=0.0;
            PearsonsCorrelation pearsonsCorrelation=new PearsonsCorrelation();
            for(IdRegion tempIdRegion:tempIdRegionSet){
                tempStart=tempIdRegion.getStart();
                tempEnd=tempIdRegion.getEnd();
                tempId=tempStart;
                while (tempId<tempEnd){
                    //peaCorCoe=getPeaCorCoe(tempId,tempId+1,candidateRawDatas);
                    peaCorCoe=pearsonsCorrelation.correlation(oneRawDataMatrix.getRow(tempId),
                            oneRawDataMatrix.getRow(tempId+1));
                    if(peaCorCoe<Parameters.pccThreshold){
                        if((tempId-tempStart+1)>=Parameters.minCNALength){
                            IdRegion validIdRegion=new IdRegion();
                            validIdRegion.setStart(tempStart);
                            validIdRegion.setEnd(tempId);
                            idRegionSet.add(validIdRegion);
                        }
                        tempStart=tempId+1;
                    }
                    tempId++;
                }
                //最后一个探针单独处理
                if ((tempId==tempEnd)&&((tempId-tempStart+1)>=Parameters.minCNALength)){
                    IdRegion validIdRegion=new IdRegion();
                    validIdRegion.setStart(tempStart);
                    validIdRegion.setEnd(tempId);
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
        public void regionMerge(){
            int cnaId=0;
            int cnaLength=0;

            //m_log.info("There are " + idRegionSet.size() + " CNA units:");
            if(idRegionSet.size()>0) {
               // m_log.info("cnaId" + "\t\t" + "start" + "\t\t" + "end" + "\t\t" + "length");

                for (IdRegion idRegion : idRegionSet) {
                    cnaLength = idRegion.getEnd() - idRegion.getStart() + 1;

                    CNARegion tempCNARegion = new CNARegion();
                    tempCNARegion.setCnaId(cnaId);
                    tempCNARegion.setIdRegion(idRegion);
                    tempCNARegion.setLength(cnaLength);

                    CNARegionSet.add(tempCNARegion);
                    lengthSet.add(cnaLength);
                    permuteProbeSize += cnaLength;
                   // m_log.info(cnaId + "\t\t" + idRegion.getStart() + "\t\t" + idRegion.getEnd() + "\t\t" + cnaLength);

                    cnaId++;
                }
            }
            lengthSet_copy=new ArrayList<>(lengthSet);
        }

        //显著性检测
        public void permuteDetection(){

            //获得CNA单元的唯一长度
            Set uniqueLength=new HashSet(lengthSet);
            uniqueLengthSet=new ArrayList(uniqueLength);
            Collections.sort(uniqueLengthSet);
            //ArrayList uniqueLengthSet_copy=new ArrayList(uniqueLengthSet);

            double[][] maxUScore=new double[Parameters.permuteNum][uniqueLengthSet.size()];

            calUScore(UScore);
            UScore_copy=new ArrayList<>(UScore);
            boolean loopFlag=true;
            ArrayList<CNARegion> permuteCNARegions=new ArrayList(CNARegionSet);
            while (loopFlag==true){
                permute(permuteCNARegions, maxUScore);
                loopFlag=SCAExclude(maxUScore,permuteCNARegions);
            }
            calPScore(PScore,maxUScore);
            getResultDatas();
        }

        //计算所有CNA单元的U值
        public void calUScore(List<Double> UScore){
            int start,end;
            int cnaLength;
            double tempUScore;

            for(IdRegion idRegion:idRegionSet){
                start=idRegion.getStart();
                end=idRegion.getEnd();
                cnaLength=end-start+1;
                tempUScore=0.0;
                for(int i=start;i<=end;i++){
//                    double[] temp=oneRawDataMatrix.getRow(i);
//                    double tempSum=getSum(temp);
                    tempUScore+=getSum(oneRawDataMatrix.getRow(i));
                }
                tempUScore=Math.abs(tempUScore/(cnaLength*(colNum-1)));
                UScore.add(tempUScore);
            }
        }

        public double getSum(double[] nums){
            double sum=0;
            for(int i=0;i<nums.length;i++){
                sum+=nums[i];
            }
            return sum;
        }

        //CNA单元交换，生成最大U值矩阵
        public void permute(List<CNARegion> permuteCNARegions,double[][] maxUScore){
            //int[][] idIndex=new int[permuteMatrix.length][sampleNum];
            for(int i=0;i<Parameters.permuteNum;i++){
                //一次交换开始
                RealMatrix randomPermuteMatrix=new BlockRealMatrix(permuteProbeSize,colNum-1);
                int maxIndex=permuteCNARegions.size();
                int index,len,cnaStart;

                //Generate sampleNum*permuteMatrix.length random numbers
                //生成 样本数*交换探针数 个随机数
                Random rn=new Random();
                for(int j=0;j<colNum-1;j++){
                    maxIndex=permuteCNARegions.size();
                    boolean[] bool=new boolean[maxIndex];
                    int randomInt=0;
                    int pStart=0;

                    //在每个样本之间进行交换
                    for(int k=0;k<permuteCNARegions.size();k++){

                        //生成不重复的随机数
                        do {
                            randomInt=rn.nextInt(maxIndex);
                        }while (bool[randomInt]);

                        bool[randomInt]=true;
                        index=randomInt;
                        len=permuteCNARegions.get(index).getLength();
                        cnaStart=permuteCNARegions.get(index).getIdRegion().getStart();

                        //根据序号和原有矩阵更新交换矩阵
//                        randomPermuteMatrix.setSubMatrix(
//                                oneRawDataMatrix.getSubMatrix(cnaStart-pStart,cnaStart+len-pStart-1,j,j).getData(),pStart,j);
                        randomPermuteMatrix.setSubMatrix(
                                  oneRawDataMatrix.getSubMatrix(cnaStart,cnaStart+len-1,j,j).getData(),pStart,j);
                        /*for(int m=pStart;m<(pStart+len);m++){
                            randomPermuteMatrix[m][j]=candidateRawDatas.get(cnaStart+m-pStart).getData()[j];
                        }*/
                        pStart+=len;
                    }
                }
                findMaxUScore(randomPermuteMatrix,maxUScore[i]);
            }
        }

        //找最大U值
        public void findMaxUScore(RealMatrix randomPermuteMatrix,double[] maxUScoreAti){
            double tempUScore;
            int probeNum=randomPermuteMatrix.getColumnDimension();
            double[] probeSum=new double[probeNum];
            int i,j,k;
            int len;
            int endPos;

            //计算每个探针的和
            for(i=0;i<probeNum;i++){
                probeSum[i]=0;
                probeSum[i]+=getSum(randomPermuteMatrix.getRow(i));
                /*for(j=0;j<colNum-1;j++){
                    probeSum[i]+=randomPermuteMatrix[i][j];
                }*/
                probeSum[i]=Math.abs(probeSum[i]);
            }

            //计算第一个唯一长度的U值
            len=uniqueLengthSet.get(0);
            tempUScore=0.0;
            maxUScoreAti[0]=0.0;
            endPos=probeNum-len;
            double[] regionSum=new double[endPos];
            for (j=0;j<endPos;j++){
                for (k=j;k<j+len;k++){
                    tempUScore+=probeSum[k];
                }
                regionSum[j]=Math.abs(tempUScore);
                tempUScore=Math.abs(tempUScore/((colNum-1)*len));
                if(tempUScore>maxUScoreAti[0]){
                    maxUScoreAti[0]=tempUScore;
                }
            }

            //计算余下唯一长度的U值
            for(i=1;i<uniqueLengthSet.size();i++){
                len=uniqueLengthSet.get(i);
                endPos=probeNum-len;
                maxUScoreAti[i]=0.0;
                for (j=0;j<endPos;j++){
                    tempUScore=regionSum[j];
                    for (k=j+uniqueLengthSet.get(i-1);k<j+len;k++){
                        tempUScore+=probeSum[k];
                    }
                    regionSum[j]=tempUScore;
                    tempUScore=tempUScore/((colNum-1)*len);
                    if(tempUScore>maxUScoreAti[i]){
                        maxUScoreAti[i]=tempUScore;
                    }
                }
            }
        }

        //排除SCAs并更新permuteCNARegions
        public boolean SCAExclude(double[][] maxUScore,ArrayList<CNARegion> permuteCNARegions){
            double sigValue;
            int index;
            boolean flag = false;
            int maxLen = Collections.max(uniqueLengthSet);//To make sure there are enough probes
            int permuteNum=Parameters.permuteNum;
            double sigValueThreshold=Parameters.sigValueThreshold;

            int permuteSize=permuteCNARegions.size();
            ArrayList<Integer> tempLengthSet=new ArrayList<>(lengthSet_copy);

            //根据长度搜索CNA单元的U值
            for (int i=0;i<permuteSize;i++){
                sigValue=0.0;
                index=uniqueLengthSet.indexOf(lengthSet_copy.get(i));

                //累积大于U值得最大U值得个数
                for(int j=0;j<permuteNum;j++){
                    if(maxUScore[j][index]>UScore.get(i)){
                        sigValue++;
                    }
                }
                sigValue/=(permuteNum+1);

                if(sigValue<sigValueThreshold){
                    flag=true;
                    if(permuteSize-tempLengthSet.get(i)>maxLen)
                        permuteProbeSize-=tempLengthSet.get(i);
                    permuteCNARegions.remove(i);
                    UScore_copy.remove(i);
                    lengthSet_copy.remove(i);
                    permuteSize=permuteCNARegions.size();
                }
            }
            return flag;
        }

        //再次计算所有CNA单元的P值
        public void calPScore(List<Double> PScore,double[][] maxUScore) {
            double pValue;
            int index;
            for (int i = 0; i < lengthSet.size(); i++) {
                pValue = 0.0;
                index = uniqueLengthSet.indexOf(lengthSet.get(i));
                for (int j = 0; j < Parameters.permuteNum; j++) {
                    if (maxUScore[j][index] > UScore.get(i)) {
                        pValue++;
                    }
                }
                pValue /= (Parameters.permuteNum + 1);
                PScore.add(pValue);
            }
        }

        //get ResultDatas and update the UScore,PScore and SCATag of the CNARegionSet
        //获得结果并更新CNARegionSet的U值，P值和SCATag
        public void getResultDatas(){
            double pScoreThreshold=Parameters.sigValueThreshold;
            Set<Region> tempResultRegions=new HashSet<>();
            //m_log.info("There are "+CNARegionSet.size()+" results:");
            //m_log.info("cnaId"+"\t\t"+"start"+"\t\t"+"end"+"\t\t"+"length"+"\t\t"+"uValue"+"\t\t\t\t"+"pValue"+"\t\t\t\t"+"SCATag");
            //update the CNARegionSet
            for (int i=0;i<CNARegionSet.size();i++){
                CNARegion tempCNARegion=new CNARegion();
                tempCNARegion=CNARegionSet.get(i);
                tempCNARegion.setuValue(UScore.get(i));
                tempCNARegion.setpValue(PScore.get(i));
                tempCNARegion.setSCATag(0);
                if(PScore.get(i)<pScoreThreshold){
                    tempCNARegion.setSCATag(1);
                    Region tempRegion=new Region();
                    tempRegion.setStartId(CNARegionSet.get(i).getIdRegion().getStart());
                    tempRegion.setEndId(CNARegionSet.get(i).getIdRegion().getEnd());
                    tempRegion.setLength(CNARegionSet.get(i).getLength());
                    tempResultRegions.add(tempRegion);
                }
                CNARegionSet.set(i,tempCNARegion);
//                m_log.info(i+"\t\t"+tempCNARegion.getIdRegion().getStart()+"\t\t"+tempCNARegion.getIdRegion().getEnd()+"\t\t"
//                        +tempCNARegion.getLength()+"\t\t"+tempCNARegion.getuValue()+"\t\t"+tempCNARegion.getpValue()+"\t\t"
//                        +tempCNARegion.getSCATag());
            }
            oneResultData.setRegionSet(tempResultRegions);
        }
    }

    class IdRegion {
        private int start;
        private int end;

        public int getStart(){
            return start;
        }
        public int getEnd(){
            return end;
        }

        public void setStart(int start){
            this.start=start;
        }
        public void setEnd(int end){
            this.end=end;
        }
    }

    class CNARegion {
        private int cnaId;
        private IdRegion idRegion;
        private int length;
        private double uValue;
        private double pValue;
        private int SCATag;

        public int getCnaId(){
            return cnaId;
        }
        public IdRegion getIdRegion(){
            return idRegion;
        }
        public int getLength(){
            return length;
        }
        public double getuValue(){
            return uValue;
        }
        public double getpValue(){
            return pValue;
        }
        public int getSCATag(){
            return SCATag;
        }

        public void setCnaId(int cnaId){
            this.cnaId=cnaId;
        }
        public void setIdRegion(IdRegion idRegion){
            this.idRegion=idRegion;
        }
        public void setLength(int length){
            this.length=length;
        }
        public void setuValue(double uValue){
            this.uValue=uValue;
        }
        public void setpValue(double pValue){
            this.pValue=pValue;
        }
        public void setSCATag(int SCATag){
            this.SCATag=SCATag;
        }
    }

    static class Parameters {
        static double pccThreshold = 0.9;
        static int permuteNum = 1000;
        static double sigValueThreshold = 0.0476;
        static int minCNALength = 6;
    }
}

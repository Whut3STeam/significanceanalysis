package edu.whut.significance.methods;

import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.ResultData;
import edu.whut.significance.methods.AbstractSig;
import edu.whut.significance.methods.GlobalParameters;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Created by SunMing on 2017/5/24.
 */
public class SAIC extends AbstractSig {
    ResultData resultData;
    RealMatrix rawMatrix;

    public void preprocess(RawData rawData){
        rawMatrix=rawData.getDataMatrix().transpose();
    }

    public void process(ResultData resultData){
        RealMatrix ampRawMatrix=new BlockRealMatrix(rawMatrix.getRowDimension(),rawMatrix.getColumnDimension());
        RealMatrix delRawMatrix=new BlockRealMatrix(rawMatrix.getRowDimension(),rawMatrix.getColumnDimension());
        ResultData ampResultData=new ResultData();
        ResultData delResultData=new ResultData();
        classify(rawMatrix,ampRawMatrix,delRawMatrix);
        permute(ampRawMatrix,ampResultData);
        permute(delRawMatrix,delResultData);
    }

    public void classify(RealMatrix rawMatrix,RealMatrix ampRawMatrix,RealMatrix delRawMatrix){
        GlobalParameters globalParameters=new GlobalParameters();
        double ampThreshold=globalParameters.getAmpThreshold();
        double delThreshold=globalParameters.getDelThreshold();
        int rowNum=rawMatrix.getRowDimension();
        int colNum=rawMatrix.getRow(0).length;

        ampRawMatrix=rawMatrix.copy();
        for(int i=0;i<rowNum;i++){
            for(int j=1;j<colNum;j++){
                if(ampRawMatrix.getEntry(i,j)<=ampThreshold)
                    ampRawMatrix.setEntry(i,j,0);
            }
        }

        delRawMatrix=rawMatrix.copy();
        for(int i=0;i<rowNum;i++){
            for(int j=1;j<colNum;j++){
                if(delRawMatrix.getEntry(i,j)>=delThreshold)
                    delRawMatrix.setEntry(i,j,0);
            }
        }
    }

    public void permute(RealMatrix rawMatrix,ResultData resultData){

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

    class Parameters {
        private double pccThreshold = 0.9;
        private int permuteNum = 1000;
        private double sigValueThreshold = 0.0476;
        private int minCNALength = 6;

        public double getPccThreshold() {
            return pccThreshold;
        }
        public int getPermuteNum() {
            return permuteNum;
        }
        public double getSigValueThreshold() {
            return sigValueThreshold;
        }
        public int getMinCNALength() {
            return minCNALength;
        }

        public void setPccThreshold(double pccThreshold) {
            this.pccThreshold = pccThreshold;
        }
        public void setPermuteNum(int permuteNum) {
            this.permuteNum = permuteNum;
        }
        public void setSigValueThreshold(double sigValueThreshold) {
            this.sigValueThreshold = sigValueThreshold;
        }
        public void setMinCNALength(int minCNALength) {
            this.minCNALength = minCNALength;
        }
    }
}

package edu.whut.significance.dataset;

import java.util.List;

/**
 * Created by SunMing on 2017/5/23.
 */
public class ExampleJ {
    List<Sample> samples;

    public List<Sample> getSamples() {
        return samples;
    }
    public void setSamples(List<Sample> samples) {
        this.samples = samples;
    }

    public static class Sample {
        String id;
        int length;//样本长度
        int sCount;//样本个数
        double baselineValue;//基线值
        double sigma;//强度噪声
        int winCount;//窗口数量
        Window[] windows;
        double[][] data;
        int breakPoint;//基线分割点
        double baselineInrcement;//第二个基线相对于第一个基线的变化值
        double baselineFluctuation;//基点的波动百分比
        int passengerCount;//乘客基因数量
        double passengerVal;//乘客的取值


        public static class Window {
            int midPos;
            int width;
            double LRRVal;
            double alpha;
            double beta;

            public int getMidPos() {
                return midPos;
            }

            public void setMidPos(int midPos) {
                this.midPos = midPos;
            }

            public int getWidth() {
                return width;
            }

            public void setWidth(int width) {
                this.width = width;
            }

            public double getLRRVal() {
                return LRRVal;
            }

            public void setLRRVal(double LRRVal) {
                this.LRRVal = LRRVal;
            }

            public double getAlpha() {
                return alpha;
            }

            public void setAlpha(double alpha) {
                this.alpha = alpha;
            }

            public double getBeta() {
                return beta;
            }

            public void setBeta(double beta) {
                this.beta = beta;
            }
        }

        public String getId() {
            return id;
        }
        public void setId(String id) {
            this.id = id;
        }
        
        public int getLength() {
            return length;
        }
        public void setLength(int length) {
            this.length = length;
        }
        
        public int getsCount() {
            return sCount;
        }
        public void setsCount(int sCount) {
            this.sCount = sCount;
        }
        
        public double getBaselineValue() {
            return baselineValue;
        }
        public void setBaselineValue(double baselineValue) {
            this.baselineValue = baselineValue;
        }
        
        public double getSigma() {
            return sigma;
        }
        public void setSigma(double sigma) {
            this.sigma = sigma;
        }
        public int getWinCount() {
            return winCount;
        }
        public void setWinCount(int winCount) {
            this.winCount = winCount;
        }

        public Window[] getWindows() {
            return windows;
        }
        public void setWindows(Window[] windows) {
            this.windows = windows;
        }
        
        public double[][] getData() {
            return data;
        }
        public void setData(double[][] data) {
            this.data = data;
        }


        public int getBreakPoint() {
            return breakPoint;
        }
        public void setBreakPoint(int breakPoint) {
            this.breakPoint = breakPoint;
        }

        public double getBaselineInrcement() {
            return baselineInrcement;
        }
        public void setBaselineInrcement(double baselineInrcement) {
            this.baselineInrcement = baselineInrcement;
        }

        public double getBaselineFluctuation() {
            return baselineFluctuation;
        }
        public void setBaselineFluctuation(double baselineFluctuation) {
            this.baselineFluctuation = baselineFluctuation;
        }

        public int getPassengerCount() {
            return passengerCount;
        }
        public void setPassengerCount(int passengerCount) {
            this.passengerCount = passengerCount;
        }
    }
}

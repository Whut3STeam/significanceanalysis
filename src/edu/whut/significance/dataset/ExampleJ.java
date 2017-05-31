package edu.whut.significance.dataset;

import java.util.List;

/**
 * Created by SunMing on 2017/5/23.
 */
public class ExampleJ {
    List<Sample> samples;

    public static class Sample {
        String id;
        int length;
        double value;
        double sigma;
        Windows windows;
        int count;
        double[][] data;

        public static class Windows {
            int midPos;
            int width;
            double value;
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

            public double getValue() {
                return value;
            }

            public void setValue(double value) {
                this.value = value;
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

        public double getValue() {
            return value;
        }

        public void setValue(double value) {
            this.value = value;
        }

        public double getSigma() {
            return sigma;
        }

        public void setSigma(double sigma) {
            this.sigma = sigma;
        }

        public Windows getWindows() {
            return windows;
        }

        public void setWindows(Windows windows) {
            this.windows = windows;
        }

        public int getCount() {
            return count;
        }

        public void setCount(int count) {
            this.count = count;
        }

        public double[][] getData() {
            return data;
        }

        public void setData(double[][] data) {
            this.data = data;
        }
    }

    public List<Sample> getSamples() {
        return samples;
    }

    public void setSamples(List<Sample> samples) {
        this.samples = samples;
    }
}

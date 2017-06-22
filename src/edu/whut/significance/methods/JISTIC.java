package edu.whut.significance.methods;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import edu.whut.significance.dataset.RawData;
import edu.whut.significance.dataset.Region;
import edu.whut.significance.dataset.ResultData;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.*;

/**
 * Created by SunMing on 2017/5/24.
 */
public class JISTIC extends AbstractSig {
    ResultData resultData;
    RealMatrix rawMatrix;
    int rowNum;
    int colNum;

    @Override
    public void preprocess(RawData rawData) {
        rawMatrix = rawData.getDataMatrix().getSubMatrix(1,
                rawData.getDataMatrix().getRowDimension()-1, 0, rawData.getDataMatrix().getColumnDimension()-1).transpose();
        rowNum = rawMatrix.getRowDimension();
        colNum = rawMatrix.getColumnDimension();
    }

    @Override
    public void process(ResultData resultData) {
        Distribution.run(rawMatrix);
        List<SignificantRegionJ> ampRegions = Distribution.regions.get(0);
        List<SignificantRegionJ> delRegions = Distribution.regions.get(1);

        mergeResult(resultData,ampRegions,delRegions);
//        System.out.println("OK");
    }

    public void mergeResult(ResultData resultData, List<SignificantRegionJ> ampRegions, List<SignificantRegionJ> delRegions){
        RangeSet<Integer> rangeSet= TreeRangeSet.create();
        for(SignificantRegionJ region:ampRegions){
            rangeSet.add(Range.closed(region.start,region.end));
        }
        for(SignificantRegionJ region:delRegions){
            rangeSet.add(Range.closed(region.start,region.end));
        }

        int start,end;
        for(Range range:rangeSet.asRanges()){
            start=(int)range.lowerEndpoint();//[422..596]的形式
            end=(int)range.upperEndpoint();
            Region tempRegion=new Region(start,end);
            resultData.getRegionSet().add(tempRegion);
        }
    }
}

class Distribution {
    static class sparse_cell implements Comparable<sparse_cell> {
        int gscore;
        int count;//位点总数
        double logcount;//位点总数取对数

        private static sparse_cell singletonref;

        public static sparse_cell getSingletonObject(int gscore) {
            if (singletonref == null)
                singletonref = new sparse_cell(gscore);
            singletonref.gscore = gscore;
            singletonref.count = 1;
            singletonref.logcount = 0;
            return singletonref;
        }

        public static sparse_cell getSingletonObject(sparse_cell c1, sparse_cell c2) {
            if (singletonref == null)
                // it's ok, we can call this constructor
                singletonref = new sparse_cell(c1.gscore + c2.gscore);
            singletonref.gscore = c1.gscore + c2.gscore;
            singletonref.count = -1;
            singletonref.logcount = c1.logcount + c2.logcount;
            return singletonref;
        }

        sparse_cell(int i) {
            gscore = i;
            count = 1;
            logcount = Double.NaN;
        }

        sparse_cell(sparse_cell c1, sparse_cell c2) {
            gscore = c1.gscore + c2.gscore;
            count = -1;
            logcount = c1.logcount + c2.logcount;
        }

        sparse_cell(sparse_cell c1) {
            gscore = c1.gscore;
            count = c1.count;
            logcount = c1.logcount;
        }


        void add(sparse_cell o) {
            if (count > 0)
                count += o.count;
            else {
                logcount = addLogs(logcount, o.logcount);
            }
        }

        static double addLogs(double a, double b) {
            double larger = Math.max(a, b);
            double smaller = Math.min(a, b);
            return Math.log1p(Math.exp(smaller - larger)) + larger;
        }

        public boolean equals(Object o) {
            return (o instanceof sparse_cell) && ((sparse_cell) o).gscore == gscore;
        }

        public int compareTo(sparse_cell o) {
            return o.gscore - gscore;
        }

        public String toString() {
            return String.valueOf(gscore) + ':' + Math.round(logcount);
        }

        public void computeLog() {
            logcount = Math.log(count);
        }
    }

    List<TreeSet<sparse_cell>> sample_distributions;
    TreeSet<sparse_cell> gscore_distribution;
    TreeMap<Integer, Double> gscore_to_qvalue = null;
    SNPdata data;
    static List<List<SignificantRegionJ>> regions;

    public static Distribution[] distributions = null;
    double log_total_perms;
    int num_markers;
    Marker.aberration type;


    static TreeSet<sparse_cell> convolution(TreeSet<sparse_cell> d1, TreeSet<sparse_cell> d2) {
        TreeSet<sparse_cell> result = new TreeSet<sparse_cell>();
        for (Iterator<sparse_cell> iter1 = d1.iterator(); iter1.hasNext(); ) {
            sparse_cell c1 = iter1.next();
            for (Iterator<sparse_cell> iter2 = d2.iterator(); iter2.hasNext(); ) {
                // obtain cell with sum of gscores and product of counts
                sparse_cell new_cell = sparse_cell.getSingletonObject(c1, iter2.next());
                // if cell exists add counts, if not add new cell
                sparse_cell old_cell = result.floor(new_cell);
                if (new_cell.equals(old_cell))
                    old_cell.add(new_cell);
                else
                    result.add(new sparse_cell(new_cell));
            }
        }
        return result;
    }

    public TreeSet<sparse_cell> convolution(List<TreeSet<sparse_cell>> dataconv) {
        TreeSet<sparse_cell> result = dataconv.get(0);
        for (ListIterator<TreeSet<sparse_cell>> iter = dataconv.listIterator(1); iter.hasNext(); ) {
            TreeSet<sparse_cell> nextitem = iter.next();
            result = convolution(result, nextitem);
        }
        return result;
    }

    public Distribution(SNPdata data, Marker.aberration type) {
        this.type = type;
        this.data = data;
        sample_distributions = data.sample_distributions.get(type.ordinal());
        int num_samples = sample_distributions.size();
        num_markers = data.cp_num_markers;
        // compute log for each cell
        for (Iterator<TreeSet<sparse_cell>> distr_iter = sample_distributions.iterator(); distr_iter.hasNext(); )
            for (Iterator<sparse_cell> iter = distr_iter.next().iterator(); iter.hasNext(); )
                iter.next().computeLog();
        // convolution across samples样本卷积
        gscore_distribution = convolution(sample_distributions);

        Iterator<sparse_cell> cell_iter = gscore_distribution.iterator();
        sparse_cell previous_cell = cell_iter.next();
        while (cell_iter.hasNext()) {
            sparse_cell next_cell = cell_iter.next();
            next_cell.add(previous_cell);
            previous_cell = next_cell;
        }
        log_total_perms = Math.log(num_markers) * num_samples;
    }

    static void run(RealMatrix rawMatrix) {
        SNPdata data = new SNPdata(rawMatrix);
        regions = new ArrayList<>();

        //计算空分布
        Marker.aberration[] types = Marker.aberration.values();
        distributions = new Distribution[types.length];
        Arrays.fill(Distribution.distributions, null);

        for (Marker.aberration ab : types) {
            distributions[ab.ordinal()] = new Distribution(data, ab);
            distributions[ab.ordinal()].calculateSignificance(data);//完成gscore_to_qvalue（Map）的赋值，给Gsig赋值
            if (SignificantRegionJ.peeloff_gscore_thres[ab.ordinal()] == 0 &&
                    !Double.isNaN(SignificantRegionJ.peeloff_qval_thres))
                if (Math.abs(SignificantRegionJ.peeloff_qval_thres - Parameters.qThres) < .0001)
                    SignificantRegionJ.peeloff_gscore_thres[ab.ordinal()] = Parameters.gThres;
                else {
                    Map.Entry<Integer, Double> entry;
                    Iterator<Map.Entry<Integer, Double>> iter = distributions[ab.ordinal()].gscore_to_qvalue.entrySet().iterator();
                    do
                        SignificantRegionJ.peeloff_gscore_thres[ab.ordinal()] = (entry = iter.next()).getKey().intValue();
                    while (entry.getValue().doubleValue() > SignificantRegionJ.peeloff_qval_thres);
                }

            //找显著性区域
            List<SignificantRegionJ> region = new Vector<SignificantRegionJ>();
            region.addAll(SignificantRegionJ.findRegions(distributions[ab.ordinal()], data.markers));
            regions.add(region);
        }
    }

    void calculateSignificance(SNPdata data) {
        gscore_to_qvalue = new TreeMap<Integer, Double>();
        List<Integer> byScore = data.byScore.get(type.ordinal());//byScore = new Vector<List<Integer>>();
        ListIterator<Integer> iter = byScore.listIterator();//byScore的个数是位点的个数

        //计算q值
        int significant_markers = 0;
        while (significant_markers < byScore.size()) {
            sparse_cell gscore_cell = gscore_distribution.floor(sparse_cell.getSingletonObject(iter.next()));
            double pval = Math.exp(gscore_cell.logcount - log_total_perms);
            while (significant_markers < byScore.size() &&
                    iter.previous().equals(byScore.get(significant_markers))) {
                significant_markers++;
                iter.next();
            }
            if (significant_markers == byScore.size())
                iter.previous();
            double qval = (pval * num_markers) / significant_markers;
            Integer m = null;
            while (iter.nextIndex() < significant_markers)//保存的是相同G值得最后一个位置
                m = iter.next();
            gscore_to_qvalue.put(m.intValue(), qval);
        }

        double previous_qval = 1;
        Map.Entry<Integer, Double> entry;
        for (Iterator<Map.Entry<Integer, Double>> map_iter = gscore_to_qvalue.entrySet().iterator(); map_iter.hasNext(); )
            if ((entry = map_iter.next()).getValue().doubleValue() > previous_qval)
                entry.setValue(previous_qval);
            else if ((previous_qval = entry.getValue().doubleValue()) < Parameters.qThres && Parameters.gThres == 0)
                Parameters.gThres = entry.getKey().intValue();


        for (ListIterator<Marker> marker_iter = data.listIterator(); marker_iter.hasNext(); ) {
            Marker marker = marker_iter.next();
            if (marker.gScore[type.ordinal()] >= 0) {
                entry = Distribution.distributions[type.ordinal()].gscore_to_qvalue.floorEntry(marker.gScore[type.ordinal()]);
                marker.qScore[type.ordinal()] = entry == null ? 1. : entry.getValue().doubleValue();
            }
        }
    }
}

class SNPdata extends Vector<Marker> {
    List<List<Integer>> byScore;
    List<List<TreeSet<Distribution.sparse_cell>>> sample_distributions;
    public int cp_num_markers = 0;
    List<Marker> markers;
    GlobalParameters globalParameters = new GlobalParameters();

    SNPdata(RealMatrix rawMatrix) {
        cp_num_markers = rawMatrix.getRowDimension();
        markers = new ArrayList<>();
        byScore = new Vector<List<Integer>>();
        sample_distributions = new Vector<List<TreeSet<Distribution.sparse_cell>>>();
        for (int i = 0; i < Marker.aberrationtypes.length; i++) {
            Marker.aberration ab = Marker.aberrationtypes[i];
            byScore.add(new Vector<Integer>());
            sample_distributions
                    .add(new Vector<TreeSet<Distribution.sparse_cell>>());
        }

        Marker.columns = rawMatrix.getColumnDimension();
        for (int i = 0; i < rawMatrix.getColumnDimension(); i++) {
            sample_distributions.get(0).add(
                    new TreeSet<Distribution.sparse_cell>());
            sample_distributions.get(1).add(
                    new TreeSet<Distribution.sparse_cell>());
        }

        for (int m = 0; m < rawMatrix.getRowDimension(); m++) {
               /*double[] data=rawMatrix.getRow(m);
               Marker new_marker = new Marker(data,m);*/
            Marker new_marker = new Marker(rawMatrix.getRow(m), m);
            byScore.get(0).add(new_marker.gScore[0]);
            byScore.get(1).add(new_marker.gScore[1]);
            for (int i = 0; i < Marker.columns; i++)
                for (int j = 0; j < 2; j++) {
                    double val = new_marker.copy_number[i];
                    if (j == 0 && val <= globalParameters.getAmpThreshold()
                            || j == 1
                            && val > globalParameters.getDelThreshold())
                        val = 0;
                    Distribution.sparse_cell new_cell = Distribution.sparse_cell.getSingletonObject((int) Math.round(Math.abs(val) / Parameters.binSize));
                    Distribution.sparse_cell old_cell = sample_distributions
                            .get(j).get(i).floor(new_cell);
                    if (new_cell.equals(old_cell))
                        old_cell.add(new_cell);
                    else
                        sample_distributions.get(j).get(i).add(new Distribution.sparse_cell(new_cell.gscore));
                }
            markers.add(new_marker);
        }


        for (int i = 0; i < Marker.aberrationtypes.length; i++) {
            Collections.sort(byScore.get(i));
            Collections.reverse(byScore.get(i));
        }
        addAll(markers);
    }
}

class Marker extends GenomeLocation {
    int id;
    double[] copy_number;
    int[] gScore;//0存扩增，1存缺失；q值同
    double[] qScore;
    static int columns = 0;
    GlobalParameters globalParameters = new GlobalParameters();

    public static aberration[] aberrationtypes = Marker.aberration.values();

    static enum aberration {
        AMP, DEL
    }

    ;

    Marker(double[] data, int id) {
        super(id, id);
        copy_number = new double[data.length];
        gScore = new int[aberrationtypes.length];
        qScore = new double[aberrationtypes.length];

        for (int i = 1; i < data.length; i++) {
            copy_number[i - 1] = data[i];
        }

        //计算G值
        Arrays.fill(gScore, -1);
        Arrays.fill(qScore, Double.NaN);
        int amp_sum = 0, del_sum = 0;
        for (int i = 0; i < copy_number.length; i++) {
            if (copy_number[i] > globalParameters.getAmpThreshold()) {
                amp_sum += Math.round(copy_number[i] / Parameters.binSize);//离散化
            } else if (copy_number[i] < globalParameters.getDelThreshold())
                del_sum -= Math.round(copy_number[i] / Parameters.binSize);
        }
        gScore[aberration.AMP.ordinal()] = amp_sum;
        gScore[aberration.DEL.ordinal()] = del_sum;

        //计算q值
        for (int i = 0; i < Marker.aberrationtypes.length; i++) {
            aberration type = Marker.aberrationtypes[i];
            if (gScore[type.ordinal()] >= 0 && Distribution.distributions != null) {
                Map.Entry<Integer, Double> entry = Distribution.distributions[type.ordinal()] == null ? null
                        : Distribution.distributions[type.ordinal()].gscore_to_qvalue.floorEntry(gScore[type.ordinal()]);
                qScore[type.ordinal()] = entry == null ? 1. : entry.getValue().doubleValue();
            }
        }
    }

    int aberrationStrength(aberration type, int sample) {
        switch (type) {
            case AMP:
                if (copy_number != null)
                    if (copy_number[sample] > globalParameters.getAmpThreshold())
                        return 1;
                break;

            case DEL:
                if (copy_number != null)
                    if (copy_number[sample] < globalParameters.getDelThreshold())
                        return 1;
                break;
        }
        return 0;
    }

    int GscoreContribution(aberration type, int sample) {
        return aberrationStrength(type, sample) == 0 ? 0 : (int) Math.abs(Math.round(copy_number[sample] / Parameters.binSize));
    }
}

class SignificantRegionJ extends RegionJ {
    List<RegionJ> subregions;
    List<PeakRegionJ> peaks;

    static int sparam = 11;//s值
    static int[] peeloff_gscore_thres = {0, 0};
    static double peeloff_qval_thres = Float.NaN;

    SignificantRegionJ(Distribution distribution, List<Marker> markers) {
        super(distribution.type, markers);
        subregions = new Vector<RegionJ>();
        int region_start = 0;
        while (region_start < markers.size()) {
            int region_end = region_start;
            while (region_end < markers.size()
                    && markers.get(region_end).gScore[type.ordinal()] == markers
                    .get(region_start).gScore[type.ordinal()])
                region_end++;
            subregions.add(new RegionJ(type, markers.subList(region_start, region_end)));
            region_start = region_end;
        }
        peaks = new Vector<PeakRegionJ>();
    }

    static List<SignificantRegionJ> findRegions(Distribution distribution, List<Marker> markers) {
        int num_samples = markers.get(0).copy_number.length;
        Marker.aberration type = distribution.type;
        List<SignificantRegionJ> regions = new Vector<SignificantRegionJ>();
        List<Set<Integer>> peeled_off = new Vector<Set<Integer>>();
        for (int i = 0; i < markers.size(); i++)
            peeled_off.add(new HashSet<Integer>());
        int[] Gscores = new int[markers.size()];
        double[] qvalues = new double[markers.size()];
        int[] peel_off_Gscores = new int[markers.size()];
        double[] peel_off_qvalues = new double[markers.size()];
        double peak_qvalue = 1;
        //求最小q值peak_qvalue
        for (Iterator<Marker> iter = markers.iterator(); iter.hasNext(); )
            peak_qvalue = Math.min(peak_qvalue, iter.next().qScore[type.ordinal()]);

        //Gscores，peel_off_Gscores，qvalues，peel_off_qvalues赋初值
        for (ListIterator<Marker> iter = markers.listIterator(); iter.hasNext(); ) {
            Marker marker = iter.next();
            Gscores[iter.previousIndex()] = peel_off_Gscores[iter.previousIndex()] = marker.gScore[type.ordinal()];
            qvalues[iter.previousIndex()] = peel_off_qvalues[iter.previousIndex()] = marker.qScore[type.ordinal()];
        }

        //表示存在峰值
        peak_loop:
        while (peak_qvalue < Parameters.qThres) {
            int raw_peak_start = 0;
            while (peel_off_qvalues[raw_peak_start] > peak_qvalue)//指针raw_peak_start移到第一个q值最小处，即第一个峰值的开始处
                raw_peak_start++;
            int raw_peak_end = raw_peak_start;
            while (raw_peak_end < markers.size()
                    && peel_off_qvalues[raw_peak_end] <= peak_qvalue)//继续移指针raw_peak_end到最后一个q值最小处，即第一个峰值的结束处
                raw_peak_end++;
            int peak_start = raw_peak_start, peak_end = raw_peak_end;
            for (int sample = 0; sample < num_samples; sample++) {
                if (!peeled_off.get(peak_start).contains(sample)) {
                    int one_out_peak_gscore = 0;
                    for (int marker = peak_start; marker < peak_end; marker++)//求临界Gn值，保存在one_out_peak_gscore中
                        one_out_peak_gscore = Math.max(one_out_peak_gscore,
                                peel_off_Gscores[marker] - markers.get(marker).GscoreContribution(type, sample));

                    //G值小于等于0时，退出程序
                    if (one_out_peak_gscore <= 0) {
                        RegionJ peak = new RegionJ(type, markers.subList(peak_start, peak_end));
                        System.exit(1);
                    }

                    //左右扩展区域到次峰，为了什么？？？？
                    int a = peel_off_Gscores[peak_start - 1];
                    int b = markers.get(peak_start - 1).GscoreContribution(type, sample);
//                    System.out.println("ok");
                    while (peak_start > 0 && peel_off_Gscores[peak_start - 1]
                            - markers.get(peak_start - 1).GscoreContribution(type, sample) >= one_out_peak_gscore
                            && peel_off_qvalues[peak_start - 1] <= Parameters.qThres)//向左扩展区域
                        peak_start--;
                    //subList()不包含结束位置
                    while (peak_end < markers.size() && peel_off_Gscores[peak_end]
                            - markers.get(peak_end).GscoreContribution(type, sample) >= one_out_peak_gscore
                            && peel_off_qvalues[peak_end] <= Parameters.qThres)//向右扩展区域
                        peak_end++;
                }
            }
            PeakRegionJ peak = new PeakRegionJ(type, markers.subList(peak_start, peak_end), peak_qvalue);
            SignificantRegionJ region = null;
            //查询峰是否在区域集中
            int region_index = Collections.binarySearch(regions, peak);
            if (region_index < 0)
                region_index = -region_index - 2;
            if (region_index >= 0 && regions.get(region_index).overlaps(peak))
                region = regions.get(region_index);
                //峰不在区域中，则生成新区域，并加到区域集中
            else {
                //继续扩展区域到最大，包含峰值周围所有原始q值大于qThres的所有位点
                int start = peak_start, end = peak_end;
                while (start > 0 && !(qvalues[start - 1] >= Parameters.qThres))
                    start--;
                while (Double.isNaN(qvalues[start]))
                    start++;
                while (end < markers.size() && !(qvalues[end] >= Parameters.qThres))
                    end++;
                while (Double.isNaN(qvalues[end - 1]))
                    end--;
                regions.add(region_index + 1, region = new SignificantRegionJ(
                        distribution, markers.subList(start, end)));
            }

            region.peaks.add(-Collections.binarySearch(region.peaks, peak) - 1,
                    peak);

            // limits for peeloff生成剥离的限制条件，即Gn值
            int max_peeloff[] = new int[2];

            // aberration remaining from peak for each sample
            int[] remain = new int[num_samples];
            // aberration added on top of remaining for each sample
            int[] added = new int[num_samples];
            int totalremain = 0;
            int totaladded = 0;

            // obtaining forward limit
            ListIterator<Marker> iter = markers.listIterator(peak_end - 1);
            max_peeloff[1] = iter.nextIndex();

            // initialize with peak end marker初始化峰结束的位点
            Marker marker = iter.next();
            for (int sample = 0; sample < num_samples; sample++) {
                remain[sample] = 0;
                if (!peeled_off.get(max_peeloff[1]).contains(sample))
                    remain[sample] = marker.GscoreContribution(type, sample);
                totalremain += remain[sample];
            }

            // iterate forward while under threshold
            while (iter.hasNext()) {
                totalremain = 0;
                totaladded = 0;
                max_peeloff[1] = iter.nextIndex();
                marker = iter.next();
                for (int sample = 0; sample < num_samples; sample++) {
                    int markerCont = 0;
                    if (!peeled_off.get(max_peeloff[1]).contains(sample)) {
                        markerCont = marker.GscoreContribution(type, sample);
                        if (markerCont < remain[sample]) {
                            int end = Math.min(max_peeloff[1] + sparam, markers.size());//避免超出整个样本的右边界，sparam是s值
                            for (int lfmid = max_peeloff[1]; lfmid < end; lfmid++) {
                                //求Gr值，向右
                                markerCont = Math.min(Math.max(markerCont, markers.get(lfmid).GscoreContribution(type, sample)), remain[sample]);
                            }
                        }
                    }

                    added[sample] = Math.max(markerCont - remain[sample], 0);//Gn值
                    remain[sample] = markerCont - added[sample];
                    totalremain += remain[sample];
                    totaladded += Math.max(markerCont - remain[sample], 0);
                }

                // finished if cutoff or no peak remaining
                if (totalremain == 0 || totaladded >= peeloff_gscore_thres[type.ordinal()]) {
                    break;
                }

            }

            // obtaining backward limit
            iter = markers.listIterator(peak_start + 1);

            // initialize with peak start marker初始化峰开始的位点
            max_peeloff[0] = iter.previousIndex();
            marker = iter.previous();
            totalremain = 0;
            for (int sample = 0; sample < num_samples; sample++) {
                remain[sample] = 0;
                if (!peeled_off.get(max_peeloff[0]).contains(sample))
                    remain[sample] = marker.GscoreContribution(type, sample);
                totalremain += remain[sample];
            }

            // iterate backwards while under threshold
            while (iter.hasPrevious()) {
                totalremain = 0;
                totaladded = 0;
                max_peeloff[0] = iter.previousIndex();
                marker = iter.previous();
                for (int sample = 0; sample < num_samples; sample++) {

                    int markerCont = 0;
                    if (!peeled_off.get(max_peeloff[0]).contains(sample)) {
                        markerCont = marker.GscoreContribution(type, sample);
                        if (markerCont < remain[sample]) {
                            int end = Math.max(max_peeloff[0] - sparam, 0);//避免超出整个样本的左边界，sparam是s值
                            for (int lfmid = max_peeloff[0]; lfmid >= end; lfmid--) {
                                //求Gr值，向左
                                markerCont = Math.min(Math.max(markerCont, markers.get(lfmid).GscoreContribution(type, sample)), remain[sample]);
                            }
                        }
                    }
                    added[sample] = Math.max(markerCont - remain[sample], 0);
                    remain[sample] = markerCont - added[sample];
                    totalremain += remain[sample];
                    totaladded += Math.max(markerCont - remain[sample], 0);
                }

                // finished if cutoff or no peak remaining
                if (totalremain == 0 || totaladded >= peeloff_gscore_thres[type.ordinal()]) {
                    break;
                }

            }

            peel_off_loop:
            for (int sample = 0; sample < num_samples; sample++) {
                for (ListIterator<Marker> iter1 = markers
                        .listIterator(peak_start); iter1.nextIndex() < peak_end; ) {
                    //如果样本对位点G值有贡献并且该位点的剥离样本里不包含此样本
                    if (iter1.next().aberrationStrength(type, sample) > 0
                            && !peeled_off.get(iter1.previousIndex()).contains(sample)) {
                        int peeloff_limits[];

                        peeloff_limits = new int[2];
                        while (iter1.hasPrevious() && iter1.previous().aberrationStrength(type, sample) > 0)
                            ;
                        if (iter1.hasNext() && iter1.next().aberrationStrength(type, sample) > 0)
                            iter1.previous();
                        peeloff_limits[0] = iter1.nextIndex();

                        while (iter1.hasNext() && iter1.next().aberrationStrength(type, sample) > 0)
                            ;
                        if (iter1.hasPrevious() && iter1.previous().aberrationStrength(type, sample) > 0)
                            iter1.next();
                        peeloff_limits[1] = iter1.previousIndex();

//                        peeloff_limits[0] = Math.max(peeloff_limits[0], max_peeloff[0]);
//                        peeloff_limits[1] = Math.min(peeloff_limits[1], max_peeloff[1]);
                        peeloff_limits[0] = Math.max(peeloff_limits[0], max_peeloff[0]+1);
                        peeloff_limits[1] = Math.min(peeloff_limits[1], max_peeloff[1]-1);

                        //这一步是要将被剥离的位点的G值置0
                        for (int marker1 = peeloff_limits[0]; marker1 <= peeloff_limits[1]; marker1++)
                            if (peeled_off.get(marker1).add(sample))//去掉样本对位点的G值的贡献
                                peel_off_Gscores[marker1] -= markers.get(marker1)
                                        .GscoreContribution(type, sample);
                        continue peel_off_loop;
                    }
                }
            }

            //重新更新下一次迭代的初始条件：最小P值，要进行剥离的位点的P值
            peak_qvalue = 1;
            for (int i = 0; i < markers.size(); i++) {
                Map.Entry<Integer, Double> entry = distribution.gscore_to_qvalue
                        .floorEntry(peel_off_Gscores[i]);
                peel_off_qvalues[i] = entry == null ? 1. : entry.getValue()
                        .doubleValue();
                peak_qvalue = Math.min(peak_qvalue, peel_off_qvalues[i]);//下一次迭代的最小q值为余下位点的最小q值
            }
        }
        return regions;
    }
}

class RegionJ extends GenomeLocation {
    Marker.aberration type;
    List<Marker> markers;
    double max_qScore;
    double min_qScore;


    RegionJ(Marker.aberration type, List<Marker> markers) {
        super(markers.get(0).start, markers.get(markers.size() - 1).end);
        this.type = type;
        this.markers = markers;
        max_qScore = 0;
        min_qScore = 1;
        for (Iterator<Marker> iter = markers.iterator(); iter.hasNext(); ) {
            Marker marker = iter.next();
            max_qScore = Math.max(max_qScore, marker.qScore[type.ordinal()]);
            min_qScore = Math.min(min_qScore, marker.qScore[type.ordinal()]);
        }

    }
}

class PeakRegionJ extends RegionJ {

    double peak_qvalue;
    boolean has_focal = true;

    PeakRegionJ(Marker.aberration type, List<Marker> markers, double peak_qvalue) {
        super(type, markers);
        this.peak_qvalue = peak_qvalue;
    }
}

class GenomeLocation implements Comparable<GenomeLocation> {
    int start;
    int end;

    GenomeLocation(int pos1, int pos2) {
        start = Math.min(pos1, pos2);
        end = Math.max(pos1, pos2);
    }

    @Override
    public int compareTo(GenomeLocation o) {
        return start == o.start ? ((end - o.end) >= 0 ? 0 : -1) : start - o.start;
        //return start == o.start ? end - o.end: start - o.start;
    }

    //完全相同的一个段或其子类，则返回true，否则返回false
    public boolean equals(Object o) {
        return (o instanceof GenomeLocation) &&
                same_location((GenomeLocation) o);
    }

    //完全相同的一个段则返回true，否则返回false
    public boolean same_location(GenomeLocation o) {
        return o.start == start
                && o.end == end;
    }

    //是否重叠
    public boolean overlaps(GenomeLocation o) {
        return start <= o.end && o.start <= end;
    }
}

class SigRegion {
    int id;
    RangeSet<Integer> region;
    int length;
}

class Parameters {
    static int gThres = 0;
    static double qThres = 0.25;
    static double binSize = 0.001;
}


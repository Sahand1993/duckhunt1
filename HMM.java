import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.Math;
import java.rmi.ServerError;
import java.util.ArrayList;
import java.util.List;

public class HMM {

    public BufferedReader br;
    public double[][] A;
    public double[][] B;
    public double[][] pi;
    private List<int[]> O;
    public double[][] alpha;
    public double[][] beta; // i, t
    private double[][][] digamma; // i, j, t
    private double[][] gamma; // i, t
    private int maxIters = 1000;
    public double[] colSums;
    private double logProb;
    public double[][] delta;
    private int[][] deltaIndex;

    public HMM() {
        O = new ArrayList<>();
    }

    public static void normalize(double[][] v) {
        double sum = 0;
        for(int i = 0; i < v.length; i++){
            for(int j = 0; j < v[0].length; j++) {
                sum += v[i][j];
            }
        }
        for(int i = 0; i < v.length; i++){
            for(int j = 0; j < v[0].length; j++) {
                v[i][j] /= sum;
            }
        }
    }

    public double fit(int[] O){
        estimateParams(O);
        double newLogProb = computeLogP();
        int iters = 0;
        do {
            //System.err.printf("logProb: %f\n", newLogProb);
            logProb = newLogProb;
            estimateParams(O);
            newLogProb = computeLogP();
            iters++;
        } while(iters < maxIters && (logProb - newLogProb) < -0.000001);
        /*System.err.printf("Stopped after %d iterations. with pi = ", iters);
        for(int i = 0; i < pi[0].length; i++){
            System.err.printf("%f ", pi[0][i]);
        }
        System.err.println();
        */
        return newLogProb;
    }

    public double computeLogP() {
        double logProb = 0.0;
        for(int i = 0; i < alpha[0].length; i++) {
            logProb += Math.log(1 / colSums[i]); // Does it matter which log-base is used?
        }
        logProb = -logProb;
        return logProb;
    }

    /**
     *
     * @param O The index of the O sequence we want to use for training.
     */
    public void estimateParams(int[] O) {
        fillAlpha(O);
        fillBeta(O);
        fillGammas(O);

        for(int i = 0; i < A.length; i++) {
            pi[0][i] = gamma[i][0];
        }

        double numer;
        double denom;
        for(int i = 0; i < A.length; i++) {
            for(int j = 0; j < A.length; j++) {
                numer = 0.0;
                denom = 0.0;
                for(int t = 0; t < O.length - 1; t++) {
                    numer += digamma[i][j][t];
                    denom += gamma[i][t];
                }
                numer += 0.000001;
                A[i][j] = numer / (denom + 0.000001 * A.length);
            }
        }

        for(int i = 0; i < B.length; i++) {
            for(int j = 0; j < B[0].length; j++) {
                numer = 0.0;
                denom = 0.0;
                for(int t = 0; t < O.length; t++) {
                    if(O[t] == j) {
                        numer += gamma[i][t];
                    }
                    denom += gamma[i][t];
                }
                numer += 0.000001;
                B[i][j] = numer / (denom + 0.000001 * B[0].length);
            }
        }
    }

    public void fillAlpha(int[] O) {
        alpha = new double[A.length][O.length];
        colSums = new double[O.length];
        // First column of alpha
        double colSum = 0;
        for (int i = 0; i < alpha.length; i++) {
            alpha[i][0] = pi[0][i] * B[i][O[0]];
            colSum += alpha[i][0];
            if(Double.isNaN(alpha[i][0])){
                System.err.printf("Calculated NaN in first column, alpha[%d][%d]\n", i, 0);
                System.err.println("pi:");
                printMatrix(pi);
                System.exit(1);
            }
        }
        colSums[0] = colSum;
        if(Double.isNaN(colSum)){
            System.err.println("Encountered colSum which is NaN in first column");
            System.exit(1);
        }
        normalizeCol(alpha, 0, colSum);

        double[][] alphaOld;
        double[][] newColAlpha;
        for (int i = 1; i < O.length; i++) {
            alphaOld = extractColumn(alpha, i - 1);
            newColAlpha = matrixMul(transpose(alphaOld), A);
            newColAlpha = vectorMul(transpose(newColAlpha), extractColumn(B, O[i]));
            colSum = 0;
            for(int j = 0; j < alpha.length; j++) {
                alpha[j][i] = newColAlpha[j][0];
                if(Double.isNaN(alpha[j][i])){
                    System.err.printf("Calculated NaN in alpha[%d][%d]\n", j, i);
                    System.exit(1);
                }
                colSum += alpha[j][i];
            }
            colSums[i] = colSum;
            if(Double.isNaN(colSum)){
                System.err.printf("Encountered colSum which is NaN in column %d\n", i);
                System.err.print("newColAlpha: ");
                for(int j = 0; j < alpha.length; j++) {
                    System.err.printf("%f ", newColAlpha[j][0]);
                }
                System.exit(1);
            }
            normalizeCol(alpha, i, colSum);
        }
    }

    private void fillBeta(int[] O) {
        beta = new double[A.length][O.length];
        // Last col of beta
        for (int i = 0; i < beta.length; i++) {
            beta[i][beta[0].length - 1] = 1 / colSums[O.length - 1];
            if(Double.isNaN(beta[i][beta[0].length - 1])){
                System.err.println("beta had NaN at middle:");
                System.err.printf("colSums[O.length - 1]: %d", colSums[O.length - 1]);
                System.exit(1);
            }
            if(Double.isInfinite(beta[i][beta[0].length - 1])) {
                System.err.println("beta was Infinity at middle:");
                System.err.printf("colSums[O.length - 1]: %d", colSums[O.length - 1]);
                System.exit(1);
            }

        }
        for (int t = O.length - 2; t >= 0; t--) {
            for (int i = 0; i < A.length; i++) {
                double sum = 0;
                for (int j = 0; j < A.length; j++) {
                    sum += beta[j][t + 1] * B[j][O[t + 1]] * A[i][j];
                    if(Double.isNaN(sum)){
                        System.err.println("sum was NaN");
                        System.err.printf("beta[j][t + 1], B[j][O[t + 1]], A[i][j]: %f, %f, %f\n", beta[j][t + 1], B[j][O[t + 1]], A[i][j]);
                        System.exit(1);
                    }
                }
                beta[i][t] = sum / colSums[t];
                if(Double.isNaN(beta[i][t])){
                    System.err.println("beta had NaN at end:");
                    System.err.printf("sum, colSums[t]: %f, %f\n", sum, colSums[t]);
                    System.exit(1);
                }
                if(Double.isInfinite(beta[i][t])){
                    System.err.println("beta was infinity at end:");
                    System.err.printf("sum, colSums[%d]: %f, %f\n", t, sum, colSums[t]);
                    printMatrix(alpha);
                    printVector(colSums);
                    System.err.printf("T: %d\n", alpha[0].length);
                    System.exit(1);
                }
            }
        }
    }

    private void printVector(double[] v) {
        for(int i = 0; i < v.length; i++) {
            System.err.print(v[i] + " ");
        }
        System.err.println();
    }

    /**
     * Fills delta and deltaIndex matrices with probabilities.
     */
    public void fillDelta(int[] O){
        delta = new double[A.length][O.length];
        deltaIndex = new int[A.length][O.length];

        firstColDelta(O);

        double[] probs;
        // For each timestep t
        for(int t = 1; t < O.length; t++) {
            // For each possible state at t
            for(int i = 0; i < A.length; i++) {
                probs = new double[A.length];
                // For each possible state at t - 1
                for(int j = 0; j < A.length; j++) {
                    probs[j] = delta[j][t - 1] * A[j][i] * B[i][O[t]];
                }
                DoubleInt doubleInt = max(probs);
                delta[i][t] = doubleInt.getMax(); // TODO: add 0.000001 to avoid underflow?
                deltaIndex[i][t] = doubleInt.getArgmax();
            }
        }
    }


    /**
     *  Fills the first column of delta
     */
    public void firstColDelta(int[] O){
        for(int i = 0; i < A.length; i++) {
            delta[i][0] = pi[0][i] * B[i][O[0]];
        }
    }

    /**
     * Return max of array
     */
    public DoubleInt max(double[] arr){
        double max = arr[0];
        int argmax = 0;
        for(int i = 1; i < arr.length; i++) {
            if (arr[i] > max) {
                max = arr[i];
                argmax = i;
            }
        }
        return new DoubleInt(max, argmax);
    }

    /**
     * Finds the most likely path by backtracking in deltaIndex.
     */
    public String findSequence(){
        String seq = "";
        double maxProb = 0;
        int iMax = -1;
        for (int i = 0; i < A.length; i++) {
            double tempProb = delta[i][delta[0].length - 1];
            if (tempProb > maxProb) {
                maxProb = tempProb;
                iMax = i;
            }
        }
        if (maxProb > 0) {
            seq = backTrack(iMax, delta[0].length - 1);
        }
        return seq;
    }

    /**
     *  Returns the sequence of states from backtracking as a String.
     * @param i is the state we're backtracking from.
     * @param t is the timestep where we have state i.
     * @return the entire backtracking sequence
     */
    public String backTrack(int i, int t){
        if (t == 0) {
            return "" + i + " ";
        }
        return backTrack(deltaIndex[i][t], t - 1) + i + " ";
    }

    /**
     * Returns last most likely state from viterbi algorithm (delta)
     * @return
     */
    public int lastState() {
        int state = -1;
        double maxProb = 0;
        double newProb;
        for(int i = 0; i < delta.length; i++) {
            newProb = delta[i][delta[0].length - 1];
            if(newProb > maxProb) {
                maxProb = newProb;
                state = i;
            }
        }
        return state;
    }


    public class DoubleInt {
        private double max;
        private int argmax;
        public DoubleInt(double max, int argmax) {
            this.max = max;
            this.argmax = argmax;
        }

        public double getMax() {
            return max;
        }

        public void setMax(double max) {
            this.max = max;
        }

        public int getArgmax() {
            return argmax;
        }

        public void setArgmax(int argmax) {
            this.argmax = argmax;
        }
    }

    private void fillGammas(int[] O) {
        gamma = new double[A.length][O.length];
        digamma = new double[A.length][A.length][O.length - 1];
        for (int t = 0; t < O.length - 1; t++) {
            for (int i = 0; i < A.length; i++) {
                gamma[i][t] = 0;
                for (int j = 0; j < A.length; j++) {
                    digamma[i][j][t] = alpha[i][t] * A[i][j] * B[j][O[t + 1]] * beta[j][t + 1];
                    gamma[i][t] += digamma[i][j][t];
                    if(Double.isNaN(gamma[i][t])) {
                        System.err.printf("encountered NaN in gamma[%d][%d]\n", i, t);
                        System.err.printf("%f, %f, %f, %f", alpha[i][t], A[i][j], B[j][O[t + 1]], beta[j][t + 1]);
                        System.exit(1);
                    }
                }
            }
        }
        for (int j = 0; j < A.length; j++) {
            gamma[j][O.length - 1] = alpha[j][O.length - 1];
            if(Double.isNaN(gamma[j][O.length - 1])) {
                System.err.printf("encountered NaN in gamma[%d][%d] at end of fillGammas\n", j, O.length - 1);
            }
        }
    }

    public static void printMatrix(double[][][] m) {
        for(int t = 0; t < m[0][0].length; t++) {
            for(int i = 0; i < m.length; i++) {
                for(int j = 0; j < m[0].length; j++){
                    System.err.printf("%f ", m[i][j][t]);
                }
                System.err.println();
            }
            System.err.println();
        }
        System.err.println();
    }

    private void normalizeCol(double[][] m, int colNo, double colSum) {
        for (int i = 0; i < m.length; i++) {
            m[i][colNo] /= colSum;
        }
    }

    private void printCol(double[][] m, int colNo){
        for (int i = 0; i < m.length; i++) {
            System.err.printf("%f\n", m[i][colNo]);
        }
    }


    public double[][] vectorMul(double[][] a, double[][] b){
        double[][] res = new double[a.length][1];
        for(int i = 0; i < a.length; i++) {
            res[i][0] = a[i][0] * b[i][0];
        }
        return res;
    }

    public static double[][] extractColumn(double[][] m, int colNo){
        double[][] res = new double[m.length][1];
        for(int i = 0; i < m.length; i++) {
            res[i][0] = m[i][colNo];
        }
        return res;
    }

    public static double[][] transpose(double[][] trans) {
        double[][] res = new double[trans[0].length][trans.length];
        for (int i = 0; i < trans.length; i++) {
            for (int j = 0; j < trans[0].length; j++) {
                res[j][i] = trans[i][j];
            }
        }
        return res;
    }
    public double sumMatrix(double[][] matrix) {
        double sum = 0;
        for (int i = 0; i < matrix.length; i++) {
            sum += matrix[i][0];
        }
        return sum;
    }
    public static void printMatrix(double[][] m){
        for(int i = 0; i < m.length; i++) {
            for( int j = 0; j < m[0].length; j++) {
                System.err.printf("%2f ", m[i][j]);
            }
            System.err.println();
        }
    }

    public static void printMatrix(int[][] m) {
        for(int i = 0; i < m.length; i++) {
            for( int j = 0; j < m[0].length; j++) {
                System.err.print(m[i][j] + " ");
            }
            System.err.println();
        }
        System.err.println();
    }

    public static void printVector(int[] v) {
        for(int i = 0; i < v.length; i++) {
            System.err.print(v[i] + " ");
        }
        System.err.println();
    }

    public static double[][] matrixMul(double[][] a, double[][] b) {
        double[][] res = new double[a.length][b[0].length];
        for(int i = 0; i < a.length; i++) {

            for(int j = 0; j < b[0].length; j++ ) {
                for (int k = 0; k < a[0].length; k++) {
                    res[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return res;
    }

    /**
     * Generates a random row of length n centered around 0 with uniform deviation from 0 up to 0.05.
     * @param n length of returned vector
     * @return random vector
     */
    private double[] randomRow(int n, double dev) {
        double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            dev = (Math.random() - 0.5) * dev;
            res[i] = 1 / n + dev;
        }
        return res;
    }

    /**
     * Creates and returns row-stochastic matrix
     * @param rows: number of rows
     * @param cols: number of columns
     * @return row-stochastic matrix
     */
    private double[][] randomMatrix(int rows, int cols, int DEV) {
        double[][] res = new double[rows][cols];
        double[] randomRow;
        double sum;
        for (int i = 0; i < rows; i++) {
            randomRow = randomRow(cols, DEV);
            sum = 0;
            for(int j = 0; j < cols; j++) {
                res[i][j] = 1.0 / cols + randomRow[j];
                sum += res[i][j];
            }
            for(int j = 0; j < cols; j++) {
                res[i][j] /= sum;
            }
        }
        return res;
    }

    /**
     *
     * @param N: number of hidden states
     * @param M: number of observation symbols
     */
    public double trainHMM(int N, int M, int[] O) {


        return fit(O);
    }

    public void randomizeParams(int N, int M, int DEV) {
        // Create random pi
        pi = randomMatrix(1, N, DEV);

        // Create random A
        A = randomMatrix(N, N, DEV);

        // Create random B
        B = randomMatrix(N, M, DEV);
    }

    public void initializeParams(int N, int M, int DEV) {
        pi = randomMatrix(1, N, DEV);
        A = identityMatrix(N);
        B = randomMatrix(N, M, DEV);
    }

    public void initializeParamsGuess(int N, int M, int DEV) {
        pi = randomMatrix(1, N, DEV); //TODO: change to uniform
        pi = new double[][]{{0.2, 0.2, 0.2, 0.2, 0.2}};
        A = identityMatrix(N);
        B = randomMatrix(N, M, DEV);
    }

    private double[][] identityMatrix(int n) {
        double[][] res = new double[n][n];
        for(int i = 0; i < n; i++){
            res[i][i] = 1;
        }
        return res;
    }

    public void addObservations(int[] OIn) {
        O.add(OIn);
    }

    public int sizeO() {
        return O.size();
    }

}

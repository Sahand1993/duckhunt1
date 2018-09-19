// Run with  and java Main verbose > player2server < server2player
import java.rmi.ServerError;
import java.util.*;
import java.lang.Integer;

class Player {

    private static final int N = 5;
    private static final int M = 9;
    Map<Integer, HMM> models = new HashMap<>(); // Species number mapping to HMM model.

    public Player() {
    }

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to get the best action.
         * This skeleton never shoots.
         */
        if (pState.getRound() < 1) {
            return cDontShoot;
        }

        int[][] observations = getObservations(pState);
        int[] O;
        int[] newO;
        double[][] scores = new double[Constants.COUNT_MOVE][Constants.COUNT_SPECIES];
        MaxScoreObj[] maxScores = new MaxScoreObj[pState.getNumBirds()];
        for(int i = 0; i < observations.length; i++){
            if(pState.getBird(i).isDead()){
                maxScores[i] = new MaxScoreObj(Double.NEGATIVE_INFINITY, Double.NaN, -1, -1, -1);
                continue;
            }
            O = observations[i];
            for(int j = 0; j < Constants.COUNT_MOVE; j++) {
                newO = addToArray(O, j);
                for(Integer key : models.keySet()) {
                    HMM model = models.get(key);
                    model.fillAlpha(newO);
                    scores[j][key] = model.computeLogP();
                }
            }
            MaxScoreObj maxScoreObj = createMaxScoreObj(scores, i);
            maxScores[i] = maxScoreObj;
        }

        MaxScoreObj totalMax = max(maxScores);
        if(totalMax.getStorkScore() > 0){
           // System.err.printf("NOT shooting since stork score was: %f\n", totalMax.getStorkScore());
            return cDontShoot;
        }else if(totalMax.getSpecies() == Constants.SPECIES_BLACK_STORK){
            //System.err.println("NOT shooting since stork gave highest probability");
            return cDontShoot;
        }
        else{
            //System.err.printf("\nshooting since stork score was: %f with logP: %f\n", totalMax.getStorkScore(), totalMax.getMaxScore() / pState.getBird(0).getSeqLength());
            //System.err.printf("birdNo = %d and move = %d\n", totalMax.getBirdNo(), totalMax.getMove());
            return cDontShoot; //new Action(totalMax.getBirdNo(), totalMax.getMove());
        }


        // This line would predict that bird 0 will move right and shoot at it.
        // return Action(0, MOVE_RIGHT);
    }

    private MaxScoreObj max(MaxScoreObj[] maxScores) {
        MaxScoreObj maxScore = maxScores[0];
        MaxScoreObj newMaxScore;
        for(int i = 1; i < maxScores.length; i++) {
            newMaxScore = maxScores[i];
            if(newMaxScore.compareTo(maxScore) > 0) {
                maxScore = newMaxScore;
            }
        }
        return maxScore;
    }

    private int[][] getObservations(GameState pState) {
        int[][] res = new int[pState.getNumBirds()][pState.getRound() + 1];
        for(int i = 0; i < pState.getNumBirds(); i++) {
            Bird bird = pState.getBird(i);
            int[] O = new int[bird.getSeqLength()];
            for (int j = 0; j < bird.getSeqLength(); j++) {
                O[j] = bird.getObservation(j);
            }
            res[i] = O;
        }
        return res;
    }

    public MaxScoreObj createMaxScoreObj(double[][] scores, int birdNo) {
        HMM.DoubleInt res;
        int rightMove = -1;
        int species = -1;
        double maxScore = Double.NEGATIVE_INFINITY;
        double storkScore;
        for(int i = 0; i < scores.length; i++) {
            for(Integer j : models.keySet()) {
                if(scores[i][j] > maxScore) {
                    maxScore = scores[i][j];
                    rightMove = i;
                    species = j;
                }
            }
        }
        if(models.containsKey(Constants.SPECIES_BLACK_STORK)) {
            storkScore = scores[rightMove][Constants.SPECIES_BLACK_STORK];
        }else {
            storkScore = Double.NaN;
        }
        return new MaxScoreObj(maxScore, storkScore, rightMove, birdNo, species);
    }

    private class MaxScoreObj implements Comparable<MaxScoreObj> {
        private double maxScore;
        private double storkScore;
        private int move;
        private int birdNo;
        private int species;

        @Override
        public int compareTo(MaxScoreObj o) {
            if(Double.isNaN(this.maxScore))
                return -1;
            return Double.compare(this.maxScore, o.getMaxScore());
        }

        public MaxScoreObj(double d1, double d2, int moveNo, int birdNo, int species) {
            maxScore = d1;
            storkScore = d2;
            move = moveNo;
            this.birdNo = birdNo;
            this.species = species;
        }

        public int getBirdNo() {
            return birdNo;
        }

        public void setBirdNo(int birdNo) {
            this.birdNo = birdNo;
        }

        public int getMove() {
            return move;
        }

        public double getMaxScore() {
            return maxScore;
        }

        public double getStorkScore() {
            return storkScore;
        }

        public String toString() {
            return String.format("move %d has maxProb %f and storkProb %f", move, maxScore, storkScore);
        }

        public int getSpecies() {
            return species;
        }
    }

    private int[] addToArray(int[] a, int b) {
        int[] res = new int[a.length + 1];
        int i;
        for(i = 0; i < a.length; i++) {
            res[i] = a[i];
        }
        res[i] = b;
        return res;
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to guess the species of
         * each bird. This skeleton makes no guesses, better safe than sorry!
         */

        int[] lGuess = new int[pState.getNumBirds()];
        if(pState.getRound() == 0) {
            for(int i = 0; i < pState.getNumBirds(); i++) {
                lGuess[i] = Constants.SPECIES_RAVEN;
            }
            return lGuess;
        }

        int[][] observations = getObservations(pState);
        HMM model;
        int[] O;
        double newScore, maxScore;

        for(int i = 0; i < pState.getNumBirds(); i++) {
            maxScore = Double.NEGATIVE_INFINITY;
            for(Integer species : models.keySet()) {
                model = models.get(species);
                O = observations[i];
                O = cutIfDead(O);
                model.fillAlpha(O);
                newScore = model.computeLogP();
                if (newScore > maxScore) {
                    lGuess[i] = species;
                    maxScore = newScore;
                }
            }
        }
        for(int i = 0; i < lGuess.length; i++){
            System.err.println(lGuess[i]);
        }
/*
        for (int i = 0; i < pState.getNumBirds(); ++i) {
            lGuess[i] = Constants.SPECIES_UNKNOWN;
        }*/

        return lGuess;
    }

    private int[] cutIfDead(int[] O) {
        for(int j = 0; j < O.length; j++){
            if(O[j] == -1) {
                O = Arrays.copyOfRange(O, 0, j);
                break;
            }
        }
        return O;
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
        for(int i = 0; i < pSpecies.length; i++) {
            System.err.printf("correct answer for guess %d: %d\n", i, pSpecies[i]);
        }
        List<Integer> unique = uniqueSpecies(pSpecies);
        for(int i = 0; i < unique.size(); i++) {
            for(int j = 0; j < pState.getNumBirds(); j++){
                //models.get(unique.get(i)).printMatrix(models.get(i).pi);
                if(pSpecies[j] == -1){
                    continue;
                }
                if(pSpecies[j] == unique.get(i)) {
                    // This is if the bird in pState is the same species as unique.get(i)
                    Bird bird = pState.getBird(j);
                    int[] O = new int[bird.getSeqLength()];
                    for(int k = 0; k < bird.getSeqLength(); k++) {
                        O[k] = bird.getObservation(k);
                    }
                    O = cutIfDead(O);
                    // if model for species exists then fit to new observation sequence, if absent: randomize params and fit
                    HMM model;
                    if(models.get(unique.get(i)) == null) {
                        models.put(unique.get(i), new HMM());
                        model = models.get(unique.get(i));
                        model.randomizeParams(N, M);
                    } else {
                        model = models.get(unique.get(i));
                    }
                    model.fit(O);
                }
            }
        }
    }

    public List<Integer> uniqueSpecies(int[] pSpecies) {
        HashMap<Integer, Boolean> hm = new HashMap<>();
        for(int i = 0; i < pSpecies.length; i++) {
            hm.putIfAbsent(pSpecies[i], true);
        }
        List<Integer> unique = new ArrayList<Integer>();
        for(Integer key : hm.keySet()) {
            unique.add(key);
        }
        Collections.sort(unique);
        return unique;
    }

    public static final Action cDontShoot = new Action(-1, -1);
}

// Run with  and java Main verbose > player2server < server2player
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
        for(int i = 0; i < 1; i++) {
            Bird bird = pState.getBird(i);
            for(int j = 0; j < bird.getSeqLength(); j++) {
                System.err.printf("%d ", bird.getObservation(j));
            }
            System.err.println();
        }

        int[][] observations = getObservations(pState);
        int[] O;
        int[] newO;
        double[][] scores = new double[Constants.COUNT_MOVE][Constants.COUNT_SPECIES];
        for(int i = 0; i < observations.length; i++){
            O = observations[i];
            for(int j = 0; j < Constants.COUNT_MOVE; j++) {
                newO = addToArray(O, j);
                for(Integer key : models.keySet()) {
                    HMM model = models.get(key);
                    model.fillAlpha(newO);
                    scores[j][key] = model.computeLogP(newO);
                    for (int k = 0; k < model.colSums.length; k++) {
                        System.err.printf("%f ", model.colSums[k]);
                    }
                    System.err.println();
                    System.err.printf("\nLogP: %f\n", model.computeLogP(newO));
                }
            }
        }

        MaxScoreObj maxScoreObj = createMaxScoreObj(scores);
        System.err.println(maxScoreObj.toString());
        return cDontShoot;

        // This line would predict that bird 0 will move right and shoot at it.
        // return Action(0, MOVE_RIGHT);
    }

    public MaxScoreObj createMaxScoreObj(double[][] scores) {
        HMM.DoubleInt res;
        int rightMove = -1;
        double maxScore = Double.NEGATIVE_INFINITY;
        double storkScore;
        for(int i = 0; i < scores.length; i++) {
            for(int j = 0; j < scores[0].length; j++) {
                if(scores[i][j] > maxScore) {
                    maxScore = scores[i][j];
                    rightMove = i;
                }
            }
        }
        storkScore = scores[rightMove][Constants.SPECIES_BLACK_STORK];
        return new MaxScoreObj(maxScore, storkScore, rightMove);
    }

    private class MaxScoreObj {
        private double maxScore;
        private double storkScore;
        private int move;

        public MaxScoreObj(double d1, double d2, int i) {
            maxScore = d1;
            storkScore = d2;
            move = i;
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
                lGuess[i] = Constants.SPECIES_PIGEON;
            }
            return lGuess;
        }

        for (int i = 0; i < pState.getNumBirds(); ++i) {
            lGuess[i] = Constants.SPECIES_UNKNOWN;
        }
        return lGuess;
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
            System.err.printf("%d ", unique.get(i));
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
                    // if model for species exists then fit to new observation sequence, if absent: randomize params and fit
                    HMM model;
                    if(models.get(unique.get(i)) == null) {
                        models.put(unique.get(i), new HMM());
                        model = models.get(unique.get(i));
                        model.randomizeParams(N, M);
                    } else {
                        System.err.println("there was already a model with pi: ");
                        model = models.get(unique.get(i));
                        model.printMatrix(model.pi);
                    }
                    model.fit(O);
                    System.err.println("pi from outside:");
                    model.printMatrix(model.pi);
                }
            }
        }
        System.err.println();

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

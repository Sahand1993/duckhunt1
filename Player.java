// Run with  and java Main verbose > player2server < server2player
import java.util.*;
import java.lang.Integer;

class Player {

    private static final int N = 5;
    private static final int M = 9;
    private static final double SHOOT_THRESHOLD = 0.6;
    Map<Integer, HMM> models = new HashMap<>(); // Species number mapping to HMM model.
    int[][] observationsPerBird;
    HMM[] shootModels;
    int[] likeliestLastState;
    Map<Integer, MoveScore> bird2MoveScore;
    List<ArrayList<HMM>> speciesModels; // speciesModels[species][modelNumber]

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

        // Shoot after 40 timesteps. Train one hmm per round and bird.
        // Calculate the probability of being at each state x for each hmm at time t.
        // Then for t+1, get the probability of each observation by taking
        // the probability of each state times the probability of transitioning to each other
        // state, times the probability of emitting the observation.

        // if it is after 40 timesteps
        if(pState.getBird(0).getSeqLength() < 50)
            return cDontShoot;

        // TODO: duck aversion?

        // create 1 model per bird in game and train it on the current data
        // TODO: if not getting good results, see if we should save the models instead of creating new ones every timestep
        shootModels = new HMM[pState.getNumBirds()];
        bird2MoveScore = new HashMap<>();
        //likeliestLastState = new int[shootModels.length];
        int lastState;
        double[][] stateVector;
        double[][] emissionVector;
        MoveScore moveScore;

        outerloop:
        for(int i = 0; i < pState.getNumBirds(); i++) {

            shootModels[i] = new HMM();
            shootModels[i].randomizeParams(N, M); // TODO: check if this initialization is optimal for shooting. Uniform pi?
            shootModels[i].pi = new double[][]{{0.2, 0.2, 0.2, 0.2, 0.2}};
            int[] O = new int[pState.getBird(i).getSeqLength()];
            for(int j = 0; j < pState.getBird(i).getSeqLength(); j++) {
                int nextMove = pState.getBird(i).getObservation(j);
                if(nextMove == -1) break outerloop; // if the bird is dead, we don't care about any more of its observations
                O[j] = pState.getBird(i).getObservation(j);
            }
            shootModels[i].fit(O);
            // for each model, calculate most likely state for bird with viterbi
            shootModels[i].fillDelta(O);
            // Get likeliest last state
            lastState = shootModels[i].lastState();
            // create a vector with 1 on that state and 0 everywhere else
            stateVector = new double[1][N];
            stateVector[0][lastState] = 1.0;
            // multiply that vector with A to get state probabilities at next timestep.
            stateVector = HMM.matrixMul(stateVector, shootModels[i].A);
            // multiply the state probabilities at next timestep with B to get observation probabilities
            emissionVector = HMM.matrixMul(stateVector, shootModels[i].B);
            // Shoot on the first emission probability that is over 0.75.
            double emissionProb;
            double maxProb = 0;
            int move = -1;
            for(int j = 0; j < emissionVector[0].length; j++) {
                emissionProb = emissionVector[0][j];
                if(emissionProb > SHOOT_THRESHOLD){
                    System.err.println("shot");
                    return new Action(i, j);
                }
            }
            // put the move and the probability in a move per bird vector and shoot on the highest score after this loop

        }



        // create vector with 1 on that state and 0 everywhere else
        // multiply that vector with A to get state probabilities at next timestep.
        // multiply the state probabilities at next timestep with B to get observation probabilities

        // take max out of observation probabilities vector. If it's over 0.75, shoot.

        return cDontShoot;


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
        System.err.printf("round %d in guess\n", pState.getRound());
        if(pState.getRound() == 0){
            for(int i = 0; i < lGuess.length; i++){
                lGuess[i] = Constants.SPECIES_PIGEON;
            }
            return lGuess;
        }

        double logP;
        int speciesGuess;
        // TODO:
        // for each bird in pState,
        for(int i = 0; i < pState.getNumBirds(); i++){
            // TODO: get the observation sequence of that bird
            int[] O = new int[pState.getBird(i).getSeqLength()];
            for(int j = 0; j < pState.getBird(i).getSeqLength(); j++){
                int newObs = pState.getBird(i).getObservation(j);
                if(newObs == -1)
                    break;
                O[j] = newObs;
            }
            // for each HMM in speciesModels[][], do a logP score.
            speciesGuess = -1;
            double maxLogP = Double.NEGATIVE_INFINITY;
            for(int j = 0; j < speciesModels.size(); j++){ // j is the species
                for(int k = 0; k < speciesModels.get(j).size(); k++) { // k is not interesting for the guessing
                    speciesModels.get(j).get(k).fillAlpha(O);
                    logP = speciesModels.get(j).get(k).computeLogP();
                    // is this log p greater than the greatest logP?
                    if(logP > maxLogP) {
                        maxLogP = logP;
                        speciesGuess = j;
                    }
                }
            }
            // The highest score will get its species guessed
            lGuess[i] = speciesGuess;
        }
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
        System.err.printf("round %d in reveal\n", pState.getRound());
        if(pState.getRound() == 0){
            speciesModels = new ArrayList<>();
            // add an empty list for each species
            for(int i = 0; i < pSpecies.length; i++){
                speciesModels.add(new ArrayList<>());
            }
        }

        for(int i = 0; i < pSpecies.length; i++) {
            System.err.printf("correct answer for guess %d: %d\n", i, pSpecies[i]);
        }

        // For each bird in pSpecies
        for(int i = 0; i < pSpecies.length; i++) {
            // generate a new model
            HMM model = new HMM();
            model.initializeParamsGuess(N, M);
            // train the model on the observation data from this round
            if(pState.getBird(i).isAlive()) {
                int[] O = new int[pState.getBird(i).getSeqLength()];
                for (int j = 0; j < pState.getBird(i).getSeqLength(); j++) {
                    O[j] = pState.getBird(i).getObservation(j);
                }
                model.fit(O);
            }// TODO: should we have an else for broken observation sequences due to dead bird?
            // add the model to speciesModels.get(species)
            System.err.printf("pSpecies[i]: %d\n", pSpecies[i]);
            speciesModels.get(pSpecies[i]).add(model);
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

    private class MoveScore {
        int move;
        double score;
        public MoveScore(int move, double score) {
            this.move = move;
            this.score = score;
        }

        public double getScore() {
            return score;
        }

        public int getMove() {
            return move;
        }
    }
}

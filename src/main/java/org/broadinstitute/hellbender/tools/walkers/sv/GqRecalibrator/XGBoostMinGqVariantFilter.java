package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import htsjdk.variant.vcf.VCFConstants;
import ml.dmlc.xgboost4j.java.*;
import ml.dmlc.xgboost4j.java.util.BigDenseMatrix;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

@CommandLineProgramProperties(
        summary = "Extract matrix of properties for each variant. Also extract, num_variants x num_trios x 3 tensors of" +
                "allele count and genotype quality. These data will be used to train a variant filter based on min GQ" +
                "(and stratified by other variant properties) that maximizes the admission of variants with Mendelian" +
                "inheritance pattern while omitting non-Mendelian variants." +
                "Derived class must implement abstract method train_filter()",
        oneLineSummary = "Extract data for training min GQ variant filter from .vcf and .ped file with trios.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class XGBoostMinGqVariantFilter extends MinGqVariantFilterBase {
    @Argument(fullName="max-training-rounds", shortName="tr", doc="Maximum number of rounds of training", optional=true, minValue=1)
    public int maxTrainingRounds = 100;
    @Argument(fullName="early-stopping-rounds", shortName="e", doc="Stop training if no improvement is made in validation set for this many rounds. Set <= 0 to disable.", optional=true)
    public int earlyStoppingRounds = 10;
    @Argument(fullName="learning-rate", shortName="lr", doc="Learning rate for xgboost", optional=true)
    public double eta = 0.1;
    @Argument(fullName="max-depth", shortName="d", doc="Max depth of boosted decision tree", optional=true, minValue=1)
    public int maxDepth = 10;
    @Argument(fullName="gamma", doc="Regularization factor for xgboost", optional=true)
    public double gamma = 1.0e-3;
    @Argument(fullName="subsample", doc="Proportion of data selected for each tree", optional=true)
    public double subsample = 0.9;
    @Argument(fullName="colsample-by-tree", doc="Proportion of columns selected for each tree", optional=true)
    public double colsampleByTree = 0.9;
    @Argument(fullName="colsample-by-level", doc="Proportion of columns selected for each level of each tree", optional=true)
    public double colsampleByLevel = 1.0;
    @Argument(fullName="min-child-weight", doc="Minimum sum of hessian weight for each node of tree", optional=true, minValue=0.0)
    public double minChildWeight = 10000;

    private Booster booster = null;
    private int gqPropertyIndex = -1;

    private static final String TRAIN_MAT_KEY = "train";
    private static final String VALIDATION_MAT_KEY = "validation";


    protected float[] getRowMajorSampleProperties(int[] variantIndices) {
        // Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
        final int numRows = getNumTrainableSampleVariants(variantIndices);

        final float[] rowMajorVariantProperties = new float[numRows * getNumProperties()];
        // Loop over variants and filterable samples. Store properties for each sample in a single flat array
        final int numSamples = getNumSamples();
        int flatIndex = 0;
        for(final int variantIndex : variantIndices) {
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                flatIndex = propertiesTable.copyPropertiesRow(
                    rowMajorVariantProperties, flatIndex, variantIndex, sampleIndex, needsNormalizedProperties()
                );
            }
        }

        return rowMajorVariantProperties;
    }

    protected BigDenseMatrix getVariantSamplePropertiesMatrix(int[] variantIndices) {
        // Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
        final int numRows = getNumTrainableSampleVariants(variantIndices);

        final int numProperties = getNumProperties();
        final BigDenseMatrix variantSampleProperties = new BigDenseMatrix(numRows, numProperties);
        // Loop over variants and filterable samples. Store properties for each sample in a single flat array
        final int numSamples = getNumSamples();
        int row = 0;
        for(final int variantIndex : variantIndices) {
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                final float[] sampleProperties = propertiesTable.getPropertiesRow(variantIndex, sampleIndex,
                                                                                  needsNormalizedProperties());
                for(int column = 0; column < sampleProperties.length; ++column) {
                    variantSampleProperties.set(row, column, sampleProperties[column]);
                }
                ++row;
            }
        }

        return variantSampleProperties;
    }

    private DMatrix getDMatrix(final int[] variantIndices) {
        if(progressVerbosity > 1) {
            System.out.println("\tgetDMatrix");
        }
//        final float[] propertiesValues = getRowMajorSampleProperties(variantIndices);
//        if(progressVerbosity > 1) {
//            System.out.format("\t\tgot %d property entries\n", propertiesValues.length);
//        }
        final BigDenseMatrix propertiesValues = getVariantSamplePropertiesMatrix(variantIndices);
        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d property entries\n",
                    propertiesValues.nrow * (long)propertiesValues.ncol);
        }
        final boolean[] sampleVariantTruth = getSampleVariantTruth(variantIndices);
        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d sampleVariantTruth status\n", sampleVariantTruth.length);
        }
        // balance true / false categories
        final int numRows = sampleVariantTruth.length;
        final int numTrue = (int)IntStream.range(0, numRows).filter(i -> sampleVariantTruth[i]).count();
        final int numFalse = numRows - numTrue;
        final double trueWeight = (1 + numFalse) / (double)(numRows + 2);
        final double falseWeight = (1 + numTrue) / (double)(numRows + 2);

        final float[] sampleVariantWeights = new float[numRows];
        final int numSamples = getNumSamples();        int row = 0;
        for(int variantIndex : variantIndices) {
            final int propertyBin = propertyBins[variantIndex];
            final Set<Integer> inheritanceTrainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
            final Set<Integer> truthTrainableSampleIndices = new HashSet<>(
                goodSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet())
            );
            truthTrainableSampleIndices.addAll(
                badSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet())
            );

            final float minGqWeight = (float)minBinWeight; // everyone gets some weight for being partially based on minGq
            final float truthWeight = minGqWeight + (float)propertyBinTruthWeights[propertyBin];
            final float inheritanceWeight = minGqWeight + (float)propertyBinInheritanceWeights[propertyBin];
            final float bothWeights = minGqWeight + (float)(propertyBinTruthWeights[propertyBin] +
                                                            propertyBinInheritanceWeights[propertyBin]);

            for(int sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
                if(getSampleVariantIsTrainable(variantIndex, sampleNum)) {
                    if(inheritanceTrainableSampleIndices.contains(sampleNum)) {
                        if(truthTrainableSampleIndices.contains(sampleNum)) {
                            sampleVariantWeights[row] = bothWeights;
                        } else {
                            sampleVariantWeights[row] = inheritanceWeight;
                        }
                    } else if(truthTrainableSampleIndices.contains(sampleNum)) {
                        sampleVariantWeights[row] = truthWeight;
                    } else {  // this sample is only trainable because we've set minGQ. Take the mean of the two weights
                        sampleVariantWeights[row] = minGqWeight;
                    }
                    sampleVariantWeights[row] *= sampleVariantTruth[row] ? trueWeight : falseWeight;
                    ++row;
                }
            }
        }
        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d row weights\n", row);
        }

        return getDMatrix(propertiesValues, sampleVariantTruth, sampleVariantWeights, numRows, getNumProperties());
    }

    private DMatrix getDMatrixForPrediction(final float[] sampleVariantProperties, final int numRows, final int numColumns) {
        final boolean[] sampleIsGood = new boolean[numRows];
        final float[] variantWeight = new float[numRows];
        return getDMatrix(sampleVariantProperties, sampleIsGood, variantWeight, numRows, numColumns);
    }

    private DMatrix getDMatrix(final float[] propertiesArr, final boolean[] sampleVariantTruth,
                               final float[] variantWeights, final int numRows, final int numColumns) {
        try {
            for(final float val : propertiesArr) {
                if(!Float.isFinite(val)) {
                    throw new GATKException("rowMajorVariantProperties contains a non-finite value (" + val + ")");
                }
            }
            final DMatrix dMatrix = new DMatrix(
                    propertiesArr, numRows, numColumns, Float.NaN
            );
            if(progressVerbosity > 2) {
                System.out.format("\t\tcreated DMatrix object\n");
            }
            // Set baseline (initial prediction for min GQ)
            if(gqPropertyIndex < 0) {
                gqPropertyIndex = propertiesTable.getPropertyNames().indexOf(VCFConstants.GENOTYPE_QUALITY_KEY);
            }
            int baselineIndex = gqPropertyIndex;
            final float[] baseline = new float[numRows];
            final float[] labels = new float[numRows];
            for(int idx = 0; idx < numRows; ++idx) {
//                final float gq = propertiesArr[baselineIndex];
//                baseline[idx] = probToLogits(
//                    FastMath.max(0.1F, FastMath.min((float)phredToProb(gq), 0.9F))
//                );
                baseline[idx] = 0F;
                baselineIndex += numColumns;
                labels[idx] = sampleVariantTruth[idx] ? 1F : -1F;
            }
            dMatrix.setLabel(labels);
            dMatrix.setBaseMargin(baseline);
            dMatrix.setWeight(variantWeights);
            if(progressVerbosity > 2) {
                System.out.println("\tgetDMatrix complete");
            }
            return dMatrix;
        }
        catch(XGBoostError xgBoostError) {
            throw new GATKException("Error forming DMatrix", xgBoostError);
        }
    }

    private DMatrix getDMatrix(final BigDenseMatrix propertiesMatrix, final boolean[] sampleVariantTruth,
                               final float[] variantWeights, final int numRows, final int numColumns) {
        try {
            final DMatrix dMatrix = new DMatrix(propertiesMatrix, Float.NaN);
            if(progressVerbosity > 2) {
                System.out.format("\t\tcreated DMatrix object\n");
            }
            // Set baseline (initial prediction for min GQ)
            if(gqPropertyIndex < 0) {
                gqPropertyIndex = propertiesTable.getPropertyNames().indexOf(VCFConstants.GENOTYPE_QUALITY_KEY);
            }
            int baselineIndex = gqPropertyIndex;
            final float[] baseline = new float[numRows];
            final float[] labels = new float[numRows];
            for(int idx = 0; idx < numRows; ++idx) {
//                final float gq = propertiesArr[baselineIndex];
//                baseline[idx] = probToLogits(
//                    FastMath.max(0.1F, FastMath.min((float)phredToProb(gq), 0.9F))
//                );
                baseline[idx] = 0F;
                baselineIndex += numColumns;
                labels[idx] = sampleVariantTruth[idx] ? 1F : -1F;
            }
            dMatrix.setLabel(labels);
            dMatrix.setBaseMargin(baseline);
            dMatrix.setWeight(variantWeights);
            if(progressVerbosity > 2) {
                System.out.println("\tgetDMatrix complete");
            }
            return dMatrix;
        }
        catch(XGBoostError xgBoostError) {
            throw new GATKException("Error forming DMatrix", xgBoostError);
        }
    }

    private Map<String, Object> getXgboostParams() {
        return new HashMap<String, Object>() {
            private static final long serialVersionUID = 0L;
            {
                put("eta", eta);
                put("max_depth", maxDepth);
                put("gamma", gamma);
                put("subsample", subsample);
                put("colsample_bytree", colsampleByTree);
                put("colsample_bylevel", colsampleByLevel);
                put("min_child_weight", minChildWeight);
                put("validate_parameters", true);
                put("objective", "binary:logistic");
                put("eval_metric", "logloss");
            }
        };
    }

    private class DataSubset {
        final String name;
        final int[] variantIndices;
        final boolean isTrainingSet;
        final int numTrainableSampleVariants;

        final List<Float> scores;
        private int bestScoreInd;
        private float bestScore;

        final DMatrix dMatrix;
        final float[] pSampleVariantGood;
        final float[] d1Loss;
        final float[] d2Loss;

        DataSubset(final String name, final int[] variantIndices, final boolean isTrainingSet) {
            System.out.println("Building DataSubset " + name);
            this.name = name;
            this.variantIndices = variantIndices;
            this.isTrainingSet = isTrainingSet;
            numTrainableSampleVariants = getNumTrainableSampleVariants(variantIndices);
            if(progressVerbosity > 1) {
                System.out.format("\t%d variants, %d sample x variants, %d properties, %d entries\n",
                                   variantIndices.length, numTrainableSampleVariants, getNumProperties(),
                                   numTrainableSampleVariants * (long)getNumProperties());
            }

            scores = new ArrayList<>();
            bestScore = Float.POSITIVE_INFINITY;
            bestScoreInd = 0;

            this.dMatrix = getDMatrix(variantIndices);
            pSampleVariantGood = new float[numTrainableSampleVariants];

            // pre-allocate derivative arrays if they are needed (i.e. this is the training set)
            d1Loss = isTrainingSet ? new float[numTrainableSampleVariants] : null;
            d2Loss = isTrainingSet ? new float[numTrainableSampleVariants] : null;

            final boolean[] sampleVariantTruth = getSampleVariantTruth(variantIndices);
            final float[] tempTruthProbs = new float[numTrainableSampleVariants];
            for(int idx = 0; idx < numTrainableSampleVariants; ++idx) {
                tempTruthProbs[idx] = sampleVariantTruth[idx] ? 1F : 0F;
            }
            System.out.println("Best possible " + name + " loss:\n\t" +
                               getLoss(tempTruthProbs, variantIndices).toString().replaceAll("\n", "\n\t"));
        }

        public int size() { return numTrainableSampleVariants; }

        private float[][] getRawPredictions(final Booster booster, final DMatrix dMatrix) {
            try {
                return booster.predict(dMatrix, true, 0);
            } catch(XGBoostError xgBoostError) {
                throw new GATKException("In " + name + " DataSubset: Predict error", xgBoostError);
            }
        }

        private void displayQuantilesFloat(final float[] arr, final String description) {
            DoubleStream.Builder builder = DoubleStream.builder();
            for(int idx = 0; idx < arr.length; ++idx) {
                builder.add(arr[idx]);
            }
            displayPercentiles(description, builder.build());
        }

        @SuppressWarnings("SameParameterValue")
        private void displayFloatHistogram(final float[] arr, final String description, final int numBins) {
            double minValue = arr[0];
            double maxValue = arr[0];
            DoubleStream.Builder builder = DoubleStream.builder();
            final Set<Double> nonFinite = new HashSet<>();
            for(final float val : arr) {
                if(val < minValue) {
                    minValue = val;
                } else if(val > maxValue) {
                    maxValue = val;
                } else if(!Double.isFinite(val)){
                    nonFinite.add((double)val);
                }
                builder.add(val);
            }
            if(maxValue - minValue < FastMath.min(1e-6, 1e-3 * FastMath.abs(maxValue + minValue) / 2)) {
                System.out.println(description);
                System.out.format("\t%f: 100%%\n", minValue);
            } else {
                displayHistogram(description, builder.build(), numBins, minValue, maxValue);
            }
            if(!nonFinite.isEmpty()) {
                System.out.format("%s has non-finite values:", description);
                for(final Double value : nonFinite) {
                    System.out.format(" %f", value);
                }
                System.out.println();
            }
        }

        private FilterLoss getRoundLoss(final Booster booster, final boolean training) {
            final float[][] rawPredictions = getRawPredictions(booster, dMatrix);
            for(int idx = 0; idx < rawPredictions.length; ++idx) {
                pSampleVariantGood[idx] = scaled_logits_to_p(rawPredictions[idx][0]);
            }

            if(isTrainingSet && training) {
                final FilterLoss lossForTraining = getTrainingLoss(pSampleVariantGood, d1Loss, d2Loss, variantIndices);
                System.out.format("\t%s\n\t\t%s\n", "optimizer objective", lossForTraining.toString().replaceAll("\n", "\n\t\t"));
                if(progressVerbosity > 0) {
                    System.out.format("\tBoosting round %d on %s\n", 1 + getRound(), name);
                }
                if(progressVerbosity > 1) {
                    displayQuantilesFloat(pSampleVariantGood, "pSampleVariantIsGood");
                    displayQuantilesFloat(d1Loss, "d1Loss");
                    displayQuantilesFloat(d2Loss, "d2Loss");
                }
                try {
                    booster.boost(dMatrix, d1Loss, d2Loss);
                } catch(XGBoostError xgBoostError) {
                    throw new GATKException("In " + name + " DataSubset: Boost error", xgBoostError);
                }
            }

            final FilterLoss lossForDiagnostics = getLoss(pSampleVariantGood, variantIndices);


            if(progressVerbosity > 0) {
                System.out.format("\t%s\n\t\t%s\n", name, lossForDiagnostics.toString().replaceAll("\n", "\n\t\t"));
            }
            return lossForDiagnostics;
        }

        public DoubleStream streamPredictions(final Booster booster) {
            return Arrays.stream(getRawPredictions(booster, dMatrix))
                .mapToDouble(rawPrediction -> scaled_logits_to_p(rawPrediction[0]));
        }


        public float getLastScore() { return scores.get(scores.size() - 1); }

        public void appendScore(final float score) {
            scores.add(score);
            if(score < bestScore) {
                bestScore = score;
                bestScoreInd = getRound();
            }
        }

        public int getRound() {
            return scores.size() - 1;
        }

        public boolean isBestScore() {
            return bestScoreInd == getRound();
        }

        public boolean stop() {
            return getRound() == maxTrainingRounds || stopEarly();
        }

        public boolean stopEarly() {
            return earlyStoppingRounds > 0 && getRound() - bestScoreInd > earlyStoppingRounds;
        }
    }

    @Override
    protected boolean needsNormalizedProperties() { return false; }

    @Override
    protected void predictBatch(final float[] variantProperties, final int numRows, final int numColumns,
                                final float[] outputProbabilities) {
        try {
            final float[][] rawPredictions = booster.predict(
                getDMatrixForPrediction(variantProperties, numRows, numColumns), true, 0
            );
            for(int row = 0; row < numRows; ++ row) {
                outputProbabilities[row] = scaled_logits_to_p(rawPredictions[row][0]);
            }
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error predicting", xgBoostError);
        }
    }

    @Override
    protected void loadModel(final InputStream inputStream) {
        try {
            booster = XGBoost.loadModel(inputStream);
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error loading XGBoost model", xgBoostError);
        }
    }

    @Override
    protected void saveModel(final OutputStream outputStream) {
        try {
            booster.saveModel(outputStream);
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error saving XGBoost model", xgBoostError);
        }
    }

    private Booster initializeBooster(final List<DataSubset> dataSubsets) {
        final DataSubset trainingSubset = dataSubsets.stream()
            .filter(dataSubset -> dataSubset.isTrainingSet)
            .findFirst()
            .orElseThrow(() -> new GATKException("dataSubsets does not contain a training set"));
        final Map<String, DMatrix> watches = new HashMap<>();
        // add watches for any DataSubset that a) is not the main training set, and
        //                                     b) does not use mini batches (it can comfortably hang out in memory)
        for(final DataSubset dataSubset : dataSubsets) {
            if(!dataSubset.equals(trainingSubset)) {
                watches.put(dataSubset.name, dataSubset.dMatrix);
            }
        }

        try {
            // Do 0 rounds of training on the first mini-batch (if using, otherwise whole dMatrix) of training DataSubset
            return XGBoost.train(trainingSubset.dMatrix, getXgboostParams(), 0, watches, null, null);
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error creating Booster", xgBoostError);
        }
    }

    private boolean trainOneRound(final Booster booster, final List<DataSubset> dataSubsets) {
        // evaluate booster on all data sets, calculate derivatives on training set
        if(progressVerbosity > 0) {
            System.out.println("Training round " + (dataSubsets.get(0).getRound() + 1));
        }

        boolean needSaveCheckpoint = !dataSubsets.isEmpty();
        boolean needStop = !dataSubsets.isEmpty();
        for(final DataSubset dataSubset : dataSubsets) {
            final FilterLoss roundLoss = dataSubset.getRoundLoss(booster, true);
            dataSubset.appendScore(roundLoss.toFloat());
            if(!dataSubset.isTrainingSet) {
                // check if need to stop iterating or save model checkpoint
                needSaveCheckpoint = needSaveCheckpoint && dataSubset.isBestScore();
                needStop = needStop && dataSubset.stop();
            }
        }

        // check if booster needs to be saved, or if early stopping is necessary
        if(needSaveCheckpoint) {
            saveModelCheckpoint();
        }
        return !needStop;
    }

    void displayFeatureImportance(final Booster booster) {
        final List<String> propertyNames = propertiesTable.getPropertyNames();
        final Map<String, Integer> featureScore;
        try {
            featureScore = booster.getFeatureScore((String) null).entrySet().stream()
                .collect(
                    Collectors.toMap(
                        entry -> propertyNames.get(Integer.parseInt(entry.getKey().substring(1))),
                        Map.Entry::getValue
                    )
                );
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error getting feature score", xgBoostError);
        }

        System.out.println("Feature importance:");
        featureScore.entrySet().stream()
            .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
            .forEachOrdered(entry -> System.out.format("\t%s: %d\n", entry.getKey(), entry.getValue()));
    }

    @Override
    protected void trainFilter() {
        final List<DataSubset> dataSubsets = new ArrayList<DataSubset>() {
            private static final long serialVersionUID = 0L;
            {
                add(new DataSubset(TRAIN_MAT_KEY, getTrainingIndices(), true));
                add(new DataSubset(VALIDATION_MAT_KEY, getValidationIndices(), false));
            }
        };

        clearUnneededProperties();

        booster = initializeBooster(dataSubsets);
        //noinspection StatementWithEmptyBody
        while (trainOneRound(booster, dataSubsets))
            ;

        loadModelCheckpoint();
        if(progressVerbosity > 0) {
            System.out.println("Losses of final trained classifier:");
            for (final DataSubset dataSubset : dataSubsets) {
                dataSubset.getRoundLoss(booster, false);
            }
        }
        displayHistogram("Final training adjusted GQ histogram",
                          dataSubsets.get(0).streamPredictions(booster).mapToInt(this::p_to_scaled_logits),true);
        displayHistogram("Final training probability histogram",
                          dataSubsets.get(0).streamPredictions(booster),20, 0.0, 1.0);
        displayFeatureImportance(booster);
    }
}

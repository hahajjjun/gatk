package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Retrieves PE/SR evidence and performs breakpoint refinement
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 *     <li>
 *         Mean depth table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk AggregatePairedEndAndSplitReadEvidence
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Adds PE/SR evidence to structural variant records",
        oneLineSummary = "Adds PE/SR evidence to structural variant records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class AggregatePairedEndAndSplitReadEvidence extends TwoPassVariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String PE_INNER_WINDOW_LONG_NAME = "pe-inner-window";
    public static final String PE_OUTER_WINDOW_LONG_NAME = "pe-outer-window";
    public static final String SR_WINDOW_LONG_NAME = "sr-window";
    public static final String X_CHROMOSOME_LONG_NAME = "x-chromosome-name";
    public static final String Y_CHROMOSOME_LONG_NAME = "y-chromosome-name";

    @Argument(
            doc = "Split reads evidence file",
            fullName = SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file",
            fullName = DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first column and expected per-base coverage per sample " +
                    "in the second column. Required if a split read evidence file is provided.",
            fullName = SAMPLE_COVERAGE_LONG_NAME,
            optional = true
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = PE_INNER_WINDOW_LONG_NAME,
            optional = true
    )
    private int innerWindow = 50;

    @Argument(
            doc = "Outer discordant pair window size (bp)",
            fullName = PE_OUTER_WINDOW_LONG_NAME,
            optional = true
    )
    private int outerWindow = 500;

    @Argument(
            doc = "Split read window size (bp)",
            fullName = SR_WINDOW_LONG_NAME,
            optional = true
    )
    private int splitReadWindow = 200;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLES and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv). Required only if the input VCF contains allosomal records.",
            fullName = SVCluster.PLOIDY_TABLE_LONG_NAME,
            optional = true
    )
    private GATKPath ploidyTablePath;

    @Argument(
            doc = "X chromosome name",
            fullName = X_CHROMOSOME_LONG_NAME,
            optional = true
    )
    private String xChromosomeName = "chrX";

    @Argument(
            doc = "Y chromosome name",
            fullName = Y_CHROMOSOME_LONG_NAME,
            optional = true
    )
    private String yChromosomeName = "chrY";

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private BreakpointRefiner breakpointRefiner;
    private DiscordantPairEvidenceTester discordantPairEvidenceTester;
    private DiscordantPairEvidenceAggregator discordantPairCollector;
    private SplitReadEvidenceAggregator startSplitCollector;
    private SplitReadEvidenceAggregator endSplitCollector;
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private VCFHeader header;
    private Set<String> maleSamples;
    private Set<String> femaleSamples;

    private Collection<SimpleInterval> discordantPairIntervals;
    private Collection<SimpleInterval> splitReadIntervals;

    private Collection<VariantContext> outputBuffer;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        if (!splitReadCollectionEnabled() && !discordantPairCollectionEnabled()) {
            throw new UserException.BadInput("At least one evidence file must be provided");
        }
        loadSampleCoverage();
        if (discordantPairCollectionEnabled()) {
            initializeDiscordantPairCollection();
        }
        if (splitReadCollectionEnabled()) {
            initializeSplitReadCollection();
        }
        discordantPairIntervals = new ArrayList<>();
        splitReadIntervals = new ArrayList<>();
        outputBuffer = new ArrayList<>();
        writer = createVCFWriter(Paths.get(outputFile));
        header = getVCFHeader();
        writer.writeHeader(header);
        if (ploidyTablePath != null) {
            initializeSampleSexSets();
        }
    }

    private void initializeSampleSexSets() {
        final PloidyTable table = new PloidyTable(ploidyTablePath.toPath());
        maleSamples = header.getGenotypeSamples().stream()
                .filter(s -> table.get(s, xChromosomeName) == 1 && table.get(s, yChromosomeName) == 1)
                .collect(Collectors.toSet());
        femaleSamples = header.getGenotypeSamples().stream()
                .filter(s -> table.get(s, xChromosomeName) == 2 && table.get(s, yChromosomeName) == 0)
                .collect(Collectors.toSet());
    }

    private void initializeDiscordantPairCollection() {
        initializeDiscordantPairDataSource();
        discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, innerWindow, outerWindow);
        discordantPairEvidenceTester = new DiscordantPairEvidenceTester(sampleCoverageMap, dictionary);
    }

    private void initializeSplitReadCollection() {
        if (sampleCoverageFile == null) {
            throw new UserException.BadInput("Sample coverage table is required when a split read evidence file is present");
        }
        initializeSplitReadEvidenceDataSource();
        startSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, true);
        endSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, splitReadWindow, false);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap, dictionary);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile.toString(),
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile1",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        final String fileString = sampleCoverageFile.toString();
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(fileString, e);
        }
    }

    @Override
    public void closeTool() {
        if (discordantPairSource != null) {
            discordantPairSource.close();
        }
        if (splitReadSource != null) {
            splitReadSource.close();
        }
        if (writer != null) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecordUtils.create(variant);
        discordantPairIntervals.add(discordantPairCollector.getEvidenceQueryInterval(call));
        splitReadIntervals.add(startSplitCollector.getEvidenceQueryInterval(call));
    }

    @Override
    public void afterFirstPass() {
        discordantPairCollector.setCacheIntervals(discordantPairIntervals);
        startSplitCollector.setCacheIntervals(splitReadIntervals);
    }

    /**
     * Returns set of samples to exclude for evidence stat calculations by sex.
     * @param record
     * @return set of sample ids or null
     */
    public Set<String> getSamplesToExcludeForStatsBySex(final SVCallRecord record) {
        // TODO paired BND records may cause problems here
        if (record.getContigA().equals(xChromosomeName)) {
            Utils.validate(femaleSamples != null, "Ploidy table must be provided to call on X chromosome");
            final Set<String> carrierSamples = record.getCarrierSampleSet();
            final Collection<String> femaleCarriers = Sets.intersection(carrierSamples, femaleSamples);
            if (femaleCarriers.isEmpty()) {
                // Exclude no samples for X chromosome when there are no female carriers
                return Collections.emptySet();
            } else {
                // If there are female carriers on X, exclude males
                return maleSamples;
            }
        } else if (record.getContigA().equals(yChromosomeName)) {
            Utils.validate(maleSamples != null, "Ploidy table must be provided to call on Y chromosome");
            // Always exclude females for Y chromosome
            return femaleSamples;
        } else {
            // Default autosome
            return Collections.emptySet();
        }
    }

    @Override
    public void secondPassApply(final VariantContext variant, final ReadsContext readsContext,
                                final ReferenceContext referenceContext, final FeatureContext featureContext) {
        SVCallRecord record = SVCallRecordUtils.create(variant);
        final Set<String> excludedSamples = getSamplesToExcludeForStatsBySex(record);
        flushOutputBuffer(record.getPositionAInterval());
        if (discordantPairCollectionEnabled()) {
            final List<DiscordantPairEvidence> discordantPairEvidence = discordantPairCollector.collectEvidence(record);
            final EvidenceStatUtils.PoissonTestResult result = discordantPairEvidenceTester.poissonTestRecord(record, discordantPairEvidence, excludedSamples);
            final double p = result == null ? 1. : result.getP();
            final Integer q = EvidenceStatUtils.probToQual(p, (byte) 99);
            final Integer carrierSignal = result == null ? null :
                    EvidenceStatUtils.carrierSignalFraction(result.getCarrierSignal(), result.getBackgroundSignal());
            final Map<String, Object> attributes = new HashMap<>();
            attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, q);
            attributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, carrierSignal);
            record = SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
            record = SVCallRecordUtils.assignDiscordantPairCountsToGenotypes(record, discordantPairEvidence);
        }
        if (splitReadCollectionEnabled()) {
            final List<SplitReadEvidence> startSplitReadEvidence = startSplitCollector.collectEvidence(record);
            final List<SplitReadEvidence> endSplitReadEvidence = endSplitCollector.collectEvidence(record);
            record = breakpointRefiner.refineCall(record, startSplitReadEvidence, endSplitReadEvidence, excludedSamples);
        }
        outputBuffer.add(SVCallRecordUtils.getVariantBuilder(record).make());
    }

    @Override
    public Object onTraversalSuccess() {
        outputBuffer.stream().sorted(IntervalUtils.getDictionaryOrderComparator(dictionary)).forEach(writer::add);
        outputBuffer.clear();
        return super.onTraversalSuccess();
    }

    private void flushOutputBuffer(final SimpleInterval currentLocus) {
        outputBuffer.stream()
                .filter(v -> !variantIsActive(v, currentLocus))
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .forEach(writer::add);
        outputBuffer = outputBuffer.stream()
                .filter(v -> variantIsActive(v, currentLocus))
                .collect(Collectors.toList());
    }

    private boolean variantIsActive(final VariantContext variant, final SimpleInterval currentLocus) {
        return variant.getContig().equals(currentLocus.getContig()) && variant.getStart() >= currentLocus.getStart() - splitReadWindow;
    }

    private VCFHeader getVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        if (splitReadCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.START_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read quality at start of variant"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END_SPLIT_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read quality at end of variant"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.START_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample split read signal at start of variant"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample split read signal at end of variant"));
        }
        if (discordantPairCollectionEnabled()) {
            header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair quality"));
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_CARRIER_SIGNAL_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Carrier sample discordant pair signal"));
        }
        return header;
    }

    private boolean splitReadCollectionEnabled() {
        return splitReadsFile != null;
    }

    private boolean discordantPairCollectionEnabled() {
        return discordantPairsFile != null;
    }
}

package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class DiscordantPairEvidenceAggregatorTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/printevidence/test_hg38.pe.txt.gz";

    @DataProvider(name = "getDiscordantPairIntervalTestData")
    public Object[][] getDiscordantPairIntervalTestData() {
        return new Object[][]{
                {1100, true, 200, 500, 600, 1300},   // positive strand
                {1100, false, 200, 500, 900, 1600},  // negative strand
        };
    }

    @Test(dataProvider= "getDiscordantPairIntervalTestData")
    public void getDiscordantPairIntervalTest(final int pos, final boolean strand,
                                              final int innerWindow, final int outerWindow,
                                              final int expectedStart, final int expectedEnd) {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY,
                innerWindow, outerWindow);
        final String contig = "chr1";

        final SVCallRecord startRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos, strand, contig, pos + 100, true);
        final SimpleInterval startWindow = aggregator.getDiscordantPairStartInterval(startRecord);
        Assert.assertEquals(startWindow.getContig(), contig);
        Assert.assertEquals(startWindow.getStart(), expectedStart);
        Assert.assertEquals(startWindow.getEnd(), expectedEnd);

        final SVCallRecord endRecord = SVTestUtils.newBndCallRecordWithPositionAndStrands(contig, pos - 100, true, contig, pos, strand);
        final SimpleInterval endWindow = aggregator.getDiscordantPairEndInterval(endRecord);
        Assert.assertEquals(endWindow.getContig(), contig);
        Assert.assertEquals(endWindow.getStart(), expectedStart);
        Assert.assertEquals(endWindow.getEnd(), expectedEnd);
    }

    @Test
    public void testGetters() {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);
        final int innerWindow = 100;
        final int outerWindow = 200;
        final DiscordantPairEvidenceAggregator aggregator = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow);
        Assert.assertEquals(aggregator.getInnerWindow(), innerWindow);
        Assert.assertEquals(aggregator.getOuterWindow(), outerWindow);
    }

    @DataProvider(name = "testCollectEvidenceData")
    public Object[][] testCollectEvidenceData() {
        return new Object[][] {
                // 1-based coordinates

                // Single pair, +- strands
                {25047743, true, 25048449, false, 0, 0, 1},
                {25047743, true, 25048449, true, 0, 0, 0},  // wrong strands
                {25047743, false, 25048449, true, 0, 0, 0},  // wrong strands
                {25047743, false, 25048449, false, 0, 0, 0},  // wrong strands
                {25047744, true, 25048448, false, 0, 1, 1},  // outer window +1
                {25047742, true, 25048450, false, 0, 1, 0},
                {25047742, true, 25048450, false, 1, 0, 1},  // inner window +1
                {25047744, true, 25048448, false, 1, 0, 0},

                // Single pair, -- strands
                {25052409, false, 25052644, false, 0, 0, 1},
                {25052409, true, 25052644, false, 0, 0, 0},  // wrong strands
        };
    }

    @Test(dataProvider= "testCollectEvidenceData")
    public void testCollectEvidence(final int posA, final boolean strandA, final int posB, final boolean strandB,
                                    final int innerWindow, final int outerWindow, final int expectedCount) {
        final FeatureDataSource<DiscordantPairEvidence> source = new FeatureDataSource<>(TEST_EVIDENCE);

        final SVCallRecord record = SVTestUtils.newBndCallRecordWithPositionAndStrands("chr21", posA, strandA,
                "chr21", posB, strandB);

        final Collection<SimpleInterval> cacheIntervals = new ArrayList<>();
        cacheIntervals.add(new SimpleInterval("chr21", 25004918, 25005918));
        if (strandA) {
            cacheIntervals.add(new SimpleInterval("chr21", posA - outerWindow, posA + innerWindow));
        } else {
            cacheIntervals.add(new SimpleInterval("chr21", posA - innerWindow, posA + outerWindow));
        }

        // No caching
        final List<DiscordantPairEvidence> testStart = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow)
                .collectEvidence(record);
        Assert.assertEquals(testStart.size(), expectedCount);

        // With caching
        final DiscordantPairEvidenceAggregator cachedAggregatorStart = new DiscordantPairEvidenceAggregator(source, DICTIONARY, innerWindow, outerWindow);
        cachedAggregatorStart.setCacheIntervals(cacheIntervals);
        final List<DiscordantPairEvidence> testStartCached = cachedAggregatorStart.collectEvidence(record);
        Assert.assertEquals(testStartCached.size(), expectedCount);
    }

}
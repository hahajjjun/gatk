package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointRefinerTest extends GATKBaseTest {

    private static final SAMSequenceDictionary DICTIONARY = SVTestUtils.hg38Dict;
    private static final String TEST_EVIDENCE = toolsTestDir + "/walkers/sv/printevidence/test_hg38.sr.txt.gz";
    private static final double ERROR_TOL = 1e-6;

    @DataProvider(name = "refineSplitReadSiteTestData")
    public Object[][] refineSplitReadSiteTestData() {
        return new Object[][]{
                // Empty site
                {
                    Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        Collections.emptyMap(),
                        30,
                        1,
                        new SplitReadSite(1,
                                Collections.emptyMap(),
                                null)
                },
                // 0 carrier / 1 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.emptyList(),
                        Collections.singletonList("sample1"),
                        Collections.singletonMap("sample1", 22.),
                        30,
                        1,
                        new SplitReadSite(1,
                                Collections.emptyMap(),
                                null)
                },
                // 1 carrier / 0 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.emptyList(),
                        Collections.singletonList(new HashMap.SimpleEntry<>("sample1", 22.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Collections.singletonMap("sample1", 1),
                                new EvidenceStatUtils.PoissonTestResult(0.2557291599131005, 1.3636363636363638, 0.0))
                },
                // 1 carrier / 1 background sample
                {
                        Collections.singletonList(new SplitReadEvidence("sample1", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Collections.singletonMap("sample1", 1),
                                new EvidenceStatUtils.PoissonTestResult(0.2557291599131005, 1.3636363636363638, 0.0))
                },
                // 1 carrier / 1 background sample with evidence
                {
                        Lists.newArrayList(new SplitReadEvidence("sample1", "chr21", 10, 1, true),
                                new SplitReadEvidence("sample2", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 1), new HashMap.SimpleEntry<>("sample2", 1)).stream()
                                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                                new EvidenceStatUtils.PoissonTestResult(0.8422154564080215, 1.3636363636363638, 1.5))
                },
                // 1 carrier / 1 background; multiple loci with evidence
                {
                        Lists.newArrayList(new SplitReadEvidence("sample1", "chr21", 5, 1, true),
                                new SplitReadEvidence("sample2", "chr21", 5, 1, true),
                                new SplitReadEvidence("sample1", "chr21", 10, 6, true),
                                new SplitReadEvidence("sample2", "chr21", 10, 1, true)),
                        Collections.singletonList("sample1"),
                        Collections.singletonList("sample2"),
                        Lists.newArrayList(new AbstractMap.SimpleEntry<>("sample1", 22.), new AbstractMap.SimpleEntry<>("sample2", 20.)).stream()
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                        30,
                        1,
                        new SplitReadSite(10,
                                Lists.newArrayList(new HashMap.SimpleEntry<>("sample1", 6), new HashMap.SimpleEntry<>("sample2", 1)).stream()
                                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)),
                                new EvidenceStatUtils.PoissonTestResult(0.011929713129970071, 8.181818181818182, 1.5))
                },
        };
    }

    @Test(dataProvider= "refineSplitReadSiteTestData")
    public void refineSplitReadSiteTest(final List<SplitReadEvidence> sortedEvidence,
                                        final Collection<String> carrierSamples,
                                        final Collection<String> backgroundSamples,
                                        final Map<String, Double> sampleCoverageMap,
                                        final int representativeDepth,
                                        final int defaultPosition,
                                        final SplitReadSite expected) {
        final SplitReadSite test = BreakpointRefiner.refineSplitReadSite(sortedEvidence, carrierSamples,
                backgroundSamples, sampleCoverageMap, representativeDepth, defaultPosition);
        Assert.assertEquals(test.getPosition(), expected.getPosition());
        final Set<String> samples = new HashSet(carrierSamples);
        samples.addAll(backgroundSamples);
        for (final String s : samples) {
            Assert.assertEquals(test.getCount(s), expected.getCount(s));
        }
        Assert.assertTrue((test.getP() == null && expected.getP() == null)
                || Math.abs(test.getP() - expected.getP()) <= ERROR_TOL);
        Assert.assertTrue((test.getCarrierSignal() == null && expected.getCarrierSignal() == null)
                || Math.abs(test.getCarrierSignal() - expected.getCarrierSignal()) <= ERROR_TOL);
        Assert.assertTrue((test.getBackgroundSignal() == null && expected.getBackgroundSignal() == null)
                || Math.abs(test.getBackgroundSignal() - expected.getBackgroundSignal()) <= ERROR_TOL);
    }

}
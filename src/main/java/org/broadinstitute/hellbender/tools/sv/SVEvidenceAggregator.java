package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public abstract class SVEvidenceAggregator<T extends Feature> {

    private final FeatureDataSource<T> source;
    private SimpleInterval cacheInterval;
    private Deque<T> cacheEvidence;
    private OverlapDetector<SimpleInterval> cacheIntervalTree;
    protected final SAMSequenceDictionary dictionary;

    public SVEvidenceAggregator(final FeatureDataSource<T> source,
                                final SAMSequenceDictionary dictionary) {
        Utils.nonNull(source);
        Utils.nonNull(dictionary);
        this.source = source;
        this.dictionary = dictionary;
    }

    public void setCacheIntervals(final Collection<SimpleInterval> intervals) {
        cacheIntervalTree = OverlapDetector.create(
                IntervalUtils.sortAndMergeIntervalsToStream(intervals, dictionary, IntervalMergingRule.ALL)
                        .collect(Collectors.toList())
        );
    }

    abstract public SimpleInterval getEvidenceQueryInterval(final SVCallRecord record);
    public boolean evidenceFilter(final SVCallRecord record, final T evidence) { return true; }

    private SimpleInterval getRegionInterval(final SimpleInterval interval) {
        final Set<SimpleInterval> evidenceIntervals = cacheIntervalTree.getOverlaps(interval);
        Utils.validate(evidenceIntervals.size() == 1, "Expected exactly 1 evidence interval but " +
                "found " + evidenceIntervals.size());
        return evidenceIntervals.iterator().next();
    }

    public List<T> collectEvidence(final SVCallRecord call) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getEvidenceQueryInterval(call);
        final Collection<T> rawEvidence;
        if (cacheIntervalTree == null) {
            rawEvidence = source.queryAndPrefetch(callInterval);
        } else {
            final SimpleInterval regionInterval = getRegionInterval(callInterval);
            if (!regionInterval.equals(cacheInterval)) {
                cacheEvidence = new ArrayDeque<>(source.queryAndPrefetch(regionInterval));
                cacheInterval = regionInterval;
            }
            while (!cacheEvidence.isEmpty() && callInterval.getStart() > cacheEvidence.peek().getStart()) {
                cacheEvidence.pop();
            }
            rawEvidence = cacheEvidence;
        }
        final List<T> callEvidence = new ArrayList<>();
        for (final T evidence : rawEvidence) {
            if (callInterval.overlaps(evidence)) {
                if (evidenceFilter(call, evidence)) {
                    callEvidence.add(evidence);
                }
            } else {
                break;
            }
        }
        return callEvidence;
    }
}

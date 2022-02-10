package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

public final class DepthEvidence implements SVFeature {

    final String contig;
    final int start;
    final int end;
    final int[] counts;

    public static final String BCI_VERSION = "1.0";

    public DepthEvidence(final String contig, int start, final int end, final int[] counts) {
        Utils.nonNull(contig);
        Utils.nonNull(counts);
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.counts = counts;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public int[] getCounts() { return counts; }

    @Override
    public DepthEvidence extractSamples( final List<String> sampleNames,
                                         final SVFeaturesHeader header ) {
        final List<String> headerSamples = header.getSampleNames();
        final int[] extractedCounts = new int[counts.length];
        int countsSize = 0;
        for ( final String sampleName : sampleNames ) {
            final int idx = headerSamples.indexOf(sampleName);
            if ( idx != -1 ) {
                extractedCounts[countsSize++] = counts[idx];
            }
        }
        if ( countsSize == 0 ) return null;
        final int[] justTheCounts = Arrays.copyOfRange(extractedCounts, 0, countsSize);
        return new DepthEvidence(contig, start, end, justTheCounts);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DepthEvidence)) return false;
        DepthEvidence that = (DepthEvidence) o;
        return start == that.start &&
                end == that.end &&
                contig.equals(that.contig) &&
                Arrays.equals(counts, that.counts);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(contig, start, end);
        result = 31 * result + Arrays.hashCode(counts);
        return result;
    }
}

package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.exceptions.UserException;
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
    public DepthEvidence extractSamples( final List<String> sampleNames, final Object headerObj ) {
        if ( !(headerObj instanceof SVFeaturesHeader) ) {
            throw new UserException("DepthEvidence feature source without a header.  " +
                                    "We don't know which samples we have.");
        }
        final SVFeaturesHeader header = (SVFeaturesHeader)headerObj;
        final List<String> headerSamples = header.getSampleNames();
        final int nCounts = sampleNames.size();
        final int[] newCounts = new int[nCounts];
        for ( int idx = 0; idx != nCounts; ++idx ) {
            final int countsIdx = headerSamples.indexOf(sampleNames.get(idx));
            newCounts[idx] = countsIdx == -1 ? -1 : counts[countsIdx];
        }
        return new DepthEvidence(contig, start, end, newCounts);
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

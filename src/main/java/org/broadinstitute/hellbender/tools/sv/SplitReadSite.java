package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * Container class for split read counts for multiple samples at a specific position
 */
final class SplitReadSite {
    private final int position;
    private final Map<String,Integer> sampleCountsMap;
    private final EvidenceStatUtils.PoissonTestResult result;

    /**
     * @param position breakpoint position indicated by the split reads
     * @param sampleCountsMap map with (sample id, split read count > 0) entries
     */
    public SplitReadSite(final int position, final Map<String,Integer> sampleCountsMap, final EvidenceStatUtils.PoissonTestResult result) {
        Utils.nonNull(sampleCountsMap);
        this.position = position;
        this.sampleCountsMap = sampleCountsMap;
        this.result = result;
    }

    public int getPosition() {
        return position;
    }

    public Double getP() {
        return result == null ? null : result.getP();
    }

    public Double getCarrierSignal() {
        return result == null ? null : result.getCarrierSignal();
    }

    public Double getBackgroundSignal() {
        return result == null ? null : result.getBackgroundSignal();
    }

    public int getCount(final String sample) {
        if (sampleCountsMap.containsKey(sample)) {
            return sampleCountsMap.get(sample);
        }
        return 0;
    }

}

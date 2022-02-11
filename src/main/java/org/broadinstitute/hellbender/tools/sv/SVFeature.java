package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;

import java.util.List;

public interface SVFeature extends Feature {
    SVFeature extractSamples( final List<String> sampleNames, final Object header );
}

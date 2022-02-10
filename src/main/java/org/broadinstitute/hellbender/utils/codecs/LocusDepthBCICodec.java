package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.LocusDepth;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.IOException;
import java.util.List;

public class LocusDepthBCICodec extends AbstractBCICodec<LocusDepth> {
    private boolean versionChecked = false;
    private SAMSequenceDictionary dict;
    private static final String LD_BCI_FILE_EXTENSION = ".ld.bci";

    @Override
    public LocusDepth decode( final Reader<LocusDepth> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !LocusDepth.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("bci file has wrong version: expected " +
                        LocusDepth.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        return new LocusDepth(reader.getStream(), reader.getDictionary());
    }

    @Override
    public Class<LocusDepth> getFeatureType() { return LocusDepth.class; }

    @Override
    public boolean canDecode( final String path ) {
        return path.toLowerCase().endsWith(LD_BCI_FILE_EXTENSION);
    }

    @Override
    public Writer<LocusDepth> makeSink( final GATKPath path,
                                        final SAMSequenceDictionary dict,
                                        final List<String> sampleNames,
                                        final int compressionLevel ) {
        if ( sampleNames.size() != 1 ) {
            throw new UserException("LocusDepth records do not encode their sample, and must all " +
                                    "refer to a single sample, but the list of sample names is of " +
                                    "size=" + sampleNames.size());
        }
        this.dict = dict;
        final String className = LocusDepth.class.getSimpleName();
        return new Writer<>(path,
                            new SVFeaturesHeader(className, LocusDepth.BCI_VERSION, dict, sampleNames),
                            this::encode,
                            compressionLevel);
    }

    @Override
    public void encode( final LocusDepth locusDepth, final Writer<LocusDepth> writer )
            throws IOException {
        locusDepth.write(writer.getStream(), dict);
    }
}

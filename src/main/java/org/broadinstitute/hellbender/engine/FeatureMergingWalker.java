package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;

import java.util.*;

/**
 * A FeatureMergingWalker is a tool that presents one {@link Feature} at a time in sorted order from
 * multiple sources of Features.
 *
 * To use this walker you need only implement the abstract apply method in a class that declares
 * a list of FeatureInputs as an argument.
 */
public abstract class FeatureMergingWalker<F extends Feature> extends WalkerBase {

    SAMSequenceDictionary dictionary;
    final Set<String> samples = new TreeSet<>();

    @Override
    public boolean requiresFeatures(){
        return true;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "features"; }

    @Override
    void initializeFeatures() {
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                            cloudPrefetchBuffer, cloudIndexPrefetchBuffer, getGenomicsDBOptions());
    }

    /**
     * Operations performed just prior to the start of traversal.
     */
    @Override
    public void onTraversalStart() {
        setDictionaryAndSamples();
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of Feature-based traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class
     * in the engine package that extends this class. It is not meant to be overridden by tools
     * outside the engine package.
     */
    @Override
    public void traverse() {
        final MergingIterator<F> iterator = new MergingIterator<>(dictionary, features, userIntervals);
        while ( iterator.hasNext() ) {
            final PQEntry<F> entry = iterator.next();
            final F feature = entry.getFeature();
            apply(feature, entry.getHeader());
            progressMeter.update(feature);
        }
    }

    /**
     * Process an individual feature.
     * In general, subclasses should simply stream their output from apply(), and maintain as little
     * internal state as possible.
     *
     * @param feature Current Feature being processed.
     * @param header Object for the source from which the feature was drawn (may be null)
     */
    public abstract void apply( final F feature, final Object header );

    /**
     * Get the dictionary we settled on
     */
    public SAMSequenceDictionary getDictionary() { return dictionary; }

    /**
     * Get the list of sample names we accumulated
     */
    public Set<String> getSampleNames() { return samples; }

    private void setDictionaryAndSamples() {
        dictionary = getMasterSequenceDictionary();
        if ( hasReference() ) {
            dictionary = bestDictionary(reference.getSequenceDictionary(), dictionary);
        }
        if ( hasReads() ) {
            dictionary = bestDictionary(reads.getSequenceDictionary(), dictionary);
        }
        for ( final FeatureInput<? extends Feature> input : features.getAllInputs() ) {
            final Object header = features.getHeader(input);
            if ( header instanceof SVFeaturesHeader ) {
                final SVFeaturesHeader svFeaturesHeader = (SVFeaturesHeader)header;
                dictionary = bestDictionary(svFeaturesHeader.getDictionary(), dictionary);
                final List<String> sampleNames = svFeaturesHeader.getSampleNames();
                if ( sampleNames != null ) {
                    samples.addAll(svFeaturesHeader.getSampleNames());
                }
            } else if (header instanceof VCFHeader ) {
                final VCFHeader vcfHeader = (VCFHeader)header;
                dictionary = bestDictionary(vcfHeader.getSequenceDictionary(), dictionary);
                samples.addAll(vcfHeader.getSampleNamesInOrder());
            }
        }
        if ( dictionary == null ) {
            throw new UserException("No dictionary found.  Provide one as --" +
                    StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + " or --" +
                    StandardArgumentDefinitions.REFERENCE_LONG_NAME + ".");
        }
    }

    private SAMSequenceDictionary bestDictionary( final SAMSequenceDictionary newDict,
                                                    final SAMSequenceDictionary curDict ) {
        if ( curDict == null ) return newDict;
        if ( newDict == null ) return curDict;
        final SAMSequenceDictionary smallDict;
        final SAMSequenceDictionary largeDict;
        if ( newDict.size() <= curDict.size() ) {
            smallDict = newDict;
            largeDict = curDict;
        } else {
            smallDict = curDict;
            largeDict = newDict;
        }
        int lastIdx = -1;
        for ( final SAMSequenceRecord rec : smallDict.getSequences() ) {
            final int newIdx = largeDict.getSequenceIndex(rec.getContig());
            if ( newIdx == -1 ) {
                throw new UserException("Contig " + rec.getContig() +
                                        " not found in the larger dictionary");
            }
            if ( newIdx <= lastIdx ) {
                throw new UserException("Contig " + rec.getContig() +
                                        " not in same order as in larger dictionary");
            }
        }
        return largeDict;
    }

    public static final class MergingIterator<F extends Feature> implements Iterator<PQEntry<F>> {
        final SAMSequenceDictionary dict;
        final PriorityQueue<PQEntry<F>> priorityQueue;

        @SuppressWarnings("unchecked")
        public MergingIterator( final SAMSequenceDictionary dict,
                                final FeatureManager featureManager,
                                final List<SimpleInterval> intervals ) {
            final Set<FeatureInput<? extends Feature>> inputs = featureManager.getAllInputs();
            this.dict = dict;
            this.priorityQueue = new PriorityQueue<>(inputs.size());
            for ( final FeatureInput<? extends Feature> input : inputs ) {
                final Iterator<? extends Feature> iterator =
                            featureManager.getFeatureIterator(input, intervals);
                final Object header = featureManager.getHeader(input);
                if ( iterator.hasNext() ) {
                    addEntry((Iterator<F>)iterator, header);
                }
            }
        }

        @Override
        public boolean hasNext() {
            return !priorityQueue.isEmpty();
        }

        @Override
        public PQEntry<F> next() {
            final PQEntry<F> entry = priorityQueue.poll();
            if ( entry == null ) {
                throw new NoSuchElementException("iterator is exhausted");
            }
            final Iterator<F> iterator = entry.getIterator();
            final Object header = entry.getHeader();
            if ( iterator.hasNext() ) {
                addEntry(iterator, header);
            }
            return entry;
        }

        private void addEntry( final Iterator<F> iterator, final Object header ) {
            final F feature = iterator.next();
            final int seqIdx = dict.getSequenceIndex(feature.getContig());
            if ( seqIdx == -1 ) {
                throw new UserException("dictionary has no entry for " + feature.getContig());
            }
            priorityQueue.add(new PQEntry<>(iterator, feature, seqIdx, header));
        }
    }

    public static final class PQEntry<F extends Feature> implements Comparable<PQEntry<F>> {
        private final Iterator<F> iterator;
        private final F feature;
        private final int seqIdx;
        private final Object header;

        public PQEntry( final Iterator<F> iterator, final F feature,
                        final int seqIdx, final Object header ) {
            this.iterator = iterator;
            this.feature = feature;
            this.seqIdx = seqIdx;
            this.header = header;
        }

        public Iterator<F> getIterator() { return iterator; }
        public F getFeature() { return feature; }
        public Object getHeader() { return header; }

        @Override
        public int compareTo( PQEntry<F> entry ) {
            int result = Integer.compare(seqIdx, entry.seqIdx);
            if ( result == 0 ) {
                result = Integer.compare(feature.getStart(), entry.feature.getStart());
                if ( result == 0 ) {
                    result = Integer.compare(feature.getEnd(), entry.feature.getEnd());
                }
            }
            return result;
        }
    }
}

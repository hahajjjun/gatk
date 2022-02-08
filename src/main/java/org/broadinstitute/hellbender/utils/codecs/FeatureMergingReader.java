package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.util.*;

public final class FeatureMergingReader implements FeatureReader<Feature> {
    private final SAMSequenceDictionary dict;
    private final List<FeatureReader<Feature>> readers;

    public FeatureMergingReader( final SAMSequenceDictionary dict,
                                 final List<FeatureReader<Feature>> readers ) {
        this.dict = dict;
        this.readers = readers;
    }

    @Override
    public CloseableTribbleIterator<Feature> query( String chr, int start, int end ) throws IOException {
        final List<CloseableTribbleIterator<Feature>> iterators = new ArrayList<>(readers.size());
        for ( final FeatureReader<Feature> reader : readers ) {
            iterators.add(reader.query(chr, start, end));
        }
        return new MergingIterator(dict, iterators);
    }

    @Override
    public CloseableTribbleIterator<Feature> iterator() throws IOException {
        final List<CloseableTribbleIterator<Feature>> iterators = new ArrayList<>(readers.size());
        for ( final FeatureReader<Feature> reader : readers ) {
            iterators.add(reader.iterator());
        }
        return new MergingIterator(dict, iterators);
    }

    @Override
    public void close() throws IOException {
        for ( final FeatureReader<Feature> reader : readers ) {
            reader.close();
        }
        readers.clear();
    }

    @Override
    public List<String> getSequenceNames() {
        final List<String> names = new ArrayList<>(dict.size());
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            names.add(rec.getSequenceName());
        }
        return names;
    }

    @Override public Object getHeader() {
        return null;
    }

    @Override public boolean isQueryable() {
        for ( final FeatureReader<Feature> reader : readers ) {
            if ( !reader.isQueryable() ) return false;
        }
        return true;
    }

    public static final class PQEntry implements Comparable<PQEntry> {
        private final CloseableTribbleIterator<Feature> iterator;
        private final Feature feature;
        private final int seqIdx;

        public PQEntry( final CloseableTribbleIterator<Feature> iterator,
                        final Feature feature, final int seqIdx ) {
            this.iterator = iterator;
            this.feature = feature;
            this.seqIdx = seqIdx;
        }

        public CloseableTribbleIterator<Feature> getIterator() { return iterator; }
        public Feature getFeature() { return feature; }

        @Override
        public int compareTo( PQEntry entry ) {
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

    public static final class MergingIterator implements CloseableTribbleIterator<Feature> {
        final SAMSequenceDictionary dict;
        final PriorityQueue<PQEntry> priorityQueue;

        public MergingIterator( final SAMSequenceDictionary dict,
                                final List<CloseableTribbleIterator<Feature>> iterators ) {
            this.dict = dict;
            this.priorityQueue = new PriorityQueue<>(iterators.size());
            for ( final CloseableTribbleIterator<Feature> iterator : iterators ) {
                if ( !iterator.hasNext() ) {
                    iterator.close();
                    continue;
                }
                addEntry(iterator);
            }
        }

        @Override
        public boolean hasNext() {
            return !priorityQueue.isEmpty();
        }

        @Override
        public Feature next() {
            final PQEntry entry = priorityQueue.poll();
            if ( entry == null ) {
                throw new NoSuchElementException("iterator is exhausted");
            }
            final Feature feature = entry.getFeature();
            final CloseableTribbleIterator<Feature> iterator = entry.getIterator();
            if ( iterator.hasNext() ) {
                addEntry(iterator);
            } else {
                iterator.close();
            }
            return feature;
        }

        @Override
        public void close() {
            for ( final PQEntry entry : priorityQueue ) {
                entry.getIterator().close();
            }
            priorityQueue.clear();
        }

        @Override
        public Iterator<Feature> iterator() {
            throw new GATKException("fake iterable pattern not implemented");
        }

        private void addEntry( final CloseableTribbleIterator<Feature> iterator ) {
            final Feature feature = iterator.next();
            final int seqIdx = dict.getSequenceIndex(feature.getContig());
            if ( seqIdx == -1 ) {
                throw new UserException("dictionary has no entry for " + feature.getContig());
            }
            priorityQueue.add(new PQEntry(iterator, feature, seqIdx));
        }
    }
}

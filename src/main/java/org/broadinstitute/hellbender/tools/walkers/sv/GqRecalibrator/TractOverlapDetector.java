package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableReaderOptions;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;

public class TractOverlapDetector {
    final Path tractPath;
    final SVIntervalTree<TractInterval> primaryOverlapDetector;
    final SVIntervalTree<TractInterval> otherOverlapDetector;


    public TractOverlapDetector(final String tractPath) throws IOException {
        this(Paths.get(tractPath));
    }

    public TractOverlapDetector(final Path tractPath) throws IOException {
        this.tractPath = tractPath;
        primaryOverlapDetector = new SVIntervalTree<>();
        final Iterator<TractIntervalAndOther> tractIntervalIterator = new TractTableReader(tractPath).iterator();
        final TractIntervalAndOther firstTractIntervalAndOther = tractIntervalIterator.next();
        primaryOverlapDetector.put(firstTractIntervalAndOther.primaryInterval.svInterval, firstTractIntervalAndOther.primaryInterval);
        System.out.format("Loading %s ...", tractPath);
        if(firstTractIntervalAndOther.otherInterval == null) {
            otherOverlapDetector = null;
            while(tractIntervalIterator.hasNext()) {
                final TractInterval primaryInterval = tractIntervalIterator.next().primaryInterval;
                primaryOverlapDetector.put(primaryInterval.svInterval, primaryInterval);
            }
        } else {
            otherOverlapDetector = new SVIntervalTree<>();
            otherOverlapDetector.put(firstTractIntervalAndOther.otherInterval.svInterval, firstTractIntervalAndOther.otherInterval);
            while(tractIntervalIterator.hasNext()) {
                final TractIntervalAndOther tractIntervalAndOther = tractIntervalIterator.next();
                primaryOverlapDetector.put(tractIntervalAndOther.primaryInterval.svInterval, tractIntervalAndOther.primaryInterval);
                otherOverlapDetector.put(tractIntervalAndOther.otherInterval.svInterval, tractIntervalAndOther.otherInterval);
            }
        }
        System.out.println("ok.");
    }

    public String getName() {
        final String fileName = tractPath.getFileName().toString();
        final int dotIndex = fileName.indexOf(".");
        return dotIndex > 0 ? fileName.substring(0, dotIndex) : fileName;
    }

    public boolean hasOther() { return otherOverlapDetector != null; }

    public double getPrimaryOverlapFraction(final Locatable location) {
        return getOverlapFraction(location, primaryOverlapDetector);
    }

    public double getOtherOverlapfraction(final Locatable location) {
        return getOverlapFraction(location, otherOverlapDetector);
    }

    public Set<Long> getPrimaryOverlapperIds(final Locatable location) {
        return getOverlapIds(location, primaryOverlapDetector);
    }

    public Set<Long> getOtherOverlapperIds(final Locatable location) {
        return getOverlapIds(location, otherOverlapDetector);
    }

    public boolean spansPrimaryAndOther(final Locatable location1, final Locatable location2) {
        if(otherOverlapDetector == null || location1 == null || location2 == null) {
            return false;
        }
        return !(
            setIntersect(getPrimaryOverlapperIds(location1), getOtherOverlapperIds(location2)).isEmpty()
            || setIntersect(getPrimaryOverlapperIds(location2), getOtherOverlapperIds(location1)).isEmpty()
        );
    }

    public static class UnionIterator implements Iterator<Locatable> {
        private final Iterator<SVIntervalTree.Entry<TractInterval>> tractOverlappers;
        private Locatable currentUnionLocatable;

        UnionIterator(Iterator<SVIntervalTree.Entry<TractInterval>> tractOverlappers) {
            this.tractOverlappers = tractOverlappers;
            currentUnionLocatable = null;
        }

        @Override
        public boolean hasNext() {
            return tractOverlappers.hasNext() || currentUnionLocatable != null;
        }

        @Override
        public Locatable next() {
            while(tractOverlappers.hasNext()) {
                final Locatable rawLocation = tractOverlappers.next().getValue();
                if(rawLocation == null) {
                    throw new GATKException("Got null Locatable");
                }
                if(currentUnionLocatable == null) {
                    currentUnionLocatable = rawLocation;
                } else if(rawLocation.overlaps(currentUnionLocatable)) {
                    // interval tree is sorted by start, then end. So subsequent overlappers from the tree will have a
                    // start position at least as high as currentUnionLocatable, eliminating the need for comparison
                    currentUnionLocatable = new SimpleInterval(
                            currentUnionLocatable.getContig(),
                            currentUnionLocatable.getStart(),
                            FastMath.max(currentUnionLocatable.getEnd(), rawLocation.getEnd())
                    );
                } else {
                    final Locatable unionLocatable = currentUnionLocatable;
                    currentUnionLocatable = rawLocation;
                    return unionLocatable;
                }

            }
            final Locatable unionLocatable = currentUnionLocatable;
            currentUnionLocatable = null;
            return unionLocatable;
        }
    }

    private static int getOverlapBasepairs(final Locatable loc1, final Locatable loc2) {
        return CoordMath.getOverlap(loc1.getStart(), loc1.getEnd(), loc2.getStart(), loc2.getEnd());
    }

    private static double getOverlapFraction(final Locatable location, final SVIntervalTree<TractInterval> detector) {
        if(location == null) {
            return 0.0;
        }
        final Iterator<Locatable> unionIterator = new UnionIterator(
            detector.overlappers(TractInterval.locatableToSVInterval(location))
        );
        long overlapBasePairs = 0;
        while(unionIterator.hasNext()) {
            overlapBasePairs += getOverlapBasepairs(unionIterator.next(), location);
        }
        return overlapBasePairs / (double) location.getLengthOnReference();
    }

    private static Set<Long> getOverlapIds(final Locatable location, final SVIntervalTree<TractInterval> detector) {
        final Set<Long> overlapIds = new HashSet<>();
        if(detector != null) {
            for (Iterator<SVIntervalTree.Entry<TractInterval>> it = detector.overlappers(TractInterval.locatableToSVInterval(location)); it.hasNext(); ) {
                overlapIds.add(it.next().getValue().uid);
            }
        }
        return overlapIds;
    }

    private static <T> Set<T> setIntersect(final Set<T> set1, final Set<T> set2) {
        final Set<T> intersection = new HashSet<>(set1);
        intersection.retainAll(set2);
        return intersection;
    }


    static class TractInterval implements Locatable, Comparable<Locatable> {
        final Long uid;
        final SVInterval svInterval;
        static final Map<String, Integer> contigMap = new HashMap<>();
        static final List<String> contigs = new ArrayList<>();

        TractInterval(final Long uid, final String contig, final int start, final int end) {
            this.uid = uid;
            this.svInterval = getSVInterval(contig, start, end);
        }

        public static SVInterval getSVInterval(final String contig, final int start, final int end) {
            final int contigCode;
            if(contigMap.containsKey(contig)) {
                contigCode = contigMap.get(contig);
            } else {
                contigCode = contigMap.size();
                contigMap.put(contig, contigCode);
                contigs.add(contig);
            }
            return new SVInterval(contigCode, start, end);
        }

        public static SVInterval locatableToSVInterval(final Locatable locatable) {
            return getSVInterval(locatable.getContig(), locatable.getStart(), locatable.getEnd());
        }

        public int getContigCode() { return svInterval.getContig(); }

        @Override
        public String getContig() { return contigs.get(svInterval.getContig()); }

        @Override
        public int getStart() { return svInterval.getStart(); }

        @Override
        public int getEnd() { return svInterval.getEnd(); }

        @Override
        public int compareTo(final Locatable other) {
            final int contigCompare = getContig().compareTo(other.getContig());
            if(contigCompare == 0) {
                final int startCompare = Integer.compare(getStart(), other.getStart());
                return startCompare == 0 ?
                    Integer.compare(getEnd(), other.getEnd()) :
                    startCompare;
            } else {
                return contigCompare;
            }
        }

    }

    static class TractIntervalAndOther {
        final TractInterval primaryInterval;
        final TractInterval otherInterval;

        TractIntervalAndOther(final TractInterval primaryInterval, final TractInterval otherInterval) {
            if(otherInterval != null && !primaryInterval.uid.equals(otherInterval.uid)) {
                throw new GATKException("primaryInterval.uid (" + primaryInterval.uid
                                        + ") not equal to otherInterval.uid (" + otherInterval.uid + ")");
            }
            this.primaryInterval = primaryInterval;
            this.otherInterval = otherInterval;
        }
    }

    static class TractTableReader extends TableReader<TractIntervalAndOther> {
        final static String UID_FIELD = "uid";
        final static String PRIMARY_CONTIG_FIELD = "chrom";
        final static String PRIMARY_START_FIELD = "chromStart";
        final static String PRIMARY_END_FIELD = "chromEnd";
        final static String OTHER_CONTIG_FIELD = "otherChrom";
        final static String OTHER_START_FIELD = "otherStart";
        final static String OTHER_END_FIELD = "otherEnd";
        final static TableReaderOptions defaultTableReaderOptions = new TableReaderOptions(
            true,
            new HashMap<String, String>() {
                private final static long serialVersionUID = 0;
                {
                    put("genoName", PRIMARY_CONTIG_FIELD);
                    put("genoStart", PRIMARY_START_FIELD);
                    put("genoEnd", PRIMARY_END_FIELD);
                }
            }
        );

        public TractTableReader(Path path, final TableReaderOptions tableReaderOptions) throws IOException {
            super(path, tableReaderOptions);
        }

        public TractTableReader(Path path) throws IOException {
            this(path, defaultTableReaderOptions);
        }

        @Override
        protected TractIntervalAndOther createRecord(DataLine dataLine) {
            try {
                final boolean hasOther = dataLine.columns().contains(OTHER_CONTIG_FIELD);
                final TractInterval primaryInterval = new TractInterval(
                        hasOther ? dataLine.getLong(UID_FIELD) : null, dataLine.get(PRIMARY_CONTIG_FIELD),
                        dataLine.getInt(PRIMARY_START_FIELD), dataLine.getInt(PRIMARY_END_FIELD)
                );
                final TractInterval otherInterval = hasOther ?
                        new TractInterval(
                                dataLine.getLong(UID_FIELD), dataLine.get(OTHER_CONTIG_FIELD),
                                dataLine.getInt(OTHER_START_FIELD), dataLine.getInt(OTHER_END_FIELD)
                        ) :
                        null;
                return new TractIntervalAndOther(primaryInterval, otherInterval);
            }
            catch(RuntimeException exception){
                System.out.println("columns: " + dataLine.columns().names());
                throw new RuntimeException("Error opening " + getSource() + ":\n" + exception.toString(), exception);
            }
        }
    }
}

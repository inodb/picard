/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.util.PhysicalLocation;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Class to downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair. In addition, this program uses the read-name
 * and extracts the position within the tile from whence the read came from. The downsampling is based on this position.
 * <p/>
 * Note 1: This is technology and read-name dependent. if your read-names do not have coordinate information, or if your
 * BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly.
 * This has been designed with Illumina MiSeq/HiSeq in mind.
 * <p/>
 * Note 2: The downsampling is _not_ random. It is deterministically dependent on the position of the read within its tile. specifically,
 * it draws out an ellipse that covers a PROBABILITY fraction of the area and each of the edges and uses this to determine whether to keep the
 * record. Since reads with the same name have the same position (mates, secondary and supplemental alignments), the decision will be the same for all of them.
 * <p/>
 * Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction.
 * In particular, since the reads may be distributed non-evenly within the lanes/tiles, the resulting downsampling percentage will not be accurately
 * determined by the input argument PROBABILITY. One should re-MarkDuplicates after downsampling in order to "expose" the duplicates who's representative has
 * been downsampled away.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Class to downsample a BAM file while respecting that we should either get rid\n" +
                "of both ends of a pair or neither end of the pair. In addition, this program uses the read-name \n" +
                "and extracts the position within the tile from whence the read came from. THe downsampling is based on this position. \n" +
                "results with the exact same input will produce the same results.\n" +
                "\n" +
                "<it>Note 1</it>: This is technology and read-name dependent. if your read-names do not have coordinate information, or if your\n" +
                "BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly. \n" +
                "This has been designed with Illumina MiSeq/HiSeq in mind.\n" +
                "<it>Note 2</it>: The downsampling is _not_ random. It is deterministically dependent on the position of the read within its tile." +
                "<it>Note 3</it>: Downsampling twice with this program is not supported." +
                "<it>Note 4</it>: You should MarkDuplicates after downsampling." +
                "\n" +
                "Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction. " +
                "In particular, since the reads may be distributed non-evenly within the lanes/tiles, the resulting downsampling percentage will not be accurately" +
                "determined by the input argument PROBABILITY.",
        usageShort = "Down-sample a SAM or BAM file to retain a subset of the reads based on the reads location in each tile in the flowcell.",
        programGroup = SamOrBam.class
)
public class PositionBasedDownsampleSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.PROBABILITY_SHORT_NAME, doc = "The (approximate) probability of keeping a read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

    @Option(doc = "Allow Downsampling again despite this being a bad idea with possibly unexpected results.", optional = true)
    public boolean ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS = false;

    @Option(doc = "Determines whether the duplicate tag should be reset since the downsampling requires re-marking duplicates.")
    public boolean REMOVE_DUPLICATE_INFORMATION = true;

    private final Log log = Log.getInstance(PositionBasedDownsampleSam.class);

    private PhysicalLocation opticalDuplicateFinder;
    private long total = 0;
    private long kept = 0;
    public static String PG_PROGRAM_NAME = "PositionBasedDownsampleSam";

    /* max-position in tile as a function of tile. We might need to
       look per-readgroup, but at this point I'm making the assumptions that I need to downsample a
       sample where all the readgroups came from the same type of flowcell. */

    private Map<Short, Coord> tileMaxCoord;

    final Map<Short, Histogram<Short>> xPositions = new HashMap<Short, Histogram<Short>>();
    final Map<Short, Histogram<Short>> yPositions = new HashMap<Short, Histogram<Short>>();

    public static void main(final String[] args) {
        new PositionBasedDownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        if (PROBABILITY < 0 || PROBABILITY > 1) {
            errors.add("PROBABILITY must be a value between 0 and 1, found: " + PROBABILITY);
        }

        checkProgramRecords();

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (errors.isEmpty()) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }

    @Override
    protected int doWork() {

        if (PROBABILITY < 0 || PROBABILITY > 1) {
            log.error("PROBABILITY must be a value between 0 and 1, found: " + PROBABILITY);
            throw new PicardException("Bad input.");
        }

        checkProgramRecords();

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        tileMaxCoord = new CollectionUtil.DefaultingMap<Short, Coord>(new CollectionUtil.DefaultingMap.Factory<Coord, Short>() {
            @Override
            public Coord make(final Short aShort) {
                return new Coord();
            }
        }, true);

        opticalDuplicateFinder = new PhysicalLocation();

        log.info("Starting first pass. Examining read distribution in tiles.");
        fillTileMaxCoord();
        log.info("First pass done.");

        log.info("Starting second pass. Outputting reads.");
        outputSamRecords();
        log.info("Second pass done. ");

        final double finalP = kept / (double) total;
        if (Math.abs(finalP - PROBABILITY) / (Math.min(finalP, PROBABILITY) + 1e-10) > .2) {
            log.warn(String.format("You've requested PROBABILITY=%g, the resulting downsampling resulted in a rate of %f. ", PROBABILITY, finalP));
        }
        log.info(String.format("Finished! Kept %d out of %d reads (P=%g).", kept, total, finalP));
        return 0;
    }

    private void outputSamRecords() {

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);


        final SAMFileHeader header = in.getFileHeader().clone();
        final AbstractMarkDuplicatesCommandLineProgram.PgIdGenerator pgIdGenerator = new AbstractMarkDuplicatesCommandLineProgram.PgIdGenerator(header);
        final SAMProgramRecord pr = new SAMProgramRecord(pgIdGenerator.getNonCollidingId(PG_PROGRAM_NAME));
        pr.setCommandLine(getCommandLine());
        pr.setProgramVersion(getVersion());
        header.addProgramRecord(pr);

        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        final CircleSelector selector = new CircleSelector(PROBABILITY);

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != 0 && total >= STOP_AFTER) break;

            total++;

            final PhysicalLocation pos = getSamRecordLocation(rec);

            if (!xPositions.containsKey(pos.getTile())) {
                xPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-xpos", "count"));
            }
            if (!yPositions.containsKey(pos.getTile())) {
                yPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-ypos", "count"));
            }

            final boolean keeper = selector.select(pos, tileMaxCoord.get(pos.getTile()));

            if (keeper) {
                if (REMOVE_DUPLICATE_INFORMATION) rec.setDuplicateReadFlag(false);
                out.addAlignment(rec);
                ++kept;
            }
            progress.record(rec);
        }

        CloserUtil.close(out);
        CloserUtil.close(in);
    }

    private void checkProgramRecords() {
        final SamReader in = SamReaderFactory
                .makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .open(INPUT);

        for (final SAMProgramRecord pg : in.getFileHeader().getProgramRecords()) {
            if (pg.getProgramName().equals(PG_PROGRAM_NAME)) {

                final PicardException e = new PicardException("Found previous Program Record that indicates that this BAM has been downsampled already with this program. Operation not supported! Previous PG: " + pg.toString());

                if (ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS) {
                    log.warn(e);
                } else {
                    log.error(e);
                    throw e;
                }
            }
        }
        CloserUtil.close(in);
    }

    // scan all the tiles and find the largest coordinate (x & y) in that tile.
    private void fillTileMaxCoord() {

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        int total = 0;

        for (final SAMRecord rec : in) {
            if (STOP_AFTER != 0 && total >= STOP_AFTER) break;

            total++;
            progress.record(rec);
            final PhysicalLocation location = getSamRecordLocation(rec);

            //Defaulting map will create a new Coord if it's not there.
            final Coord maxPos = tileMaxCoord.get(location.getTile());
            maxPos.x = Math.max(maxPos.x, location.getX());
            maxPos.y = Math.max(maxPos.y, location.getY());
            maxPos.count++;
        }

        // now that we know what the maximal number was, we should increase it a bit, to account for sampling error
        for (final Coord coord : tileMaxCoord.values()) {
            coord.x *= (coord.count + 1d) / coord.count;
            coord.y *= (coord.count + 1d) / coord.count;
        }
        CloserUtil.close(in);
    }

    private PhysicalLocation getSamRecordLocation(final SAMRecord rec) {
        final PhysicalLocation pos = new PhysicalLocation();
        opticalDuplicateFinder.addLocationInformation(rec.getReadName(), pos);
        return pos;
    }

    private class CircleSelector {

        private final double radiusSquared;
        private final double offset;
        private final boolean positiveSelection;

        CircleSelector(final double fraction) {
            assert (fraction >= 0);
            assert (fraction <= 1);

            final double p;
            if (fraction > 0.5) {
                p = 1 - fraction;
                positiveSelection = false;
            } else {
                p = fraction;
                positiveSelection = true;
            }
            radiusSquared = p / Math.PI; //thus the area is \pi r^2 = p, thus a fraction p of the unit square will be captured

            // if offset is used as the center of the circle (both x and y), this makes the overlap
            // region with each of the boundaries of the unit square have length p (and thus a fraction
            // p of the boundaries of each tile will be removed)
            offset = Math.sqrt(radiusSquared - p * p / 4);
        }

        private double roundedPart(final double x) {return x - Math.round(x);}

        // this function checks to see if the location of the read is within the masking circle
        private boolean select(final PhysicalLocation coord, final Coord tileMaxCoord) {
            // r^2 = (x-x_0)^2 + (y-y_0)^2, where both x_0 and y_0 equal offset
            final double distanceSquared =
                    Math.pow(roundedPart((coord.x / (double) tileMaxCoord.x) - offset), 2) +
                            Math.pow(roundedPart((coord.y / (double) tileMaxCoord.y) - offset), 2);

            return (distanceSquared > radiusSquared) ^ positiveSelection;
        }
    }

    private class Coord {
        public int x;
        public int y;
        public int count;

        public Coord() {
            count = 0;
            x = 0;
            y = 0;
        }
    }
}

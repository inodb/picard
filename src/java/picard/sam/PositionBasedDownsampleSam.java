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

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import static java.lang.Math.abs;
import static java.lang.Math.round;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair. In addition, this program uses the read-name
 * and extracts the position within the tile from whence the read came from. THe downsampling is based on this position.
 * <p/>
 * Caveat Emptor: This is technology and read-name dependent. if your read-names do not have coordinate information, or if your
 * BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly.
 * This has been designed with Illumina MiSeq/HiSeq in mind.
 * <p/>
 * Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction.
 * In particular, since the reads may be distributed non-evenly within the lanes/tiles and due to the small and finite number of swathes in
 * a tile, the resulting downsampling percentage will not be accurately
 * determined by the input argument PROBABILITY.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Class to randomly downsample a BAM file while respecting that we should either get rid\n" +
                "of both ends of a pair or neither end of the pair. In addition, this program uses the read-name \n" +
                "and extracts the position within the tile from whence the read came from. THe downsampling is based on this position. \n" +
                "results with the exact same input in the same order and with the same value for RANDOM_SEED will produce the same results.\n" +
                "\n" +
                "<it>Caveat Emptor</it>: This is technology and read-name dependent. if your read-names do not have coordinate information, or if your\n" +
                "BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly. \n" +
                "This has been designed with Illumina sequencing in mind.\n" +
                "Another caveat is that downsampling twice with this program can lead to surprising results. For example since this is effectively a " +
                "projection, downsampling twice with the same arguments will yield the same results as downsampling once. " +
                "\n" +
                "Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction. " +
                "In particular, since the reads may be distributed non-evenly within the lanes/tiles and due to the small and finite number of swathes in a tile, the resulting downsampling percentage will not be accurately" +
                "determined by the input argument PROBABILITY.",
        usageShort = "Down-sample a SAM or BAM file to retain a random subset of the reads based on the reads location in each tile in the flowcell.",
        programGroup = SamOrBam.class
)
public class PositionBasedDownsampleSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName = "MS", doc = "If not null, the output file for the metrics file.", optional = true)
    public File METRICS_OUTPUT = null;

    @Option(shortName = "RS", doc = "Random seed to use if reproducibility is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Long RANDOM_SEED = 1L;

    @Option(shortName = "P", doc = "The probability of keeping any individual swath in the tiles, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(doc = "Number of swathes into which to cut-up each tile in the x-direction. Should be between 25 and 100.", optional = true)
    public Short SWATHES_X = 31;

    @Option(doc = "Number of swathes into which to cut-up each tile in the y-direction. Should be between 25 and 100.", optional = true)
    public Short SWATHES_Y = 31;

    @Option(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

    private final Log log = Log.getInstance(PositionBasedDownsampleSam.class);

    //the resulting (approximate) probability of keeping a read after rounding due to finiteness of swathes
    private double effectiveP;

    public static void main(final String[] args) {
        new PositionBasedDownsampleSam().instanceMainWithExit(args);
    }

    private PhysicalLocation opticalDuplicateFinder;

    private long total = 0;
    private long kept = 0;
    Map<Key, Boolean> decisionMap;

    //max-position in tile as a function of tile. We might need to
    //look per-readgroup, but at this point I'm making the assumptions that I need to downsample a
    //sample where all the readgroups came from the same type of flowcell.
    private Map<Short, Coord> tileMaxCoord;
    final Map<Short, Histogram<Short>> xPositions = new HashMap<Short, Histogram<Short>>();
    final Map<Short, Histogram<Short>> yPositions = new HashMap<Short, Histogram<Short>>();

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (METRICS_OUTPUT != null) IOUtil.assertFileIsWritable(METRICS_OUTPUT);

        tileMaxCoord = new HashMap<Short, Coord>();
        opticalDuplicateFinder = new PhysicalLocation();

        //decide which swathes of all the tiles to keep
        log.info("Randomizing swathes.");
        decisionMap = makeDecisionMap();

        //if the relative difference between the requested probability and the effective predicted probability is greater than 20%,
        //issue a warning
        if (abs(effectiveP - PROBABILITY) / (Math.min(effectiveP, PROBABILITY) + 1e-10) > .2) {
            log.warn(String.format("You've requested PROBABILITY=%g, but due to the finiteness of the swathes in the tiles it looks like the probability will be %f. " +
                    "If this is a problem you can try increasing the SWATHES_X or SWATHES_Y parameters.", PROBABILITY, effectiveP));
        }

        log.debug("decision map made ");
        for (final Map.Entry<Key, Boolean> entry : decisionMap.entrySet()) {
            log.debug(String.format("(%d,%d) = %b", entry.getKey().x, entry.getKey().y, entry.getValue()));
        }

        log.info("Starting first pass. Examining read distribution in tiles.");
        fillTileMaxCoord();
        log.info("First pass done.");

        log.info("Starting second pass. Outputting reads.");
        outputRecs();
        log.info("Second pass done. ");

        if (METRICS_OUTPUT != null) {
            log.info("Outputting metrics");
            outputMetrics();
        }

        final double finalP = kept / (double) total;
        if (abs(finalP - PROBABILITY) / (Math.min(finalP, PROBABILITY) + 1e-10) > .2) {
            log.warn(String.format("You've requested PROBABILITY=%g, the resulting downsampling resulted in an rate of %f. ", PROBABILITY, finalP));
        }
        log.info(String.format("Finished! Kept %d out of %d reads (P=%g).", kept, total, finalP));
        return 0;
    }

    private void outputRecs() {

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);


        for (final SAMRecord rec : in) {
            if (STOP_AFTER != 0 && total >= STOP_AFTER) break;

            total++;

            final PhysicalLocation pos = getSamRecordLocation(rec);
            if (!xPositions.containsKey(pos.getTile()))
                xPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-xpos", "count"));
            if (!yPositions.containsKey(pos.getTile()))
                yPositions.put(pos.getTile(), new Histogram<Short>(pos.getTile() + "-ypos", "count"));

            final Key key = getKey(pos);
            xPositions.get(pos.getTile()).increment(key.x);
            yPositions.get(pos.getTile()).increment(key.x);

            if (!decisionMap.containsKey(key)) {
                final PicardException e = new PicardException("Missing Key in decision map: " + key);
                log.error(e);
                throw e;
            }

            final boolean keeper = decisionMap.get(key);

            if (keeper) {
                out.addAlignment(rec);
                ++kept;
            }
            progress.record(rec);
        }

        out.close();
        CloserUtil.close(in);
    }

    private void outputMetrics() {
        try {
            final BufferedWriter out = new BufferedWriter(new FileWriter(METRICS_OUTPUT));
            for (final Histogram<Short> histogram : xPositions.values()) {
                printHistogram(out, histogram);
                out.append("\n");
            }
            for (final Histogram<Short> histogram : yPositions.values()) {
                printHistogram(out, histogram);
                out.append("\n");
            }

            out.flush();
            out.close();

        } catch (final IOException e)
        {
            log.error("error while writing to metrics file");
            throw new PicardException("unknown error");
        }
    }

    private void printHistogram(final BufferedWriter out, final Histogram<Short> histogram) throws IOException {
        final String SEPARATOR = "\t";

        final FormatUtil formatter = new FormatUtil();
        // Output a header row
        out.append(StringUtil.assertCharactersNotInString(histogram.getBinLabel(), '\t', '\n'));
        out.append(SEPARATOR);
        out.append(StringUtil.assertCharactersNotInString(histogram.getValueLabel(), '\t', '\n'));
        out.newLine();

        for (final Short key : histogram.keySet()) {
            out.append(key.toString());
            final Histogram<Short>.Bin bin = histogram.get(key);
            final double value = (bin == null ? 0 : bin.getValue());

            out.append(SEPARATOR);
            out.append(formatter.format(value));
            out.newLine();
        }
    }

    private Map<Key, Boolean> makeDecisionMap() {
        final Map<Key, Boolean> decisionMap = new HashMap<Key, Boolean>();
        final Random r = RANDOM_SEED == null ? new Random() : new Random(RANDOM_SEED);

        int availableSwathes = SWATHES_X * SWATHES_Y;
        int keptSwathes = 0;

        final int requestedSwathes = (int) round(SWATHES_X * SWATHES_Y * PROBABILITY);

        // non-standard looping is on purpose since the last element is set to equal the first.
        for (short i = 0; i < SWATHES_X - 1; i++) {
            // non-standard looping is on purpose since the last element is set to equal the first.
            for (short j = 0; j < SWATHES_Y - 1; j++) {
                //we are aiming to keep exactly requestedSwathes. So the probability of keeping every
                //swath depends on the number already kept (and the amount remaining to choose from).
                //based on discussion in http://stackoverflow.com/questions/48087/

                //due to the fact that some swathes count for 2 or 4 (edges and corners) it might fail for small values of PROBABILITY

                final Boolean keepQ = r.nextDouble() <= (requestedSwathes - keptSwathes) / availableSwathes;

                decisionMap.put(new Key(i, j), keepQ);
                if (keepQ) keptSwathes++;
                availableSwathes--;
                //to emulate neighboring tiles having overlapping regions that can cause duplicates, we make sure that the
                //swaths kept on the left/top are the same as the right/bottom
                if (i == 0) {
                    decisionMap.put(new Key((short) (SWATHES_X - 1), j), keepQ);
                    if (keepQ) keptSwathes++;
                    availableSwathes--;
                    if (j == 0) {
                        decisionMap.put(new Key((short) (SWATHES_X - 1), (short) (SWATHES_Y - 1)), keepQ);
                        if (keepQ) keptSwathes++;
                        availableSwathes--;
                    }
                }
                if (j == 0) {
                    decisionMap.put(new Key(i, (short) (SWATHES_Y - 1)), keepQ);
                    if (keepQ) keptSwathes++;
                    availableSwathes--;
                }
            }
        }

        effectiveP = keptSwathes / (double) decisionMap.size();
        log.info(String.format("Swath mask has %d entries, of which %d will be kept. The effective downsampling ratio is %g", decisionMap.size(), keptSwathes, effectiveP));
        return decisionMap;
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

            if (!tileMaxCoord.containsKey(location.getTile())) {
                tileMaxCoord.put(location.getTile(), new Coord(location.getX(), location.getY()));

            } else {
                final Coord maxPos = tileMaxCoord.get(location.getTile());
                maxPos.x = Math.max(maxPos.x, location.getX());
                maxPos.y = Math.max(maxPos.y, location.getY());
                maxPos.count++;
            }
        }

        // now that we know what the maximal number was, we should increase it a bit, to account for sampling error
        for (final Coord coord : tileMaxCoord.values()) {
            coord.x *= (coord.count + 1d) / coord.count;
            coord.y *= (coord.count + 1d) / coord.count;
        }
        CloserUtil.close(in);

    }

    private Key getKey(final PhysicalLocation pos) {
        final Coord maxCoord = tileMaxCoord.get(pos.getTile());
        final short keyX = (short) Math.min(SWATHES_X - 1, pos.getX() * SWATHES_X / maxCoord.x);
        final short keyY = (short) Math.min(SWATHES_Y - 1, pos.getY() * SWATHES_Y / maxCoord.y);

        return new Key(keyX, keyY);
    }

    private PhysicalLocation getSamRecordLocation(final SAMRecord rec) {
        final PhysicalLocation pos = new PhysicalLocation();
        opticalDuplicateFinder.addLocationInformation(rec.getReadName(), pos);
        return pos;
    }

    static private class Coord {
        public int x;
        public int y;
        public int count;

        public Coord(final int x, final int y) {
            this.x = x;
            this.y = y;
            this.count = 1;
        }
    }

    static private class Key {

        private final short x;
        private final short y;

        public Key(final short x, final short y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (!(o instanceof Key)) return false;
            final Key key = (Key) o;
            return x == key.x && y == key.y;
        }

        @Override
        public int hashCode() {
            int result = x;
            result = 31 * result + y;
            return result;
        }

        @Override
        public String toString() {
            return String.format("(%d, %d)", x, y);
        }
    }
}

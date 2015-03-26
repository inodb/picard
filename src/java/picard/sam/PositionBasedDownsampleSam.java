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
 */
@CommandLineProgramProperties(
        usage = "Class to randomly downsample a BAM file while respecting that we should either get rid\n" +
                "of both ends of a pair or neither end of the pair. In addition, this program uses the read-name \n" +
                "and extracts the position within the tile from whence the read came from. THe downsampling is based on this position. \n" +
                "results with the exact same input in the same order and with the same value for RANDOM_SEED will produce the same results.\n" +
                "\n" +
                "Caveat Emptor: This is technology and read-name dependent. if your read-names do not have coordinate information, or if your\n" +
                "BAM contains reads from multiple technologies (flowcell versions, sequencing machines) this will not work properly. \n" +
                "This has been designed with Illumina MiSeq/HiSeq in mind.\n" +
                "\n" +
                "Finally, the code has been designed to simulate sequencing less as accurately as possible, not for getting an exact downsample fraction.  ",
        usageShort = "Down-sample a SAM or BAM file to retain a random subset of the reads based on the reads location in each tile in the flowcell",
        programGroup = SamOrBam.class
)
public class PositionBasedDownsampleSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName = "MS", doc = "If not null, the output file for the metrics file.")
    public File METRICS_OUTPUT;

    @Option(shortName = "RS", doc = "Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Long RANDOM_SEED = 1L;

    @Option(shortName = "P", doc = "The probability of keeping any individual swath in the tiles, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(doc = "Number of swathes into which to cut-up each tile in the x-direction.", optional = true)
    public Short SWATHES_X = 30;

    @Option(doc = "Number of swathes into which to cut-up each tile in the y-direction.", optional = true)
    public Short SWATHES_Y = 30;

    @Option(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

    private final Log log = Log.getInstance(PositionBasedDownsampleSam.class);

    public static void main(final String[] args) {
        new PositionBasedDownsampleSam().instanceMainWithExit(args);
    }

    private PhysicalLocation opticalDuplicateFinder;

    //max-position in tile as a function of tile. we might need to
    //look per-readgroup, but at this point I'm making the assumptions that I need to downsample a
    //sample where all the readgroups came from the same type of flowcell.
    private Map<Short, Coord> tileMaxCoord;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (METRICS_OUTPUT != null) IOUtil.assertFileIsWritable(METRICS_OUTPUT);

        tileMaxCoord = new HashMap<Short, Coord>();
        opticalDuplicateFinder = new PhysicalLocation();


        log.info("Randomizing swathes.");
        final Map<Key, Boolean> decisionMap = makeDecisionMap();

        log.debug("decision map made ");
        for(final Map.Entry<Key,Boolean> entry : decisionMap.entrySet()){
            log.debug(String.format("(%d,%d) = %b",entry.getKey().x,entry.getKey().y,entry.getValue()));
        }

        fillTileMaxCoord();

        log.info("Starting Second pass.");

        long total = 0;
        long kept = 0;

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);

        final Map<Short, Histogram<Short>> xPositions = new HashMap<Short, Histogram<Short>>();
        final Map<Short, Histogram<Short>> yPositions = new HashMap<Short, Histogram<Short>>();

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

            if(!decisionMap.containsKey(key)){
                final PicardException e = new PicardException("Missing Key in decision map: "+key);
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

        if (METRICS_OUTPUT != null) {
            try {
                final BufferedWriter mout = new BufferedWriter(new FileWriter(METRICS_OUTPUT));
                for (final Histogram<Short> histogram : xPositions.values()) {
                    printHistogram(mout, histogram);
                    mout.append("\n");
                }
                for (final Histogram<Short> histogram : yPositions.values()) {
                    printHistogram(mout, histogram);
                    mout.append("\n");
                }

                mout.flush();
                mout.close();

            } catch (final IOException e) {
                log.error("error while writing to metrics file");
                throw new PicardException("unknown error");
            }
        }

        log.info("Finished! Kept " + kept + " out of " + total + " reads.");

        return 0;
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

        // non-standard looping is on purpose since the last element is set to equal the first.
        for (short i = 0; i < SWATHES_X - 1; i++) {
            // non-standard looping is on purpose since the last element is set to equal the first.
            for (short j = 0; j < SWATHES_Y - 1; j++) {
                final Boolean keepQ = r.nextDouble() <= PROBABILITY;

                decisionMap.put(new Key(i, j), keepQ);
                //to emulate neighboring tiles having overlapping regions that can cause duplicates, we make sure that the
                //swaths kept on the left/top are the same as the right/bottom
                if (i == 0) {
                    decisionMap.put(new Key((short) (SWATHES_X - 1), j), keepQ);
                    if (j == 0) {
                        decisionMap.put(new Key((short) (SWATHES_X - 1), (short) (SWATHES_Y - 1)), keepQ);
                    }
                }
                if (j == 0) {
                    decisionMap.put(new Key(i, (short) (SWATHES_Y - 1)), keepQ);
                }
            }
        }
        int totalKeeps=0;
        for(final Boolean keep: decisionMap.values()){
            if(keep) totalKeeps++;
        }

        log.info(String.format("Swath mask has %d entries, of which %d will be kept. The effective downsampling ratio is %g", decisionMap.size(),totalKeeps, totalKeeps/(double)decisionMap.size()));
        return decisionMap;
    }

    private void fillTileMaxCoord() {

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        log.info("first pass. examining read distribution in tiles.");

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
        //TODO: check math logic here.
        for (final Coord coord : tileMaxCoord.values()) {
            coord.x *= (coord.count + 1d)/coord.count ;
            coord.y *= (coord.count + 1d)/coord.count ;
        }
        CloserUtil.close(in);

        log.info("first pass Done.");
    }

    private Key getKey(final PhysicalLocation pos) {
        final Coord maxCoord = tileMaxCoord.get(pos.getTile());
        final short keyX = (short) Math.min(SWATHES_X-1, pos.getX() * SWATHES_X / maxCoord.x);
        final short keyY = (short) Math.min(SWATHES_Y-1, pos.getY() * SWATHES_Y / maxCoord.y);

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
            return String.format("(%d, %d)",x , y);
        }
    }
}

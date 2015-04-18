package picard.sam;

import edu.mit.broad.picard.util.TestNGUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class PositionBasedDownsampleSamTest extends CommandLineProgramTest {
    final static String TEST_DIR = "testdata/picard/sam/PositionalDownsampleSam/";

    @Override
    public String getCommandLineProgramName() {
        return PositionBasedDownsampleSam.class.getSimpleName();
    }

    // consider using SAMRecordSetBuilder
    SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

    @DataProvider(name = "ValidArgumentsTestProvider")
    public Object[][] ValidArgumentsTestProvider() {
        final List<Object[]> objects = new ArrayList<Object[]>();
        for (double i = 0; i <= 1; i += .1) {
            final Object[] array = {i};
            objects.add(array);
        }
        return (Object[][]) objects.toArray();
    }

    // test removing some reads from a sparse, single tile bam
    @Test(dataProvider = "ValidArgumentsTestProvider")
    public void testDownsampleSingleTile(final double fraction) throws IOException {
        testDownsampleWorker(new File("goodSingleTile.bam"), fraction);
    }

    // test removing some reads from a sparse, multiple-tile bam
    @Test(dataProvider = "ValidArgumentsTestProvider")
    public void testDownsampleLane(final double fraction) throws IOException {
        testDownsampleWorker(new File("goodLane.bam"), fraction);
    }

    public void testDownsampleWorker(final File samFile, final double fraction) throws IOException {

        final File downsampled = File.createTempFile(samFile + ".", ".sam");
        downsampled.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampled.getAbsolutePath(),
                "PROBABILITY=" + fraction
        };

        // make sure results is successful
        assertEquals(runPicardCommandLine(args), 0);

        // make sure that the resulting BAM is valid.
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        validateSamFile.INPUT = downsampled;
        assertEquals(validateSamFile.doWork(), 0);

        //make sure that the total number of record in the resulting file in in the ballpark:
        TestNGUtil.assertGreaterThan(countSamTotalRecord(downsampled), fraction * .9 * countSamTotalRecord(samFile));
        TestNGUtil.assertLessThan(countSamTotalRecord(downsampled), fraction * 1.1 * countSamTotalRecord(samFile));
    }

    private long countSamTotalRecord(final File samFile) {
        final SamReader reader = SamReaderFactory.make().open(samFile);
        assert reader.hasIndex();
        long total = 0;

        for (int i = 0; i < reader.getFileHeader().getSequenceDictionary().size(); i++) {
            total += reader.indexing().getIndex().getMetaData(i).getAlignedRecordCount();
            total += reader.indexing().getIndex().getMetaData(i).getUnalignedRecordCount();
        }
        return total;
    }

    // test that downsampling twice yields a warning.
    @Test(expectedExceptions = PicardException.class)
    public void Test() throws IOException {
        final File samFile = new File(TEST_DIR, "good.bam");
        final File downsampledOnce = File.createTempFile(samFile + ".", ".1.sam");
        final double fraction = .1;
        downsampledOnce.deleteOnExit();
        final String[] args1 = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampledOnce.getAbsolutePath(),
                "PROBABILITY=" + fraction
        };
        runPicardCommandLine(args1);

        final File downsampledTwice = File.createTempFile(samFile + ".", ".2.sam");
        downsampledOnce.deleteOnExit();
        final String[] args2 = new String[]{
                "INPUT=" + downsampledOnce.getAbsolutePath(),
                "OUTPUT=" + downsampledTwice.getAbsolutePath(),
                "PROBABILITY=" + fraction
        };
        //should blow up due to bad inputs
        runPicardCommandLine(args2);
    }

    // test that program fails on p<0  or p>1
    @DataProvider(name = "InvalidArgumentsTestProvider")
    public Object[][] InvalidArgumentsTestProvider() {
        return new Object[][]{{-1.0}, {-.00001}, {-5.0}, {1.00001}, {5.0}, {50.0}, {Double.MAX_VALUE}, {Double.POSITIVE_INFINITY}, {Double.NEGATIVE_INFINITY}};
    }

    @Test(dataProvider = "InvalidArgumentsTestProvider", expectedExceptions = PicardException.class)
    public void Test(final double fraction) throws IOException {
        final File samFile = new File(TEST_DIR, "good.bam");
        final File downsampledFile = File.createTempFile(samFile + ".", ".sam");
        downsampledFile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampledFile.getAbsolutePath(),
                "PROBABILITY=" + fraction
        };
        //should blow up due to bad inputs
        runPicardCommandLine(args);
    }
}
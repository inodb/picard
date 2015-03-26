/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import picard.PicardException;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Contains class for figuring out the location of reads.
 *
 * @author Tim Fennell
 * @author Nils Homer
 * @author Yossi Farjoun
 */

/**
 * Small interface that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
 * non-zero positive integers, x and y coordinates may be negative.
 */
public class PhysicalLocation {

    public static final String DEFAULT_READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*".intern();

    private final String readNameRegex;

    public PhysicalLocation() {this(DEFAULT_READ_NAME_REGEX);}

    public PhysicalLocation(final String readNameRegExp) {this.readNameRegex = readNameRegExp;}

    private Pattern readNamePattern;

    public short tile = -1;
    public int x = -1, y = -1;


    public short getTile() { return this.tile; }

    public void setTile(final short tile) { this.tile = tile; }

    public int getX() { return this.x; }

    public void setX(final int x) { this.x = x; }

    public int getY() { return this.y; }

    public void setY(final int y) { this.y = y; }


    private final int[] tmpLocationFields = new int[10]; // for optimization of addLocationInformation

    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocation so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc      the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    public boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        // Optimized version if using the default read name regex (== used on purpose):
        if (this.readNameRegex == this.DEFAULT_READ_NAME_REGEX) {
            final int fields = getRapidDefaultReadNameRegexSplit(readName, ':', tmpLocationFields);
            if (!(fields == 5 || fields == 7)) {
                throw new PicardException(String.format("Default READ_NAME_REGEX '%s' did not match read name '%s'.  " +
                                "You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates.  " +
                                "Note that this message will not be emitted again even if other read names do not match the regex.",
                        this.readNameRegex, readName));
            }

            final int offset = fields == 7 ? 2 : 0;
            loc.setTile((short) tmpLocationFields[offset + 2]);
            loc.setX(tmpLocationFields[offset + 3]);
            loc.setY(tmpLocationFields[offset + 4]);
            return true;
        } else if (this.readNameRegex == null) {
            return false;
        } else {
            // Standard version that will use the regex
            if (this.readNamePattern == null) this.readNamePattern = Pattern.compile(this.readNameRegex);

            final Matcher m = this.readNamePattern.matcher(readName);
            if (m.matches()) {
                loc.setTile((short) Integer.parseInt(m.group(1)));
                loc.setX(Integer.parseInt(m.group(2)));
                loc.setY(Integer.parseInt(m.group(3)));
                return true;
            } else {
                throw new PicardException(String.format("READ_NAME_REGEX '%s' did not match read name '%s'.  Your regex may not be correct.  " +
                                "Note that this message will not be emitted again even if other read names do not match the regex.",
                        this.readNameRegex, readName));
            }
        }
    }

    //TODO refactor this. code duplication within OpticalDuplicateFinder

    /**
     * Single pass method to parse the read name for the default regex.  This will only insert the 2nd to the 4th
     * tokens (inclusive).  It will also stop after the fifth token has been successfully parsed.
     */
    static protected int getRapidDefaultReadNameRegexSplit(final String readName, final char delim, final int[] tokens) {
        int tokensIdx = 0;
        int prevIdx = 0;
        for (int i = 0; i < readName.length(); i++) {
            if (readName.charAt(i) == delim) {
                if (1 < tokensIdx && tokensIdx < 5)
                    tokens[tokensIdx] = rapidParseInt(readName.substring(prevIdx, i)); // only fill in 2-4 inclusive
                tokensIdx++;
                if (4 < tokensIdx) return tokensIdx; // early return, only consider the first five tokens
                prevIdx = i + 1;
            }
        }
        if (prevIdx < readName.length()) {
            if (1 < tokensIdx && tokensIdx < 5)
                tokens[tokensIdx] = rapidParseInt(readName.substring(prevIdx, readName.length())); // only fill in 2-4 inclusive
            tokensIdx++;
        }
        return tokensIdx;
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     */
    static protected int rapidParseInt(final String input) {
        final int len = input.length();
        int val = 0;
        int i = 0;
        boolean isNegative = false;

        if (0 < len && '-' == input.charAt(0)) {
            i = 1;
            isNegative = true;
        }

        for (; i < len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val * 10) + (ch - 48);
            } else {
                break;
            }
        }

        if (isNegative) val = -val;

        return val;
    }
}

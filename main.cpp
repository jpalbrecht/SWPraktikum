#include <omp.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/arg_parse.h>


#include <seqan/bam_io.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/stream.h>


using namespace seqan;


/* Struct. Defines Parsing Options and CommandLine Arguments.
 *
 * Properties:
 *   std::string inRef      - Path to ReferenceFile. Can be a Fasta or Index File. If Index File please set Option -i
 *   std::string inReads    - Path to ReadsFile.
 *   std::string outSam     - Path to Output File
 *   bool loadIndex         - Boolean indicating if inRef is Fasta or Index File
 *   int seedLen            - declare seed Length to split Reads
 *   int indexMode          - IndexMode. If 0 program will use SuffixArray, if set to 1 program will use FM-Index
 *   bool levenshtein       - Boolean indicating Distance to calculate seedExtension with. True will result in the use
 *                                  of Levenshtein Distance. False to Hamming-Distance
 *   Score<int, Simple> alignScheme  - AlignScheme for Alignment between Read and ReferenceGenome
 *
 * Functions:
 *   ReadMapperOptions  -  Constructor of Struct
 *
 * Other header files required: <seqan/arg_parse.h> <seqan/basic.h>
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
struct ReadMapperOptions {
    std::string inRef;
    std::string inReads;
    std::string outSam;
    bool loadIndex;
    bool levenshtein;
    bool fmInd;
    bool buildIndexOnly;
    int seedLen;
    int dropout;
    int threads;
    Score<int, Simple> alignScheme;
    Score<int, Simple> extendScheme;


    ReadMapperOptions() :
            inRef(""), inReads(""), outSam(""), loadIndex(false), fmInd(false), seedLen(10),
            levenshtein(false), alignScheme(1, -1, 0, -1), dropout(3), extendScheme(2, -1, -1, -2), threads(1) { }

};


/* Struct to set Distance Object in runtime */
struct Levenshtein {
    Levenshtein() { }
};


/* Struct to set Distance Object in runtime */
struct Hamming {
    Hamming() { }
};


/* Function to call the extendSeed-Function with to runtime specialized options
 *
 * Syntax:
 *   extendWrapper(Seed<Simple> & seed1 , Dna5String &refGen,
 *                               ReadMapperOptions &rmOptions, Levenshtein &lv )
 *
 * Inputs:
 *   Seed<Simple> & seed1           - Seed to extend
 *   Dna5String &refGen             - reference Gen to extend seed over
 *   TText &read                    - read to extend Seed over
 *   ReadMapperOptions &rmOptions   - Options containing AlignmentScheme and dropOut-threshold
 *   Levenshtein                    - object telling the compiler wich version of extendSeed to call
 *
 * Outputs:
 *
 * Example:
 *   extendWrapper(seed1, refGen, rmOptions, lv)
 *
 * Other header files required: none
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
template<typename TText>
inline static void extendWrapper(Seed<Simple> &seed1, Dna5String &refGen, TText &read,
                                 ReadMapperOptions &rmOptions, Levenshtein) {
    extendSeed(seed1, refGen, read, EXTEND_BOTH, rmOptions.extendScheme, 3, GappedXDrop());
}


/* Function to call the extendSeed-Function with to runtime specialized options
 *
 * Syntax:
 *   extendWrapper(Seed<Simple> & seed1 , Dna5String &refGen,
 *                               ReadMapperOptions &rmOptions, Hamming &ha )
 *
 * Inputs:
 *   Seed<Simple> & seed1           - Seed to extend
 *   Dna5String &refGen             - reference Gen to extend seed over
 *   TText &read                    - read to extend Seed over
 *   ReadMapperOptions &rmOptions   - Options containing AlignmentScheme and dropOut-threshold
 *   Hamming                        - object telling the compiler which version of extendSeed to call
 *
 * Outputs:
 *
 * Example:
 *   extendWrapper(seed1, refGen, rmOptions, ha)
 *
 * Other header files required: none
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
template<typename TText>
inline static void extendWrapper(Seed<Simple> &seed1, Dna5String &refGen, TText &read,
                                 ReadMapperOptions &rmOptions, Hamming) {
    extendSeed(seed1, refGen, read, EXTEND_BOTH, MatchExtend());
}


/* Function to parse Arguments from Command Line
 *
 * Syntax:
 *   ArgumentParser::ParseResult parseCommandLine(ReadMapperOptions &rmOptions, int argc, char const **argv)
 *
 * Inputs:
 *   ReadMapperOptions &rmOptions   - optionObject from ReadMapperOption Class.
 *                                      See Struct Documentation for more information
 *   int argc                       - Number of Arguments
 *   char const **argv              - Arguments from Command Line
 *
 * Outputs:
 *  ArgumentParser::ParseResult     - Result of parsing Arguments. Can be PARSE_OK or PARSE_ERROR.
 *
 * Example:
 *   parseCommandLine(rmOptions, 2 , argv )
 *
 * Other header files required: <seqan/arg_parse.h>, <seqan/basic.h>
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
inline ArgumentParser::ParseResult parseCommandLine(ReadMapperOptions &rmOptions, int argc, char const **argv) {
    //------------ Set Up Parser and Arguments
    ArgumentParser parser("readMapper");
    // Set short description, version, and date.
    setShortDescription(parser, "Read Mapper");
    setVersion(parser, "1.0");
    setDate(parser, "Mai 2016");
    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIInputReference\\fP] [ \"\\fIOPTIONS\\fP\"  ]");
    addDescription(parser, "This program allows to map Reads to a Reference Genome.");
    // 1 arguments required
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "InPathRef"));



    //------------ Define Options -- Section Modification Options
    addSection(parser, "Modification Options");
    addOption(parser, ArgParseOption("i", "index", "Option to load ReferenceGenome out of IndexFile."));
    addOption(parser, ArgParseOption("b", "buildIndex", "Option to just build the Index File"));
    addOption(parser, ArgParseOption("r", "readFile", "set Path to File containing Reads. FASTA or FASTQ",
                                     ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption("o", "outSam", "set Path to output File",
                                     ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption("t", "threads", "set Number of threads to use",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("s", "seedLength", "Number of Characters of seed",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("S", "suffixArray", "set the Index Modus to SuffixArray"));
    addOption(parser, ArgParseOption("F", "fmIndex", "set the Index Modus to FM-Index"));
    addOption(parser, ArgParseOption("l", "levenshtein", "set the extension Distance to Levenshtein-Distance"));
    addOption(parser, ArgParseOption("d", "dropout", "set the dropout threshold when to stop extend seeds",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ha", "hamming", "setting the extension Distance to Hamming-Distance"));
    addOption(parser, ArgParseOption("ma", "match", "MatchCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mm", "mismatch", "MisMatchCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("go", "gapOpen", "GapOpenCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ge", "gapExtend", "GapExtendCost",
                                     ArgParseArgument::INTEGER, "INT"));


    //------------ Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-i\\fP ",
                "Load ReferenceGenome out of IndexFile");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP  \\fB-s\\fP \\fI10\\fP ",
                "Set seedLength to 10(default)");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-b\\fP \\fB-S\\fP ",
                "Set Index Modus to SuffixArray. Build Suffix Array only-Modus");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-b\\fP \\fB-F\\fP ",
                "Set Index Modus to FM-Index. Build FM-Index only-Modus");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-t\\fP \\fB4\\fP",
                "Set the extendSeed Mode to Levenshtein-Distance. Use 4 threads for searching");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-d\\fP \\fI3\\fP",
                "Set the droput threshold to 3 (default value).");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-ha\\fP ",
                "Set the extendSeed Mode to Hamming-Distance(default)");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-ma\\fP \\fI2\\fP ",
                "Set MatchCost to 2");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-mm\\fP \\fI-1\\fP ",
                "Set MisMatchCost to -1");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-go\\fP \\fI-2\\fP ",
                "Set GapOpenCost to -2");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-ge\\fP \\fI-1\\fP ",
                "Set GapExtendCost to -1");


    //------------ Parse command line & read necessary Arguments and optional Values
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;
    // read Values
    getArgumentValue(rmOptions.inRef, parser, 0);
    getOptionValue(rmOptions.inReads, parser, "readFile");
    getOptionValue(rmOptions.outSam, parser, "outSam");
    rmOptions.loadIndex = isSet(parser, "index");
    getOptionValue(rmOptions.seedLen, parser, "seedLength");
    rmOptions.fmInd = isSet(parser, "fmIndex");
    rmOptions.levenshtein = isSet(parser, "levenshtein");
    getOptionValue(rmOptions.dropout, parser, "dropout");
    getOptionValue(rmOptions.threads, parser, "threads");
    rmOptions.buildIndexOnly = isSet(parser, "buildIndex");

    //------------ read & set ScoringScheme if specified
    int ma, mm, go, ge;
    getOptionValue(ma, parser, "match");
    getOptionValue(mm, parser, "mismatch");
    getOptionValue(go, parser, "gapOpen");
    getOptionValue(ge, parser, "gapExtend");
    if (isSet(parser, "match")) {
        setScoreMatch(rmOptions.alignScheme, ma);
    }
    if (isSet(parser, "mismatch")) {
        setScoreMismatch(rmOptions.alignScheme, mm);
    }
    if (isSet(parser, "gapOpen")) {
        setScoreGapOpen(rmOptions.alignScheme, go);
    }
    if (isSet(parser, "gapExtend")) {
        setScoreGapExtend(rmOptions.alignScheme, ge);
    }


    //------------ check dissent options & return
    if (!isSet(parser, "buildIndex")) {
        if (!isSet(parser, "readFile")) {
            std::cerr << "ERROR: Please specify Read-File with option r!\n";
            return seqan::ArgumentParser::PARSE_ERROR;
        }
    }
    if (isSet(parser, "suffixArray") && isSet(parser, "fmIndex")) {
        std::cerr << "ERROR: You cannot specify both suffixArray- and fmIndex-Modus!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    if (isSet(parser, "levenshtein") && isSet(parser, "hamming")) {
        std::cerr << "ERROR: You cannot specify both levenshtein- and hamming-Distance!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    if (isSet(parser, "dropout") && isSet(parser, "hamming")) {
        std::cerr << "ERROR: You cannot specify a dropout with hamming-Distance!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    if (isSet(parser, "dropout") && !isSet(parser, "levenshtein")) {
        std::cerr << "ERROR: You cannot specify a dropout without levenshtein-Distance!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    if (isSet(parser, "buildIndex")) {
        if (isSet(parser, "readFile") || isSet(parser, "levenshtein") || isSet(parser, "hamming")
            || isSet(parser, "gapExtend") || isSet(parser, "gapOpen") || isSet(parser, "mismatch")
            || isSet(parser, "match") || isSet(parser, "threads") || isSet(parser, "dropout")
            || isSet(parser, "seedLength") || isSet(parser, "index") || isSet(parser, "outSam")) {
            std::cout <<
            "WARNING: Options specified although program just creates Index-File! Input will be ignored!\n";
        }
    }
    if (!isSet(parser, "buildIndex")) {
        if (!isSet(parser, "outSam")) {
            std::cout << "WARNING: Out-File not specified! Sam-File will be written in same directory as Input-File!\n";
            rmOptions.outSam = rmOptions.inReads;
            rmOptions.outSam += ".out";
        }
    }
    return ArgumentParser::PARSE_OK;
};


/* Function to read a FastaFile in & building/saving Index-Structure of content in file.
 * Returns true if file could be read and index could be built/wrote to disk, otherwise false
 *
 * Syntax:
 *   int erg = buildIndex(CharString fileIn, CharString fileOut);
 *
 * Inputs:
 *   CharString fileIn              - complete Path to fastA file
 *   CharString fileOut             - complete Path to output index file
 *   ReadMapperOptions &rmOptions   - read mapper program option object
 *
 * Outputs:
 *  int     - true if Output-file written, false otherwise
 *
 * Example:
 *   erg = buildIndex(fileIn, fileOut);
 *
 * Other header files required: <seqan/index.h>, <seqan/basic.h>
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
static inline int buildIndex(std::string &fileIn, std::string &fileOut, ReadMapperOptions &rmOptions) {
    //------------ Read in FastaReferenceFile & print Info
    SeqFileIn sFileIn;
    std::cout << "Reading Reference-File.." << std::endl;
    if (!open(sFileIn, toCString(fileIn))) {
        std::cerr << "Could not read Reference-File!" << std::endl;
        return 0;
    }
    std::cout << "File read!" << std::endl;
    // read Reference-Sequence
    Dna5String genome;
    CharString head;
    try {
        readRecord(head, genome, sFileIn);
    } catch (IOError io) {
        std::cerr << "ERROR: " << io.what() << std::endl;
    } catch (ParseError p) {
        std::cerr << "ERROR: " << p.what() << std::endl;
    }


    //------------ build Index & save to disk & print info
    std::cout << "Building Index.." << std::endl;
    // build Index
    if (!rmOptions.fmInd) {
        Index<Dna5String, IndexSa<> > index(genome);
        // require Index for being able to save to disk
        indexRequire(index, FibreSA());
        std::cout << "Index built!" << std::endl;
        // save Index to File!
        std::cout << "Saving Index to disk.." << std::endl;
        if (!save(index, toCString(fileOut))) {
            std::cerr << "ERROR: not able to write Index-File to disk. Check for available space & authorization" <<
            std::endl;
            return 0;
        }
        std::cout << "Index saved!" << std::endl;
    } else {
        Index<Dna5String, FMIndex<> > index(genome);
        // require Index for being able to save to disk
        indexRequire(index, FibreSALF());
        std::cout << "Index built!" << std::endl;
        // save Index to File!
        std::cout << "Saving Index to disk.." << std::endl;
        if (!save(index, toCString(fileOut))) {
            std::cerr << "ERROR: not able to write Index-File to disk. Check for available space & authorization" <<
            std::endl;
            return 0;
        }
        std::cout << "Index saved!" << std::endl;
    }
    return 1;
};


/* Function to get the CIGAR-Format out of an Alignment-Object.
 *  First Gaps and Last Gaps will be ignored. The Function takes
 *  only one Input arguments, the Alignment Object over two DNA-Strings.
 *  Output will be  String< CigarElement<> > with the CIGAR-Format.
 *
 * Syntax:
 *   String< CigarElement<> > cig = getCigar(Align<Dna5String, ArrayGaps> align);
 *
 * Inputs:
 *   Align<Dna5String, ArrayGaps> align     - Alignment object of seqan-Library with 2 Sequences of Dna5String-Type
 *
 * Outputs:
 *  String< CigarElement<> >                - String of CigarElement-objects representing CIGAR-format of Alignment
 *
 * Example:
 *   String< CigarElement<> > cig = getCigar(align)
 *
 * Other header files required: <seqan/seq_io.h>, <seqan/basic.h>
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
static inline String<CigarElement<> > getCigar(Align<Dna5String, ArrayGaps> &align) {
    typedef Align<Dna5String, ArrayGaps> TAlign;
    String<CigarElement<> > cigar;
    // Use references to the rows of align.
    typedef Row<TAlign>::Type TRow;
    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);
    // Initialize the row iterators.
    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator itRow1 = begin(row1);
    TRowIterator itEndRow1 = end(row1);
    TRowIterator itRow2 = begin(row2);
    TRowIterator itEndRow2 = end(row2);
    unsigned int numChar = 0;
    unsigned long numGaps;
    // don't count beginning Gaps!
    if (isGap(itRow1) || isGap(itRow2)) {
        numGaps = countGaps(itRow1);
        numGaps += countGaps(itRow2);
        itRow1 += numGaps;
        itRow2 += numGaps;
    }
    // don't count ending Gaps! -->
    do {
        --itEndRow1;
        --itEndRow2;
    } while (isGap(itEndRow1) || isGap(itEndRow2));
    // count up last Character from while above which was not a gap
    ++itEndRow1;
    ++itEndRow2;
    // get CIGAR String
    while (itRow1 != itEndRow1) {
        // Count insertions.
        if (isGap(itRow1)) {
            numGaps = countGaps(itRow1);
            CigarElement<> c('I', (unsigned int) numGaps);
            cigar += c;
            itRow1 += numGaps;
            itRow2 += numGaps;
            continue;
        }
        // Count deletions.
        if (isGap(itRow2)) {
            numGaps = countGaps(itRow2);
            CigarElement<> c('D', (unsigned int) numGaps);
            cigar += c;
            itRow1 += numGaps;
            itRow2 += numGaps;
            continue;
        }
        // Count matches.
        while (*itRow1 == *itRow2 && itRow1 != itEndRow1) {
            ++numChar;
            ++itRow1;
            ++itRow2;
        }
        if (numChar != 0) {
            CigarElement<> c('M', numChar);
            cigar += c;
            numChar = 0;
            continue;
        }
        // Count mismatches.
        while (*itRow1 != *itRow2 && itRow1 != itEndRow1) {
            ++numChar;
            ++itRow1;
            ++itRow2;
        }
        if (numChar != 0) {
            CigarElement<> c('S', numChar);
            cigar += c;
            numChar = 0;
        }
    }
    return cigar;
};


/* Function to find a seed-DNAString in a ReferenceDNAString over a Index-Structure given the seed and the index.
 * Function returns all matches in ReferenceString in a vector containing two subvectors. First is BeginPosition of seed.
 * Second is EndPosition of seed. If no match was found Output will be an empty vector.
 *
 * Syntax:
 *   std::vector<long> findAndExtendSeed(TText &seed, Index<TText, TIndex> &index,
 *                                                         TText &read, unsigned int &endPos,
 *                                                         ReadMapperOptions &rmOptions,
 *                                                         std::vector<unsigned long> bestSeed, TDistance)
 *
 * Inputs:
 *   TText seed                     - TText Object of a seed, normally shorter than Reference, i.e. DNA5String
 *   Index<TText, TIndex> &index    - Index-object with Template Parameters to search seed in
 *   TText &read                    - Read over which seed is defined of same Type of seed
 *   unsigned int &endPos           - Ending position of seed in read
 *   ReadMapperOptions &rmOptions   - ReadMapperOptions object containing program parameters
 *   std::vector<unsigned long> bestSeed    - Vector to store position in ReferenceGenome and best Score in
 *   TDistance                      - Template object to set seedExtension Option to runtime
 *
 * Outputs:
 *  vector<unsigned long>  - Vector containing:
 *                  First is beginning Position of best hit of seed in ReferenceDNAString(i.e. Index)
 *                  Second is ending Position of best hit of seed in ReferenceDNAString(i.e. Index)
 *                  Third is score of best hit of seed
 *
 * Example:
 *   std::vector<unsigned long> score = findAndExtendSeed(seed, index, read, 10, rmOptions, bestSeed, Hamming());
 *   std::vector<unsigned long> score = findAndExtendSeed(seed, index, read, 10, rmOptions, bestSeed, Levenshtein());
 *
 * Other header files required: <seqan/index.h>, <seqan/seeds.h>, <seqan/find.h>
 *
 * Subfunctions: extSeed
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
template<typename TText, typename TIndex, typename TDistance>
static inline std::vector<long> findAndExtendSeed(TText &seed, Index<TText, TIndex> &index,
                                                  TText &read, unsigned int &endPos,
                                                  ReadMapperOptions &rmOptions,
                                                  std::vector<long> &bestSeed, TDistance) {
    // Pattern
     StringSet<Dna5String> pat;
    appendValue(pat, seed);

    typedef typename Iterator<StringSet<Dna5String> const, Rooted>::Type TPatternsIt;
    typedef typename Iterator<Index<TText, TIndex>, TopDown<> >::Type TIndexIt;
    // lambda-function. Is called each time find-function finds pattern in text
    auto delegate = [&seed, &index, &read, &endPos, &bestSeed, &rmOptions](TIndexIt const &it,
                                                                           TPatternsIt const &patternsIt,
                                                                           unsigned /*score*/) {
        auto pattern_len = length(*patternsIt);
        // loop through occurences of pattern in text
        for (auto &occ: getOccurrences(it)) {
            // calculate end positions
            long begin = getSeqOffset(occ);
            long end = begin + pattern_len;
            long score = extSeed(seed, index, read, endPos, begin, end, rmOptions, TDistance());
            if (score > bestSeed.back()) {
                bestSeed.at(0) = (long) begin;
                bestSeed.at(1) = (long) end;
                bestSeed.at(2) = score;
            }
        }
    };
    find(index, pat, 1, delegate, Backtracking<EditDistance>());
    return bestSeed;
};


/* Function extend a given seed in both directions. Function returns a vector containing the number of extended matched
 * characters between seed (i.e. String where seed can be found in) and reference-Genome.
 *
 * Syntax:
 * unsigned long extSeed(TText &seed, Index<TText, TIndex> &index, TText &read,
 *                                                unsigned int &endPos,
 *                                                unsigned long startPosition, unsigned long endPosition,
 *                                                ReadMapperOptions &rmOptions, TDistance)
 *
 * Inputs:
 *   TText seed                     - TText Object of a seed, normally shorter than Reference, i.e. DNA5String
 *   Index<TText, TIndex> &index    - Index-object with Template Parameters to search seed in
 *   TText &read                    - Read over which seed is defined of same Type of seed
 *   unsigned int &endPos           - Ending position of seed in read
 *   unsigned long startPosition    - Start Position of seed in ReferenceGenome
 *   unsigned long endPosition      - End Position of seed in ReferenceGenome
 *   ReadMapperOptions &rmOptions   - ReadMapperOptions object containing program parameters
 *   TDistance                      - Template object to set seedExtension Option to runtime
 *
 * Outputs:
 *  unsigned long - Score of Extension
 *
 * Example:
 *   Dna5String seed = "ACGT";
 *   Dna5String index = "CCCCCCCAACGTTGGGGGGG"
 *              // normally this is a index-structure, but to show whats going on this is perfect
 *   Dna5String read = "AAAAAAAAACGTTTTTTT";
 *   unsigned int endPos = 12;
 *   unsigned long startPosition = 9;
 *   unsigned long endPosition = 13;
 *   unsigned long score = extSeed(seed, index, read, endPos, startPosition, endPosition, rmOptions, Hamming());
 *   // score would be 6
 *
 * Other header files required: <seqan/index.h>, <seqan/seeds.h>, <seqan/find.h>
 *
 * Subfunctions: none
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
template<typename TText, typename TIndex, typename TDistance>
inline static long extSeed(TText &seed, Index<TText, TIndex> &index, TText &read,
                           unsigned int &endPos,
                           unsigned long startPosition, unsigned long endPosition,
                           ReadMapperOptions &rmOptions, TDistance) {
    // necessary counters and numbers
    Dna5String refGen;
    unsigned long seedLength = length(seed);
    // define seed
    Seed<Simple> seed1(seedLength, endPos - seedLength, 2 * seedLength, endPos);
    // referenceGenome in Match Region with seedLength characters additional to each side if possible
    unsigned long beginRefGen = (startPosition >= seedLength) ? startPosition - seedLength : 0;
    unsigned long endRefGen = (endPosition + seedLength < length(index)) ? endPosition + seedLength : length(index);
    // get underlying Text of Index at position specified above
    refGen = infix(indexText(index), beginRefGen, endRefGen);
    // call of wrapper Function
    extendWrapper(seed1, refGen, read, rmOptions, TDistance());
    // return score
    return (long) (endPositionH(seed1) - beginPositionH(seed1));
};


/* Function to process a set of given reads. This Function cuts each read in variable number of seeds (seed-length
 * specified in rmOptions) and searches best position of each seed in reference genome. It returns a Vector containing
 * the BamAligmentRecord of each Read with best hit to ReferenceGenome.
 *
 * Syntax:
 *   std::vector<BamAlignmentRecord> processReads(StringSet<Dna5String> &reads, ReadMapperOptions &rmOptions,
 *                                                          Index<TText, TIndex > &index, TDistance)
 *
 * Inputs:
 *   Dna5String &reads              - Dna5String-Type of a seed which can be found in read and ReferenceGenome
 *   ReadMapperOptions &rmOptions   - Options for this Subroutine
 *   Index<TText, TIndex > &index   - Index-object as ReferenceGenome
 *   TDistance                      - Template object to set seedExtension Option to runtime
 *
 * Outputs:
 *   std::vector<BamAlignmentRecord> - Vector containing all BamAlignmentRecords of all Reads
 *
 * Example:
 *   bamRecords = processReads(reads, rmOptions, saIndex, Hamming())
 *   bamRecords = processReads(reads, rmOptions, saIndex, Levenshtein())
 *
 * Other header files required: <seqan/index.h>, <seqan/seeds.h>, <seqan/find.h>
 *
 * Subfunctions: infix, findSeed, extSeed, getCigar
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */

template<typename TText, typename TIndex, typename TDistance>
inline static std::vector<std::vector<BamAlignmentRecord> > processReads(StringSet<Dna5String> &reads,
                                                                         ReadMapperOptions &rmOptions,
                                                                         Index<TText, TIndex> &index, TDistance) {


    // ------------ define objects, counters, etc.
    Dna5String seed;
    unsigned pos = 0;
    unsigned long beginRefGen;
    unsigned long endRefGen;
    unsigned long lengthRead;
    int score;
    String<CigarElement<> > cig;
    std::vector<std::vector<unsigned long> > positions;
    std::vector<unsigned long> matched;
    std::vector<std::vector<BamAlignmentRecord> > bamRecords;



    // ------------ loop over Reads
    // parallelisation
    for (unsigned i = 1; i <= rmOptions.threads; ++i) {
        std::vector<BamAlignmentRecord> bamRecord;
        bamRecords.push_back(bamRecord);
    }
    omp_set_num_threads(rmOptions.threads);

    typedef Iterator<StringSet<Dna5String> >::Type TIterator;

#pragma omp parallel for shared(bamRecords) private(seed) private(pos) private(lengthRead) private(beginRefGen) private(endRefGen) private(score) private(cig) private(positions) private(matched)
    for (unsigned indRead = 0; indRead < length(reads); ++indRead) {
#ifdef verbose
        std::cout << "processing Element: " << position(indRead, reads) << ".." << std::endl;
#endif
        // defining bestSeed Vector for each read
        std::vector<long> bestSeed(3, -1);
        pos = 0;



        // ------------ cut in small seeds
        // Loop through seeds
        typedef Iterator<Dna5String>::Type TIterator;
        for (TIterator seedIt = begin(reads[indRead]); seedIt != end(reads[indRead]);) {
            // if rest of read is bigger than 10 letters
            if ((seedIt + rmOptions.seedLen) <= end(reads[indRead])) {
                seed = infix(*seedIt, seedIt, seedIt + rmOptions.seedLen);
                seedIt += rmOptions.seedLen;
                pos += rmOptions.seedLen;
#ifdef verbose
                std::cout << "Searching for Seed: " << seed << std::endl;
#endif
                // Discard this part ?... could just be one letter... not good to search in Genome!
            } else { // if rest is not bigger than 10 letters
                break;
            }



            // ------------ find & extend & get best seed
            bestSeed = findAndExtendSeed(seed, index, reads[indRead], pos, rmOptions, bestSeed, TDistance());
        } // next seed of read
        // catch nothing found
        if (bestSeed.at(0) != -1) {



            // ------------ perform semi-global Alignment for each Read at position of best seed
            lengthRead = length(reads[indRead]);
            // referenceGenome in Match Region
            beginRefGen = (bestSeed.at(0) >= lengthRead) ? bestSeed.at(0) - lengthRead : 0;
            endRefGen = (bestSeed.at(1) + lengthRead < length(index)) ? bestSeed.at(1) + lengthRead : length(index);
            typedef Align<TText, ArrayGaps> TAlign;
            TAlign align;
            resize(rows(align), 2);
            assignSource(row(align, 0), infix(indexText(index), beginRefGen, endRefGen));
            assignSource(row(align, 1), reads[indRead]);
            // do alignment & get score
            score = globalAlignment(align, rmOptions.alignScheme);
            //score = globalAlignment(align, MyersHirschberg());
#ifdef verbose
            std::cout << align << std::endl;
            std::cout << "Score: " << score << std::endl;
            std::cout << "" << std::endl;
#endif


            // ------------ determine CIGAR Format
            cig = getCigar(align);



            // ------------ Create AlignmentRecord for BAM-File & write in File
            BamAlignmentRecord bamAlignmentRecord;
            // set recordFlags for read
            char name[10];
            sprintf(name, "%s%ld", "Read_", (long) indRead);
            bamAlignmentRecord.qName = name;                                 // Record Name
            bamAlignmentRecord.flag = (uint32_t) indRead;   // Record Number
            bamAlignmentRecord.beginPos = (int32_t) bestSeed.at(0);          // Position in RefGen
            bamAlignmentRecord.cigar = cig;                                  // Alignment CIGAR
            bamAlignmentRecord.seq = reads[indRead];                               // Read sequence
            bamRecords.at(omp_get_thread_num()).push_back(bamAlignmentRecord);
        }
    }// next read

    return bamRecords;
}


/* Main Function of ReadMapper Project
 *
 * Syntax:
 *
 * Inputs:
 *
 * Outputs:
 *
 * Example:
 *
 * Other header files required: <iostream>, <vector>, <seqan/seq_io.h>, <seqan/index.h>, <seqan/seeds.h>
 * Subfunctions: getCigar, buildIndex
 *
 * See also: getCigar, buildIndex
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
int main(int argc, char const **argv) {
    // "/home/phil/Dokumente/ALBIPraktikum/random10M.inx" "/home/phil/Dokumente/ALBIPraktikum/random10M_reads100_10k.fasta" "/home/phil/Dokumente/ALBIPraktikum/testBAM" -i
    //------------ time for statistics
    std::clock_t start;
    double duration;
    start = std::clock();
#define verbose = 1;


    //------------ Parsing Arguments from Command Line
    ReadMapperOptions rmOptions;
    ArgumentParser::ParseResult res = parseCommandLine(rmOptions, argc, argv);
    if (res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR;
    }


    //------------ Read Genome & build/save Index
    if (!rmOptions.loadIndex) {
        // read Fasta ReferenceGenome & build Index
        std::string outRef = rmOptions.inRef;
        std::string filename = ".inx";
        outRef.append(filename);
        if (!buildIndex(rmOptions.inRef, outRef, rmOptions)) {
            std::cerr << "ERROR: Error due to Message above. Exiting.." << std::endl;
            return 0;
        }
        rmOptions.inRef = outRef;
    }

    // check for buildIndexOnly parameter
    if (!rmOptions.buildIndexOnly) {
        //------------ Read Index & Reads
        std::cout << "Loading Index.." << std::endl;
        // load Index from disk
        Index<Dna5String, IndexSa<> > saIndex;
        open(saIndex, toCString(rmOptions.inRef));
        std::cout << "Index loaded!" << std::endl;



        //------------ read Reads-Records
        std::cout << "Reading Reads.." << std::endl;
        SeqFileIn rFileIn;
        if (!open(rFileIn, toCString(rmOptions.inReads))) {
            std::cerr << "ERROR: Could not read Reads-File!" << std::endl;
        }
        StringSet<Dna5String> reads;
        StringSet<CharString> heads;
        try {
            // try to read all records
            readRecords(heads, reads, rFileIn);
        } catch (IOError io) {
            std::cerr << "ERROR: " << io.what() << std::endl;
        } catch (ParseError p) {
            std::cerr << "ERROR: " << p.what() << std::endl;
        }
        std::cout << "Reads read!" << std::endl;



        //------------ Processing reads
        std::cout << "Processing.." << std::endl;
        std::vector<std::vector<BamAlignmentRecord> > bamRecords;
        if (!rmOptions.levenshtein) {
            bamRecords = processReads(reads, rmOptions, saIndex, Hamming());
        } else {
            bamRecords = processReads(reads, rmOptions, saIndex, Levenshtein());
        }
        std::cout << "Processed!" << std::endl;



        //------------  define files for writing out & write Header
        std::cout << "Writing BAM.." << std::endl;
        BamFileIn bamFileIn;
        std::ofstream bamFile;
        bamFile.open(toCString(rmOptions.outSam));
        BamFileOut bamFileOut(context(bamFileIn), bamFile, Sam());
        // Header for BAM File
        BamHeader bamHeader;
        BamHeaderRecord bamHeaderRecord;
        setTagValue("VN", "1.3", bamHeaderRecord);
        setTagValue("SO", "unsorted", bamHeaderRecord);
        setTagValue("SN", "ref", bamHeaderRecord);
        char len[10];
        sprintf(len, "%ld", length(saIndex));
        setTagValue("LN", len, bamHeaderRecord);
        assign(bamHeader, bamHeaderRecord);
        // write in File
        writeHeader(bamFileOut, bamHeader);


        // ------------ write alignment-records in File & close
        // Iterat over vectors containing vectors of BamRecords
        std::vector<std::vector<BamAlignmentRecord> >::iterator itBamRecords;
        for (itBamRecords = bamRecords.begin(); itBamRecords != bamRecords.end(); ++itBamRecords) {
            std::vector<BamAlignmentRecord>::iterator itBamRecord;
            for (itBamRecord = (*itBamRecords).begin(); itBamRecord != (*itBamRecords).end(); ++itBamRecord) {
                writeRecord(bamFileOut, *itBamRecord);
            }
        }
        bamFile.close();
        std::cout << "BAM written!" << std::endl;
    }
    std::cout << "Finished!" << std::endl;
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Runtime: " << duration << std::endl;
    sleep(1);
    return EXIT_SUCCESS;
}
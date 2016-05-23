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
 *   std::string inRef      - path to ReferenceFile. Can be a Fasta or Index File. If Index File please set Option -i
 *   std::string inReads    - path to ReadsFile.
 *   std::string outSam     - path to Output File
 *   bool loadIndex         - boolean indicating if inRef is Fasta or Index File
 *   bool levenshtein       - boolean indicating if distance is set to edit distance
 *   bool fmInd             - boolean indicating if index mode is set to FM-Index
 *   bool buildIndexOnly    - boolean indicating if program should only build index file
 *   bool matchExtend       - boolean indicating if dropout threshold is set
 *   int seedLen            - declare seed Length to split Reads
 *   int dropout            - dropout threshold. Seed Extension will stop when score falls under threshold
 *   int threads            - number of threads to use for processing reads
 *   int approxSearch       - alowed mistakes for searching seed in genome
 *   Score<int, Simple> extendScheme - alignment scheme for seed extension in seed matching region
 *   Score<int, Simple> alignScheme  - alignment scheme for Alignment between Read and ReferenceGenome
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
    bool matchExtend;
    int seedLen;
    int dropout;
    int threads;
    int approxSearch;
    Score<int, Simple> alignScheme;
    Score<int, Simple> extendScheme;


    ReadMapperOptions() :
            inRef(""), inReads(""), outSam(""), loadIndex(false), fmInd(false), seedLen(10),
            levenshtein(false), alignScheme(1, -1, 0, -1), dropout(3), approxSearch(0),extendScheme(2, -1, -1, -2), threads(1) { }

};


/* Struct to set Distance Object in runtime */
struct Levenshtein {
    Levenshtein() { }
};


/* Struct to set Distance Object in runtime */
struct Hamming {
    Hamming() { }
};

/* Struct to set seedExtesnionMode Object in runtime */
struct MatchExt {
    MatchExt() { }
};

/* Struct to set seedExtesnionMode Object in runtime */
struct XDrop {
    XDrop() { }
};


/* Function to call the extendSeed-Function with to runtime specialized options
 *
 * Syntax:
 *   extendWrapper(Seed<Simple> & seed1 , Dna5String &refGen,
 *                               ReadMapperOptions &rmOptions, XDrop )
 *
 * Inputs:
 *   Seed<Simple> & seed1           - seed to extend
 *   Dna5String &refGen             - reference Gen to extend seed over
 *   TText &read                    - read to extend Seed over
 *   ReadMapperOptions &rmOptions   - options containing AlignmentScheme and dropOut-threshold
 *   XDrop                          - object telling the compiler wich version of extendSeed to call
 *
 * Outputs:
 *
 * Example:
 *   extendWrapper(seed1, refGen, rmOptions, XDrop())
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
                                 ReadMapperOptions &rmOptions, XDrop /**/) {
    extendSeed(seed1, refGen, read, EXTEND_BOTH, rmOptions.extendScheme, 3, GappedXDrop());
}


/* Function to call the extendSeed-Function with to runtime specialized options
 *
 * Syntax:
 *   extendWrapper(Seed<Simple> & seed1 , Dna5String &refGen,
 *                               ReadMapperOptions &rmOptions, MatchExt )
 *
 * Inputs:
 *   Seed<Simple> & seed1           - seed to extend
 *   Dna5String &refGen             - reference Gen to extend seed over
 *   TText &read                    - read to extend Seed over
 *   ReadMapperOptions &rmOptions   - options containing AlignmentScheme and dropOut-threshold
 *   MatchExt                       - object telling the compiler which version of extendSeed to call
 *
 * Outputs:
 *
 * Example:
 *   extendWrapper(seed1, refGen, rmOptions, MatchExt())
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
                                 ReadMapperOptions &rmOptions, MatchExt /**/) {
    extendSeed(seed1, refGen, read, EXTEND_BOTH, MatchExtend());
}


/* Function to call the searchSeed-Function with to runtime specialized options
 *
 * Syntax:
 *   searchWrapper(TText &seed, StringSet<TText> &pat, Index<TText, TIndex> &index,
 *                                            TText &read, unsigned int &endPos,
 *                                            ReadMapperOptions &rmOptions,
 *                                            std::vector<long> &bestSeed, TExtendMode, Levenshtein) {
 *
 * Inputs:
 *   TText & seed                   - seed to extend
 *   StringSet<TText> &pat          - pattern to search for
 *   Index<TText, TIndex> &index    - index to search over
 *   TText &read                    - read in wich seed is located
 *   unsigned int &endPos           - endPosition of seed in read
 *   ReadMapperOptions &rmOptions   - options containing AlignmentScheme and dropOut-threshold etc.
 *   std::vector<long> &bestSeed    - vector containing information about best seed found so far
 *   TExtendMode                    - object telling the compiler which version of extendSeed to call
 *   Levenshtein                    - object telling the compiler which version of searchWrapper to call
 *
 * Outputs:
 *  std::vector<long>               - Vector containing:
 *                  First is beginning Position of best hit of seed in ReferenceDNAString(i.e. Index)
 *                  Second is ending Position of best hit of seed in ReferenceDNAString(i.e. Index)
 *                  Third is score of best hit of seed
 *
 * Example:
 *   searchWrapper(seed, pat, index, read, 10, rmOptions, bestSeed, XDrop(), Levenshtein())
 *
 * Other header files required: none
 * Subfunctions: extSeed
 *
 * See also:
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
template<typename TText, typename TIndex, typename TExtendMode>
inline static std::vector<long> searchWrapper(TText &seed, StringSet<TText> &pat, Index<TText, TIndex> &index,
                                              TText &read, unsigned int &endPos,
                                              ReadMapperOptions &rmOptions,
                                              std::vector<long> &bestSeed, TExtendMode /**/, Levenshtein /**/) {

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
            long score = extSeed(seed, index, read, endPos, begin, end, rmOptions, TExtendMode());
            if (score > bestSeed.back()) {
                bestSeed.at(0) = begin;
                bestSeed.at(1) = end;
                bestSeed.at(2) = score;
            }
        }
    };
    find(index, pat, 1, delegate, Backtracking<EditDistance>());
    return bestSeed;
}

/* Function to call the searchSeed-Function with to runtime specialized options
*
* Syntax:
*   searchWrapper(TText &seed, StringSet<TText> &pat, Index<TText, TIndex> &index,
*                                            TText &read, unsigned int &endPos,
*                                            ReadMapperOptions &rmOptions,
*                                            std::vector<long> &bestSeed, TExtendMode, Hamming) {
*
* Inputs:
*   TText & seed                   - seed to extend
*   StringSet<TText> &pat          - pattern to search for
*   Index<TText, TIndex> &index    - index to search over
*   TText &read                    - read in wich seed is located
*   unsigned int &endPos           - endPosition of seed in read
*   ReadMapperOptions &rmOptions   - options containing AlignmentScheme and dropOut-threshold etc.
*   std::vector<long> &bestSeed    - vector containing information about best seed found so far
*   TExtendMode                    - object telling the compiler which version of extendSeed to call
*   Hamming                        - object telling the compiler which version of searchWrapper to call
*
* Outputs:
*  std::vector<long>               - Vector containing:
*                  First is beginning Position of best hit of seed in ReferenceDNAString(i.e. Index)
*                  Second is ending Position of best hit of seed in ReferenceDNAString(i.e. Index)
*                  Third is score of best hit of seed
*
* Example:
*   searchWrapper(seed, pat, index, read, 10, rmOptions, bestSeed, XDrop(), Hamming())
*
* Other header files required: none
* Subfunctions: extSeed
*
* See also:
* Author: Jan Philipp Albrecht
* Work address:
* email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
* Website: https://github.com/jpalbrecht/SWPraktikum
*/
template<typename TText, typename TIndex, typename TExtendMode>
inline static std::vector<long> searchWrapper(TText &seed, StringSet<TText> &pat, Index<TText, TIndex> &index,
                                              TText &read, unsigned int &endPos,
                                              ReadMapperOptions &rmOptions,
                                              std::vector<long> &bestSeed, TExtendMode /**/, Hamming /**/) {
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
            long score = extSeed(seed, index, read, endPos, begin, end, rmOptions, TExtendMode());
            if (score > bestSeed.back()) {
                bestSeed.at(0) = begin;
                bestSeed.at(1) = end;
                bestSeed.at(2) = score;
            }
        }
    };
    find(index, pat, rmOptions.approxSearch, delegate, Backtracking<HammingDistance>());;
    return bestSeed;
}


/* Function to parse Arguments from Command Line
 *
 * Syntax:
 *   ArgumentParser::ParseResult parseCommandLine(ReadMapperOptions &rmOptions, int argc, char const **argv)
 *
 * Inputs:
 *   ReadMapperOptions &rmOptions   - optionObject from ReadMapperOption Class.
 *                                      See Struct Documentation for more information
 *   int argc                       - number of Arguments
 *   char const **argv              - arguments from Command Line
 *
 * Outputs:
 *  ArgumentParser::ParseResult     - result of parsing Arguments. Can be PARSE_OK or PARSE_ERROR.
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
    addOption(parser, ArgParseOption("i", "index", "Option to load ReferenceGenome out of IndexFile instead out of FASTA. Default: FASTA"));
    addOption(parser, ArgParseOption("b", "buildIndex", "Option to just build the Index File"));
    addOption(parser, ArgParseOption("r", "readFile", "set Path to File containing Reads. FASTA or FASTQ",
                                     ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption("o", "outSam", "set Path to output File. Default: readFile path",
                                     ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption("t", "threads", "set Number of threads to use. Default: 1",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("s", "seedLength", "Number of Characters of seed. Default: 10",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("S", "suffixArray", "set the Index Modus to SuffixArray. Default: SuffixArray"));
    addOption(parser, ArgParseOption("F", "fmIndex", "set the Index Modus to FM-Index. Default: SuffixArray"));
    addOption(parser, ArgParseOption("l", "levenshtein", "set the extension Distance to Levenshtein-Distance. Default: Hamming distance"));
    addOption(parser, ArgParseOption("d", "dropout", "set the dropout threshold when to stop extend seeds Default: exact search",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("a", "approximateSearch", "set the allowed mistakes to search seed with. Default: 0",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ha", "hamming", "setting the extension Distance to Hamming-Distance. Default: Hamming distance"));
    addOption(parser, ArgParseOption("ma", "match", "MatchCost. Default: 1",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mm", "mismatch", "MisMatchCost. Default: -1 ",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("go", "gapOpen", "GapOpenCost. Default: -1",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ge", "gapExtend", "GapExtendCost. Default: 0",
                                     ArgParseArgument::INTEGER, "INT"));


    //------------ Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutIndex\\fP \\fB-i\\fP ",
                "Load ReferenceGenome out of IndexFile instead out of FASTA-File.");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP  \\fB-s\\fP \\fI10\\fP ",
                "Set seedLength to 10(default)");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-b\\fP \\fB-S\\fP ",
                "Set Index Modus to SuffixArray. Build Suffix Array only-Modus. No Reads are processed.");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-b\\fP \\fB-F\\fP ",
                "Set Index Modus to FM-Index. Build FM-Index only-Modus. No Reads are processed.");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-t\\fP \\fB4\\fP",
                "Set the search mode to Levenshtein-Distance. Use 4 threads for searching");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-l\\fP \\fB-d\\fP \\fI3\\fP",
                "Set the droput threshold to 3. Found seed will be extended till score falls below threshold.");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB../pathToInPutFile\\fP \\fB-r\\fP \\fB../pathToReadsFile\\fP \\fB-ha\\fP ",
                "Set the search mode to Hamming-Distance(default)");
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
    rmOptions.matchExtend = !isSet(parser,"dropout");
    getOptionValue(rmOptions.approxSearch, parser, "approximateSearch");


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


    //------------ check dissent options
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
    if (isSet(parser, "approximateSearch")) {
        if (rmOptions.approxSearch < 0){
            std::cout << "ERROR: Cannot search for seed with less then 0 errors accepted!\n";
            return seqan::ArgumentParser::PARSE_ERROR;
        }
        if (rmOptions.approxSearch > 3) {
            std::cout << "WARNING: Searching for seed with more than 3 errors accepted might take a long time!\n";
        }
    }

    //------------ display warinings and return
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
 *   int erg = buildIndex(CharString fileIn, CharString fileOut, ReadMapperOptions rmOptions);
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
 *   Align<Dna5String, ArrayGaps> align     - alignment object of seqan-Library with 2 Sequences of Dna5String-Type
 *
 * Outputs:
 *  String< CigarElement<> >                - string of CigarElement-objects representing CIGAR-format of Alignment
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


/* Function extends a given seed in both directions. Function returns a score of extension
 *
 * Syntax:
 * unsigned long extSeed(TText &seed, Index<TText, TIndex> &index, TText &read,
 *                                                unsigned int &endPos,
 *                                                unsigned long startPosition, unsigned long endPosition,
 *                                                ReadMapperOptions &rmOptions, TExtendMode)
 *
 * Inputs:
 *   TText seed                     - TText Object of a seed, normally shorter than Reference, i.e. DNA5String
 *   Index<TText, TIndex> &index    - index-object with Template Parameters to search seed in
 *   TText &read                    - read over which seed is defined of same Type of seed
 *   unsigned int &endPos           - ending position of seed in read
 *   unsigned long startPosition    - start Position of seed in ReferenceGenome
 *   unsigned long endPosition      - end Position of seed in ReferenceGenome
 *   ReadMapperOptions &rmOptions   - readMapperOptions object containing program parameters
 *   TExtendMode                    - template object to set seedExtension Option to runtime
 *
 * Outputs:
 *  long - Score of Extension
 *
 * Example:
 *   Dna5String seed = "ACGT";
 *   Dna5String index = "CCCCCCCAACGTTGGGGGGG"
 *              // normally this is a index-structure, but to show whats going on this is perfect
 *   Dna5String read = "AAAAAAAAACGTTTTTTT";
 *   unsigned int endPos = 12;
 *   unsigned long startPosition = 9;
 *   unsigned long endPosition = 13;
 *   unsigned long score = extSeed(seed, index, read, endPos, startPosition, endPosition, rmOptions, XDrop());
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
template<typename TText, typename TIndex, typename TExtendMode>
inline static long extSeed(TText &seed, Index<TText, TIndex> &index, TText &read,
                           unsigned int &endPos,
                           unsigned long startPosition, unsigned long endPosition,
                           ReadMapperOptions &rmOptions, TExtendMode /**/) {
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
    extendWrapper(seed1, refGen, read, rmOptions, TExtendMode());
    // return score
    return (long) (endPositionH(seed1) - beginPositionH(seed1));
};


/* Function to process a set of given reads. This Function cuts each read in variable number of seeds (seed-length
 * specified in rmOptions) and searches best position of each seed in reference genome. It returns a Vector containing
 * the BamAligmentRecord of each Read with best hit to ReferenceGenome.
 *
 * Syntax:
 *   std::vector<BamAlignmentRecord> processReads(StringSet<Dna5String> &reads, ReadMapperOptions &rmOptions,
 *                                                          Index<TText, TIndex > &index, TExtendMode, TDistance)
 *
 * Inputs:
 *   Dna5String &reads              - Dna5String-Type of a seed which can be found in read and ReferenceGenome
 *   ReadMapperOptions &rmOptions   - options for this Subroutine
 *   Index<TText, TIndex > &index   - index-object as ReferenceGenome
 *   TExtendMode                    - template object to set seedExtension Option to runtime
 *   TDistance                      - template object to set distance Option to runtime
 *
 * Outputs:
 *   std::vector< std::vector<BamAlignmentRecord> > - vector containing subvectors containing all
 *                                  BamAlignmentRecords of all Reads. There are as many subvectors as threads are used.
 *
 * Example:
 *   bamRecords = processReads(reads, rmOptions, saIndex, XDrop(), Hamming())
 *   bamRecords = processReads(reads, rmOptions, saIndex, MatchExt(), Levenshtein())
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
template<typename TText, typename TIndex, typename TExtendMode, typename TDistance>
inline static std::vector<std::vector<BamAlignmentRecord> > processReads(StringSet<Dna5String> &reads,
                                                                         ReadMapperOptions &rmOptions,
                                                                         Index<TText, TIndex> &index, TExtendMode /**/,
                                                                         TDistance /**/) {

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
            // Pattern
            StringSet<Dna5String> pat;
            appendValue(pat, seed);

            // search
            bestSeed = searchWrapper(seed, pat, index, reads[indRead], pos, rmOptions, bestSeed, TExtendMode(),
                                     TDistance());

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
            bamAlignmentRecord.seq = reads[indRead];                         // Read sequence
            bamRecords.at((unsigned long)omp_get_thread_num()).push_back(bamAlignmentRecord);
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
 * See also: getCigar, buildIndex, processReads
 * Author: Jan Philipp Albrecht
 * Work address:
 * email: jan-philipp.albrecht@charite.de, j.p.albrecht@fu-berlin.de
 * Website: https://github.com/jpalbrecht/SWPraktikum
 */
int main(int argc, char const **argv) {
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
            if(!rmOptions.matchExtend){
                bamRecords = processReads(reads, rmOptions, saIndex, XDrop(), Hamming());
            }else {
                bamRecords = processReads(reads, rmOptions, saIndex, MatchExt(), Hamming());
            }

        } else {
            if(!rmOptions.matchExtend){
                bamRecords = processReads(reads, rmOptions, saIndex, XDrop(), Levenshtein());
            }else {
                bamRecords = processReads(reads, rmOptions, saIndex, MatchExt(), Levenshtein());
            }
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
    std::cout.flush();
    return EXIT_SUCCESS;
}
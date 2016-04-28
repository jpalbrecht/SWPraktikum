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
 *                                  of Levenshtein Distance. False to Edit-Distance
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
    int seedLen;
    int indexMode;
    bool levenshtein;
    Score<int, Simple> alignScheme;

    ReadMapperOptions() :
            inRef(""), inReads(""), outSam(""), loadIndex(false),seedLen(10),indexMode(0),
            levenshtein(false), alignScheme(1, -1, 0, -1)  { }

};


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
    setDate(parser, "April 2016");
    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIInputReference\\fP] \"\\fIInputReads\\fP\"  \"\\fIOutputSam\\fP\" \"\\fIOPTIONS\\fP\"  ]");
    addDescription(parser, "This program allows to map Reads to a Reference Genome.");
    // 3 arguments required
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "InPathRef"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "InPathReads"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "OutPathSam"));


    //------------ Define Options -- Section Modification Options
    addSection(parser, "Modification Options");
    addOption(parser, ArgParseOption("i", "index", "Option to load ReferenceGenome out of IndexFile."));
    addOption(parser, ArgParseOption("s", "seedLength", "Number of Characters of seed",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("m", "indexMode", "Indicating the Index Modus to use. 0:SuffixArray 1:FM-Index",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("l", "levenshtein", "setting the extension Distance to Levenshtein-Distance"));
    addOption(parser, ArgParseOption("e", "edit", "setting the extension Distance to Edit-Distance"));
    addOption(parser, ArgParseOption("ma", "match", "MatchCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("mm", "missmatch", "MissMatchCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("go", "gapOpen", "GapOpenCost",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ge", "gapExtend", "GapExtendCost",
                                     ArgParseArgument::INTEGER, "INT"));


    //------------ Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-i\\fP ",
                "Load ReferenceGenome out of IndexFile");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-s\\fP \\fI10\\fP ",
                "Set seedLength to 10(default)");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-m\\fP \\fI0\\fP ",
                "Set Index Modus to SuffixArray");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-l\\fP",
                "Set the extendSeed Mode to Levenshtein-Distance");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-e\\fP",
                "Set the extendSeed Mode to Edit-Distance");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-ma\\fP \\fI2\\fP ",
                "Set MatchCost to 2");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-mm\\fP \\fI-1\\fP ",
                "Set MissMatchCost to -1");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-go\\fP \\fI-2\\fP ",
                "Set GapOpenCost to -2");
    addListItem(parser,
                "\\fBreadMapper\\fP \\fB-ge\\fP \\fI-1\\fP ",
                "Set GapExtendCost to -1");


    //------------ Parse command line & read necessary Arguments and optional Values
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;
    // read Values
    getArgumentValue(rmOptions.inRef, parser, 0);
    getArgumentValue(rmOptions.inReads, parser, 1);
    getArgumentValue(rmOptions.outSam, parser, 2);
    rmOptions.loadIndex = isSet(parser, "index");
    getOptionValue(rmOptions.seedLen, parser, "seedLength");
    getOptionValue(rmOptions.indexMode, parser, "indexMode");
    rmOptions.levenshtein = isSet(parser, "levenshtein");


    //------------ read ScoringScheme if set
    int ma, mm, go, ge;
    getOptionValue(ma, parser, "match");
    getOptionValue(mm, parser, "missmatch");
    getOptionValue(go, parser, "gapOpen");
    getOptionValue(ge, parser, "gapExtend");
    if (isSet(parser, "match")){
        setScoreMatch(rmOptions.alignScheme, ma);
    }
    if (isSet(parser, "missmatch")){
        setScoreMismatch(rmOptions.alignScheme, mm);
    }
    if (isSet(parser, "gapOpen")){
        setScoreGapOpen(rmOptions.alignScheme, go);
    }
    if (isSet(parser, "gapExtend")){
        setScoreGapExtend(rmOptions.alignScheme, ge);
    }


    //------------ check dissent options & return
    if (isSet(parser, "levenshtein") && isSet(parser, "edit"))
    {
        std::cerr << "ERROR: You cannot specify both levenshtein- and edit-Distance!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
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
 *   CharString fileIn      - complete Path to fastA file
 *   CharString fileOut     - complete Path to output index file
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
static inline int buildIndex(std::string &fileIn, std::string &fileOut,ReadMapperOptions &rmOptions) {
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
    if (rmOptions.indexMode){
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
 *   Align<Dna5String, ArrayGaps> align  - Alignment object of seqan-Library with 2 Sequences of Dna5String-Type
 *
 * Outputs:
 *  String< CigarElement<> >    - String of CigarElement-objects representing CIGAR-format of Alignment
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
 *   std::vector<std::vector<unsigned long> > findSeed(Dna5String seed, Index < Dna5String, IndexSa<> > index );
 *
 * Inputs:
 *   Dna5String seed  - Dna5String-Type of a seed, normally shorter than Reference
 *   Index < Dna5String, IndexSa<> > index  -  Index-object to search seed in
 *
 * Outputs:
 *  std::vector<std::vector<unsigned long> >  - Vector containing two other vectors.
 *                  First is/are beginning Position(s) of hits of seed in ReferenceDNAString(i.e. Index)
 *                  Second is/are ending Position(s) of hits of seed in ReferenceDNAString(i.e. Index)
 *
 * Example:
 *   std::vector<std::vector<unsigned long> > matches = findSeed(seed, index);
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
static inline std::vector<std::vector<unsigned long> > findSeed(Dna5String &seed, Index<Dna5String, IndexSa<> > &index) {
    // return object
    std::vector<std::vector<unsigned long> > position;
    // Pattern & Finder
    Pattern<Dna5String> pat(seed);
    Finder<Index<Dna5String, IndexSa<> > > finder(index);
    // initialize Vector to store seedPosition in
    std::vector<unsigned long> begPos;
    std::vector<unsigned long> begEnd;
    while (find(finder, pat)) {
        begPos.push_back(beginPosition(finder));
        begEnd.push_back(endPosition(finder));
//                std::cout << '[' << beginPosition(finder) << '.' << endPosition(finder) << ']' << std::endl;
    }
    if (begPos.empty()) {
//                std::cout << "nothing Found!.. next seed" << std::endl;
        // retuning empty Vector
        return position;
    }
//            std::cout << "Found Seed " << begPos.size() <<" times!" << std::endl;

    position.push_back(begPos);
    position.push_back(begEnd);
    return position;
};


/* Function extend a given seed in both directions. Function returns a vector containing the number of extended matched
 * characters between seed (i.e. String where seed can be found in) and reference-Genome.
 *
 * Syntax:
 * static std::vector<unsigned long> extendSeed(Dna5String &seed, Index<Dna5String, IndexSa<> > &index, Dna5String &read,
 *                                             unsigned int &endPos,
 *                                             std::vector<std::vector<unsigned long> > &positions)
 *
 * Inputs:
 *   Dna5String seed  - Dna5String-Type of a seed which can be found in read and ReferenceGenome
 *   Index < Dna5String, IndexSa<> > index  -  Index-object as ReferenceGenome
 *   Dna5String read  -  Dna5String-Type of a read where the seed can be found in
 *   unsigned int endPos - Position in the Read where seed ends! CAUTION: ENDPOSITION!
 *   std::vector<std::vector<unsigned long> > vector containing two subvectors.
 *                  First is/are beginning Position(s) of hits of seed in ReferenceDNAString(i.e. Index)
 *                  Secound is/are ending Position(s) of hits of seed in ReferenceDNAString(i.e. Index)
 *
 * Outputs:
 *  std::vector<unsigned long> - Vector containing the total number of characters which can be aligned to
 *                  ReferenceGenome and Read by extending the Seed-Region.
 *
 * Example:
 *   Dna5String seed = "ACGT";
 *   Dna5String index = "CCCCCCCAACGTTGGGGGGG" // normally this is a index-structure, but to show whats going on this is perfect
 *   Dna5String read = "AAAAAAAAACGTTTTTTT";
 *   std::vector<unsigned long> begin; begin.push_back = 9;
 *   std::vector<unsigned long> end; end.push_back = 13;
 *   std::vector<std::vector<unsigned long> > positions; positions.push_back(begin); positions.push_back(end);
 *   std::vector<unsigned long> matches = findSeed(seed, index, read, 12, positions); // would contain number 6
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
template <typename TText, typename TIndex>
inline static std::vector<unsigned long> extendSeed(TText &seed, Index<TText, TIndex > &index, TText &read,
                                             unsigned int &endPos,
                                             std::vector<std::vector<unsigned long> > &positions,
                                                    ReadMapperOptions &rmOptions) {
    // necessary counters and numbers
    Dna5String refGen;
    unsigned long beginRefGen;
    unsigned long endRefGen;
    unsigned long seedLength = length(seed);
    // vector to store extend-length in
    std::vector<unsigned long> matched;
    // loop throught all hits of a given seed at different positions in ReferenceGenome
    for (unsigned ind = 0; ind < positions.at(0).size(); ind++) {
        // define seed
        Seed<Simple> seed1(seedLength, endPos - seedLength, 2 * seedLength, endPos);
        // referenceGenome in Match Region with seedLength characters additional to each side if possible
        beginRefGen = (positions.at(0).at(ind) >= seedLength) ? positions.at(0).at(ind) - seedLength : 0;
        endRefGen = (positions.at(1).at(ind) + seedLength < length(index)) ? positions.at(1).at(ind) + seedLength
                                                                           : length(index);
        // get underlying Text of Index at position specified above
        refGen = infix(indexText(index), beginRefGen, endRefGen);
        if (rmOptions.levenshtein) {
            // match extension LEVENSHTEIN!
            Score<int, Simple> scoringScheme(2, -1, -1, -2);
            extendSeed(seed1, refGen, read, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
        } else {
            // match extension EDIT-Distance!
            extendSeed(seed1, refGen, read, EXTEND_BOTH, MatchExtend());
        }
        matched.push_back(endPositionH(seed1) - beginPositionH(seed1));
//                std::cout << "matched: " << matched.back() << std::endl;
    }// next seed to extend
    return matched;
};

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
    // mehr Index direkt in cache wenn ich größe verkleinere
    //------------ time for statistics
    std::clock_t start;
    double duration;
    start = std::clock();



    //------------ Parsing Arguments from Command Line
    ReadMapperOptions rmOptions;
    ArgumentParser::ParseResult res = parseCommandLine(rmOptions, argc, argv);
    if (res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR;
    }



    //------------ Read Genome & build/save Index
    if (!rmOptions.loadIndex){
        // read Fasta ReferenceGenome & build Index
        std::string outRef = rmOptions.inRef;
        std::string filename = ".inx";
        outRef.append(filename);
        if (!buildIndex(rmOptions.inRef,outRef, rmOptions )) {
            std::cerr << "ERROR: Error due to Message above. Exiting.." << std::endl;
            return 0;
        }
        rmOptions.inRef = outRef;
    }



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



    //------------  define files for writing out & write Header
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



    //------------ Loop through reads
    std::cout << "Processing.." << std::endl;
    typedef Iterator<StringSet<Dna5String> >::Type TIterator;
    for (TIterator readsIt = begin(reads); readsIt != end(reads); ++readsIt) {
//        std::cout << "processing Element: " << position(readsIt, reads) << ".." << std::endl;



        // ------------ cut in small seeds
        Dna5String seed;
        unsigned pos = 0;
        // building best seed Vector & counting-Variables
        std::vector<unsigned long> bestSeed(3, 0);
        unsigned long beginRefGen;
        unsigned long endRefGen;
        // Loop through seeds
        typedef Iterator<Dna5String>::Type TIterator;
        for (TIterator seedIt = begin(*readsIt); seedIt != end(*readsIt);) {
            // if rest of read is bigger than 10 letters
            if ((seedIt + rmOptions.seedLen) <= end(*readsIt)) {
                seed = infix(*seedIt, seedIt, seedIt + rmOptions.seedLen);
                seedIt += rmOptions.seedLen;
                pos += rmOptions.seedLen;
//                std::cout << "Searching for Seed: " << seed << std::endl;
                // Discard this part ?... could just be one letter... not good to search in Genome!
            } else { // if rest is not bigger than 10 letters
                // trotzdem suchen und schauen wieviel hits ich habe und nach threshold aussortieren.
//              seed = infix(*readIt, readIt, end(*readsIt));
                seedIt = end(*readsIt);
//              std::cout << "Searching for Seed: " << seed << std::endl;
                break;
            }



            // ------------ find seed
            std::vector<std::vector<unsigned long> > positions = findSeed(seed, saIndex);
            // catch nothing found case
            if (positions.empty()) {
                // next seed
                break;
            }



            // ------------ extend seed(s)
            std::vector<unsigned long> matched = extendSeed(seed, saIndex, *readsIt, pos, positions, rmOptions);



            // ------------ save best seed
            std::vector<unsigned long>::iterator bestMatch = std::max_element(matched.begin(), matched.end());
            unsigned long posBestMatch = (unsigned long) std::distance(matched.begin(), bestMatch);
//          std::cout << "best match with " << *bestMatch << " Bases!" << std::endl;
            // save best seed if better than these before
            if (*bestMatch > bestSeed.back()) {
                bestSeed.at(0) = positions.at(0).at(posBestMatch);
                bestSeed.at(1) = positions.at(1).at(posBestMatch);
                bestSeed.at(2) = *bestMatch;
            }

        } // next seed of read



        // ------------ perform semi-global Alignment for each Read at position of best seed
        unsigned long lengthRead = length(*readsIt);
        // referenceGenome in Match Region
        beginRefGen = (bestSeed.at(0) >= lengthRead) ? bestSeed.at(0) - lengthRead : 0;
        endRefGen = (bestSeed.at(1) + lengthRead < length(saIndex)) ? bestSeed.at(1) + lengthRead : length(saIndex);
        typedef Align<Dna5String, ArrayGaps> TAlign;
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(indexText(saIndex), beginRefGen, endRefGen));
        assignSource(row(align, 1), *readsIt);
        // do alignment & get score
        int score = globalAlignment(align, rmOptions.alignScheme);
//        std::cout << align << std::endl;
//        std::cout << "Score: " << score << std::endl;
//        std::cout << "" << std::endl;



        // ------------ determine CIGAR Format
        String<CigarElement<> > cig = getCigar(align);



        // ------------ Create AlignmentRecord for BAM-File & write in File
        BamAlignmentRecord bamAlignmentRecord;
        // set recordFlags for read
        char name[10];
        sprintf(name, "%s%ld", "Read_", position(readsIt, reads));
        bamAlignmentRecord.qName = name;                                // Record Name
        bamAlignmentRecord.flag = (uint32_t) position(readsIt, reads);   // Record Number
        bamAlignmentRecord.beginPos = (int32_t) bestSeed.at(0);          // Position in RefGen
        bamAlignmentRecord.cigar = cig;                                 // Alignment CIGAR
        bamAlignmentRecord.seq = *readsIt;                             // Read sequence
        //write alignment-record in File
        writeRecord(bamFileOut, bamAlignmentRecord);


    }// next read



    bamFile.close();
    std::cout << "Finished!" << std::endl;
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Runtime: " << duration << std::endl;
    return 1;
}
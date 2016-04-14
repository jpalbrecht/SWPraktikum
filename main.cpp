#include <iostream>
#include <vector>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/bam_io.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/stream.h>


using namespace seqan;

// Function to read a FastaFile in & building/saving Index-Structure of content of file.
// Returns true if file could be read and index could be built and could be wrote to disk, otherwise false
static int buildIndex(CharString fileIn, CharString fileOut) {
    //------------ Read in FastaReferenceFile & print Info
    SeqFileIn sFileIn;
    std::cout << "Reading File.." << std::endl;
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
    Index < Dna5String, IndexSa<> > saIndex(genome);
    // require Index for being able to save to disk
    indexRequire(saIndex, FibreSA());
    std::cout << "Index built!" << std::endl;
    // save Index to File!
    std::cout << "Saving Index to disk.." << std::endl;
    if (!save(saIndex, toCString(fileOut))) {
        std::cerr << "ERROR: not able to write Index-File to disk. Check for available space & authorization" <<
        std::endl;
        return 0;
    }
    std::cout << "Index saved!" << std::endl;

    return 1;
};


int main() {
    //------------ Read Index & Reads read Genome & build/save Index
//    if (!buildIndex("/home/phil/Dokumente/ALBIPraktikum/random10M.fasta",
//    "/home/phil/Dokumente/ALBIPraktikum/random10M.inx")){
//        std::cerr << "ERROR: Error due to Message above. Exiting.." << std::endl;
//        return 0;
//    }



    //------------ Read Index & Reads
    std::cout << "Loading Index.." << std::endl;
    // load Index from disk
    Index < Dna5String, IndexSa<> > saIndex;
    open(saIndex, toCString("/home/phil/Dokumente/ALBIPraktikum/random10M.inx"));
    std::cout << "Index loaded!" << std::endl;



    //------------ read Reads-Records
    std::cout << "Reading Reads.." << std::endl;
    SeqFileIn rFileIn;
    if (!open(rFileIn, toCString("/home/phil/Dokumente/ALBIPraktikum/testRead.fasta"))) {
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



    //------------ Loop through reads
    typedef Iterator<StringSet<Dna5String> >::Type TIterator;
    for (TIterator readsIt = begin(reads); readsIt != end(reads); ++readsIt) {
        std::cout << "processing Element: " << position(readsIt, reads) << ".." << std::endl;



        // ------------ cut in small seeds
        Dna5String seed;
        // building Iterator
        unsigned pos = 0;
        typedef Iterator<Dna5String>::Type TIterator;
        // building best seed Vector & counting-Variables
        std::vector<int> bestSeed(3,0);
        int beginRefGen;
        int endRefGen;
        // Loop through seeds
        for (TIterator readIt = begin(*readsIt); readIt != end(*readsIt);) {
            // if rest of read is bigger than 10 letters
            if ((readIt + 10) <= end(*readsIt)) {
                seed = infix(*readIt, readIt, readIt + 10);
                readIt += 10;
                pos += 10;
                std::cout << "Searching for Seed: " << seed << std::endl;
                // Discard this part ?... could just be one letter... not good to search in Genome!
            } else { // if rest is not bigger than 10 letters
                //seed = infix(*readIt, readIt, end(*readsIt));
                readIt = end(*readsIt);
                //std::cout << "Searching for Seed: " << seed << std::endl;
                break;
            }



            // ------------ find seed
            Pattern<Dna5String> pat(seed);
            Finder<Index < Dna5String, IndexSa<> > > finder(saIndex);
            // reset finder
            clear(finder);
            // initialize Vector to store seedPosition in
            std::vector<int> begPos;
            std::vector<int> begEnd;
            while (find(finder, pat)) {
                begPos.push_back(beginPosition(finder));
                begEnd.push_back(endPosition(finder));
                //std::cout << '[' << beginPosition(finder) << '.' << endPosition(finder) << ']' << std::endl;
            }
            if (begPos.empty()) {
                std::cout << "nothing Found!.. next seed" << std::endl;
                break;
            }
            std::cout << "Found Seed " << begPos.size() <<" times!" << std::endl;




            // ------------ extend seed(s)
            Dna5String refGen;
            std::vector<int> matched;
            for (unsigned ind = 0; ind < begPos.size(); ind++) {
                // define seed
                Seed<Simple> seed1(10, pos - 10, 20, pos);
                // referenceGenome in Match Region with 10 characters additional to each side if possible
                beginRefGen = (begPos.at(ind) >= 10) ? begPos.at(ind) - 10 : 0;
                endRefGen = (begEnd.at(ind) + 10 < length(saIndex)) ? begEnd.at(ind) + 10 : length(saIndex);
                // get underlying Text of Index at position specified above
                refGen = infix(indexText(saIndex), beginRefGen, endRefGen);
                // match extension
                extendSeed(seed1, refGen, *readsIt, EXTEND_BOTH, MatchExtend());
                matched.push_back(endPositionH(seed1) - beginPositionH(seed1));
                //std::cout << "matched: " << matched.back() << std::endl;
            }// next seed to extend



            // ------------ save best seed
            std::vector<int>::iterator bestMatch = std::max_element(matched.begin(), matched.end());
            int posBestMatch = std::distance(matched.begin(), bestMatch);
            std::cout << "best match with " << *bestMatch << " Bases!" << std::endl;
            // save best seed if better than these before
            if (*bestMatch > bestSeed.back()) {
                bestSeed.at(0) = begPos.at(posBestMatch);
                bestSeed.at(1) = begEnd.at(posBestMatch);
                bestSeed.at(2) = *bestMatch;
            }

        } // next seed of read

        // ------------ perform semi-global Alignment for Read at position of best seed
        // referenceGenome in Match Region with 10 characters additional to each side if possible
        beginRefGen = (bestSeed.at(0) >= 60) ? bestSeed.at(0) - 60 : 0;
        endRefGen = (bestSeed.at(1) + 60 < length(saIndex)) ? bestSeed.at(1) + 60 : length(saIndex);
        typedef Align<Dna5String, ArrayGaps> TAlign;
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(indexText(saIndex), beginRefGen, endRefGen));
        assignSource(row(align, 1), *readsIt);
        // do alignment & get score
        int score = globalAlignment(align, Score<int, Simple>(1, -1, 0, -1));
        std::cout << align << std::endl;
        std::cout << "Score: " << score << std::endl;
        std::cout << "" << std::endl;


        // ------------ determine CIGAR Format
                // Use references to the rows of align.
        typedef Row<TAlign>::Type TRow;
        TRow & row1 = row(align, 0);
        TRow & row2 = row(align, 1);
        // Initialize the row iterators.
        typedef Iterator<TRow>::Type TRowIterator;
        TRowIterator itRow1 = begin(row1);
        TRowIterator itEndRow1 = end(row1);
        TRowIterator itRow2 = begin(row2);
        TRowIterator itEndRow2 = end(row2);
        std::stringstream cigar;
        int numChar = 0;
        int gapCount = 0;
        // don't count beginning Gaps!
        if (isGap(itRow1)||isGap(itRow2)){
            int numGaps = countGaps(itRow1);
            numGaps += countGaps(itRow2);
            itRow1 += numGaps;
            itRow2 += numGaps;
        }
        // don't count ending Gaps!
        int endGaps = 0;
        while(isGap(--itEndRow1)|| isGap(--itEndRow2)) {
            endGaps ++;
        }
        // count up last letter
        itEndRow1++;
        // get CIGAR String
        while ( itRow1 != itEndRow1){
            // Count insertions.
            if (isGap(itRow1)) {
                int numGaps = countGaps(itRow1);
                cigar << numGaps << "I";
                itRow1 += numGaps;
                itRow2 += numGaps;
                continue;
            }
            // Count deletions.
            if (isGap(itRow2)){
                int numGaps = countGaps(itRow2);
                cigar << numGaps << "D";
                itRow1 += numGaps;
                itRow2 += numGaps;
                continue;
            }
            // Count matches.
            while (*itRow1 == *itRow2 && itRow1 != itEndRow1)
            {
                ++numChar;
                ++itRow1;
                ++itRow2;
            }
            if (numChar != 0){
                cigar << numChar << "M";
                numChar = 0;
                continue;
            }
            // Count mismatches.
            while (*itRow1 != *itRow2 && itRow1 != itEndRow1)
            {
                ++numChar;
                ++itRow1;
                ++itRow2;
            }
            if (numChar != 0)
                cigar << numChar << "S";
            numChar = 0;
        }
        std::cout << "CIGAR Format: " << cigar.str() << std::endl;



        // ------------ write in SAM
        CharString bamFileName = "/home/phil/Dokumente/ALBIPraktikum/testBAM";
        //BamAlignmentRecord record;
        //BamFileOut bamFileOut();




    }// next read

    std::cout << " Finished";
    return 1;
}
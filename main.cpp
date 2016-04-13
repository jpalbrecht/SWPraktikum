#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/index.h>
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
    Index<Dna5String, IndexSa<> > saIndex(genome);
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
    Index<Dna5String, IndexSa<> > saIndex;
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
    typedef Iterator<StringSet<Dna5String>, Standard > ::Type TIterator;
    for (TIterator readsIt = begin(reads, Standard()); readsIt != end(reads, Standard()); ++readsIt) {
        std::cout << "processing Element: " << position(readsIt, reads) << ".." << std::endl;
        // ------------ cut in small seeds
        Dna5String seed;
        // building Iterator
        typedef Iterator<Dna5String, Standard>::Type TIterator;
        // Loop through seeds
        for (TIterator readIt = begin(*readsIt, Standard()); readIt != end(*readsIt, Standard());) {
            // if rest of read is bigger than 10 letters
            if ((readIt + 10) <= end(*readsIt, Standard())) {
                seed = infix(*readIt, readIt, readIt + 10);
                readIt += 10;
                std::cout << toCString(seed) << std::endl;
                // Discard this part ?... could just be one letter... not good to search in Genome!
            } else { // if rest is not bigger than 10 letters
                seed = infix(*readIt, readIt, end(*readsIt, Standard()));
                readIt = end(*readsIt, Standard());
                std::cout << toCString(seed) << std::endl;
            }
            // ------------ find seed
            Pattern<Dna5String> pat(seed);
            Finder<Index<Dna5String, IndexSa<> > > finder(saIndex);
            // reset finder
            clear(finder);
            while (find(finder, pat)) {
                std::cout << '[' << beginPosition(finder) << '.' << endPosition(finder) << ']' << std::endl;
            }
        }


    }

/*    Dna5String patDNA = "TAAGAG";
    Pattern<Dna5String> pat(patDNA);
    //setNeedle(pat, patDNA);
    Finder<Index<Dna5String, IndexSa<> > > finder(saIndex);
    // SA(suffixarray) ESA FM (DIES SIND TAGS) statt FAIIndex
    // INDEX AUCH AUF FESTPLATTE SPEICHERN! ANDERE APP LIEST EIN
    // indexrequire aufrufen. index wird erst gebaut wenn genutzt!
    // alles klar?
    clear(finder);
    while (find(finder, pat)) {
        std::cout << '[' << beginPosition(finder) << '.' << endPosition(finder) << ']'<< std::endl;
    }*/

    /*
    CharString meta;
    Dna5String dna;
    readRecord(meta, dna, sFileIn);
    */
    std::cout << " Finished";
    return 1;
}
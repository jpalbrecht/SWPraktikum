#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/stream.h>


using namespace seqan;





int main()
{

    // Read in FastaReferenceFile & Index
    SeqFileIn sFileIn;
    CharString filename = "/home/phil/Dokumente/ALBIPraktikum/random10M.fasta";
    //Index<Dna5String, IndexSa<> > saIndex;
    if (!open(sFileIn, toCString(filename))) {
        std::cout << "Could not read File!";
        return 1;
    }
    // read Reference-Sequence
    Dna5String genome;
    CharString head;
    readRecord(head, genome, sFileIn);
    // build Index
    Index<Dna5String, IndexSa<> > saIndex(genome);
    //indexCreate(saIndex, FibreSA());
    // require Index for being able to save to disk
    indexRequire(saIndex, FibreSA());

    // save Index to File!

    Dna5String patDNA = "TAAGAG";
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
    }

    /*
    CharString meta;
    Dna5String dna;
    readRecord(meta, dna, sFileIn);
    */
    std::cout << " Finished" ;
    return 0;
}
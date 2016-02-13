// ==========================================================================
//                                    Gustaf
// ==========================================================================
// Copyright (c) 2011-2013, Kathrin Trappe, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_
#define SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_


#include "../stellar/stellar.h"
#include "../stellar/stellar_output.h"

#include "msplazer.h"
#include "msplazer_out.h"
#include "msplazer_algorithms.h"
#include "gustaf_matepairs.h"
#include "stellar_routines.h"
#include "create_stellarmatches_from_file.h"
using namespace seqan;

// /////////////////////////////////////////////////////////////////////////////
// For benchmark
void process_mem_usage(double& vm_usage, double& resident_set)
{
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage     = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

// /////////////////////////////////////////////////////////////////////////////
// MSplazer Wrapper
int msplazer(StellarOptions & stellarOptions, MSplazerOptions & msplazerOptions)
{
    // Finder
    typedef String<Dna5> TSequence;
    //  typedef FragmentStore<void>::TReadSeqStore TSequence;

    // Database and query ID

    typedef CharString TId;

    // import query sequences using _importSequences from Stellar
    StringSet<TSequence> queries;
    StringSet<TId> queryIDs;
    String<unsigned> readJoinPositions;
    // TODO (ktrappe) distinguish between paired and single and, call appropriate
    // importSeq function (and preprocess query files)
    if (msplazerOptions.pairedEndMode)
    {
        std::cout << "Loading paired-end read sequences... ";
        if (!_importSequences(msplazerOptions.queryFile[0], msplazerOptions.queryFile[1], msplazerOptions.revCompl, queries, queryIDs, readJoinPositions))
            return 1;
    }else
    {
        std::cout << "Loading query sequences... ";
        if (!_importSequences(stellarOptions.queryFile, "query", queries, queryIDs))
            return 1;
    }
    StringSet<TId> shortQueryIDs = queryIDs;

    /*
    unsigned readLength = 0;
    for(unsigned i = 0; i < length(queries); ++i)
        readLength += length(queries[i]);
        */
    // std::cerr << "Loaded read seq: " << queries[i] << std::endl;
    std::cout << "done" << std::endl;

    // import database sequence using _importSequences from Stellar
    StringSet<TSequence> databases;
    StringSet<TId> databaseIDs;
    std::cout << "Loading reference sequences... ";
    if (!_importSequences(stellarOptions.databaseFile, "database", databases, databaseIDs))
        return 1;

    StringSet<TId> shortDatabaseIDs = databaseIDs;

    for (unsigned i = 0; i < length(databases); ++i)
    {
        std::cout << "Loaded db seq with length: " << length(databases[i]) << std::endl;
        std::cout << "Loaded db ID: " << databaseIDs[i] << std::endl;
    }

    std::cout << "done" << std::endl;


    // /////////////////////////////////////////////////////////////////////////
    // Compute Stellar Matches and their score
    // Note: Matches will be sorted when calling score function _getMatchDistanceScore
    // Note: Matches of the reverse strand will be modified in the sense that the match positions refer to the forward
    // strand (and not the reverse strand) get distance scores for all matches

    // Stellar Match Space Allocator
    typedef StringSet<QueryMatches<StellarMatch<TSequence, TId> > > TMatches;
    // FragmentStore
    /*
       TMatches matches;
       resize(matches, length(fragments.readSeqStore));
       */

    TMatches stellarMatches;
    TMatches testStMatches;
    resize(stellarMatches, length(queries));

    // Score space allocator
    typedef StellarMatch<TSequence, TId> TMatch;
    typedef String<int> TScoreAlloc;

    String<TScoreAlloc> distanceScores;
    resize(distanceScores, length(stellarMatches));


    // TODO (ktrappe) - a message from kjk
    // There are bugs in multithread mode. 
    // The result is not exactly same with "-nth -1".
    // Probably you know this situation. -eg. file reading.
    // So I disabled this for testing.
     unsigned int nth = msplazerOptions.numThreads;
     msplazerOptions.numThreads = 1;
 
    // get Stellar matches
    bool doStellar = true;
    if (msplazerOptions.stellarInputFile != "")
    {
        doStellar = false;
        std::cout << " Getting STELLAR matches from file, not calling STELLAR" << std::endl;
    }
    if (doStellar)
    {
        std::cout << "Calling STELLAR..." << std::endl;
        std::cout << "Stellar options:" << std::endl;
        _writeFileNames(stellarOptions);
        _writeSpecifiedParams(stellarOptions);
        _writeCalculatedParams(stellarOptions);
        _getStellarMatches(queries, databases, databaseIDs, stellarOptions, stellarMatches);
        std::cout << "done" << std::endl;
    }
    else
    {
        std::cout << "Importing STELLAR matches from file " << msplazerOptions.stellarInputFile << std::endl;
        double startST = sysTime();
        // TODO (ktrappe) distinguish call with queryIDs and shortQueryIDs in case of mate pairs? stellar writes out short
        // query IDs anyway...
       if (!_getStellarMatchesFromFile(queries, shortQueryIDs, databases, databaseIDs, msplazerOptions.stellarInputFile,
                                        stellarMatches, msplazerOptions.numThreads))
            return 1;

        std::cout << "done" << std::endl;
        std::cout << "TIME importing stellar matches " <<   (sysTime() - startST) << "s" << std::endl;
    }
    double vm, rss;
    process_mem_usage(vm, rss);
    std::cout << "Memory usage (after input file loading) : " << vm << "," << rss << std::endl;
    /*
    for(unsigned i = 0; i < length(queries); ++i){
        for(unsigned j = 0; j < length(stellarMatches[i].matches); ++j)
            std::cout << stellarMatches[i].matches[j] << std::endl;
    }
    */
    double startDist = sysTime();
    std::cout << "Getting match distance..." << std::endl;
    _getMatchDistanceScore(stellarMatches, distanceScores, msplazerOptions.numThreads);
    std::cout << "TIME getting match distance " <<   (sysTime() - startDist) << "s" << std::endl;


    // ///////////////////////////////////////////////////////////////////////
    // Create chains/graphs for all queries using Stellar matches


    // Graph structure for stellar matches
    typedef Graph<Directed<int> > TGraph; // TRowSize> > TGraph;
    // static_cast<Nothing>(TGraph());
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // Breakpoint property map
    typedef SparsePropertyMap<Breakpoint<TSequence, TId>, unsigned> TSparsePropertyMap;
    typedef String<TMatch> TMatchAlloc;

    // Container for graph/chain, scores and start and end vertices
    typedef MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TSparsePropertyMap, TMatchAlloc> TMSplazerChain;
    String<TMSplazerChain> queryChains;

    
    typedef Breakpoint<TSequence, TId> TBreakpoint;
    String<TBreakpoint> globalRawBreakpoints;
    String<TBreakpoint> globalBreakpoints;
    String<TBreakpoint> globalBreakends;

    // Run on parallel - minimally refactored, just for showing performance.
    omp_set_num_threads(nth);
    unsigned int jobSize = ceil(length(stellarMatches)/nth);

    String<String<TBreakpoint> > separatedBreakpoints;
    String<String<TBreakpoint> > separatedBreakends;
    resize(separatedBreakpoints, nth);
    resize(separatedBreakends, nth);

    std::cout << "Analyze with " << nth << " Threads : Start" << std::endl;
    double startTime = sysTime();
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned int i = 0; i < nth; ++i)
    {
        unsigned int startPos = i * jobSize;
        unsigned int endPos = i + (i+1)*jobSize;

        _analyze(stellarMatches, distanceScores, queryChains, queryIDs, queries, readJoinPositions, msplazerOptions, separatedBreakpoints[i], separatedBreakends[i], startPos, endPos);
    }
    omp_set_num_threads(msplazerOptions.numThreads); //disable again (msplazerOptions.numThreads=1)
    std::cout << "Analyze with " << nth << " Threads : End(" << (sysTime() - startTime) << "s)" << std::endl;

    process_mem_usage(vm, rss);
    std::cout << "Memory usage (before postprocessing) : " << vm << "," << rss << std::endl;

    std::cout << "Postprocessing : Start" << std::endl;
    // merge breakpoints
    unsigned int similarBPId = 0;
    for (unsigned int i = 0; i < nth; ++i)
    {
        _insertBreakpoints(globalRawBreakpoints, globalBreakends, separatedBreakpoints[i], msplazerOptions, similarBPId);
        //std::cout<< i << ":" << length(separatedBreakpoints[i]) << std::endl;
    }
    std::cout<< "merged global breakpoints (raw) :" << length(globalRawBreakpoints) << std::endl;

    // refine breakpoints
    similarBPId = 0;
    for (unsigned int i = 0; i < length(globalRawBreakpoints); ++i)
    {
        if (globalRawBreakpoints[i].support >= msplazerOptions.support)
        {
            if (msplazerOptions.inferComplexBP)
                _inferComplexBP(globalBreakpoints, globalRawBreakpoints[i], msplazerOptions.breakpointPosRange, similarBPId);
            else
                appendValue(globalBreakpoints, globalRawBreakpoints[i]);
        }
    }
    std::cout<< "merged global breakpoints (refined) :" << length(globalBreakpoints) << std::endl;

    // merge breakends
    for (unsigned int i = 0; i < nth; ++i)
    {
        append(globalBreakpoints, separatedBreakends[i]);
        //std::cout<< i << ":" << length(separatedBreakends[i]) << std::endl;
    }
    std::cout<< "merged global breakpoints (merged) :" << length(globalBreakpoints) << std::endl;

    // std sort in ascending order
    double startWriting = sysTime();
    std::cout << "Sorting and writing breakpoints... ";
    std::sort(begin(globalBreakpoints), end(globalBreakpoints));
    // std::sort(begin(globalStellarIndels), end(globalStellarIndels));
    _writeGlobalBreakpoints(globalBreakpoints, msplazerOptions, Gff());
    _writeGlobalBreakpoints(globalBreakpoints, databases, databaseIDs, msplazerOptions, Vcf());
    std::cout << "...done " << (sysTime() - startWriting) << "s" << std::endl;

    process_mem_usage(vm, rss);
    std::cout << "Memory usage (after postprocessing) : " << vm << "," << rss << std::endl;
 
    // ///////////////////////////////////////////////////////////////////////
    // Write dot files
    //
    // if(length(queries) < 200)
    if (msplazerOptions.dotOut)
        _writeDotfiles(stellarMatches, queries, queryIDs, queryChains, msplazerOptions);
    return 0;
}

#endif  // #ifndef SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_

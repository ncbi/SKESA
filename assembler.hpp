/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#ifndef _DBGAssembler_
#define _DBGAssembler_

#include <random>
#include "DBGraph.hpp"
#include "counter.hpp"
#include "graphdigger.hpp"

namespace DeBruijn {

    /******************************
    General description

    CDBGAssembler implements the SKESA assembling algorithm.

    1. It uses the counts for kmers with the minimal kmer length specified (default 21 bp) to estimate the maximal kmer length
       (starting from average mate length) that has sufficient coverage requested in maxkmercount. If reads are paired and
       insert size isn't specified, it estimates the insert size by assembling between mates for a sample of the reads.

    2. It assembles iteratively starting from minimal to maximal kmer length in a specified number of steps. Each step builds a
       de Bruijn graph for the kmer size for that iteration and uses it to improve previously assembled contigs. After each
       assembly iteration, the reads already used in the contigs are removed from further consideration.

    3. If reads are paired, it uses the reads that are not marked as used and the set of de Bruijn graphs built in 2) to connect
       the mate pairs. 

    4. Using the paired reads connected in 3), it performs three additional assembly iterations with the kmer size up
       to the insert size.
    *******************************/

    template<class DBGraph>
    class CDBGAssembler {
    public:
        // fraction - Maximal noise to signal ratio of counts acceptable for extension
        // jump - minimal length of accepted dead ends; i.e. dead ends shorter than this length are ignored
        // low_count - minimal count for kmers in a contig
        // steps - number of assembly iterations from minimal to maximal kmer size in reads
        // min_count - minimal kmer count to be included in a de Bruijn graph
        // min_kmer - the minimal kmer size for the main steps
        // max_kmer_paired - insert size (0 if not known)
        // maxkmercount - the minimal average count for estimating the maximal kmer
        // memory - the upper bound for memory use (GB)
        // ncores - number of threads
        // raw_reads - reads (for effective multithreading, number of elements in the list should be >= ncores)
        typedef typename DBGraph::Node Node;
        using GraphDigger = CDBGraphDigger<DBGraph>;        
        
        template<typename... GraphArgs>
        CDBGAssembler(double fraction, int jump, int low_count, int steps, int min_count, int min_kmer, int max_kmer, bool forcesinglereads,
                      int max_kmer_paired, int maxkmercount, int ncores, list<array<CReadHolder,2>>& raw_reads, TStrList seeds, 
                      bool allow_snps, bool estimate_min_count, GraphArgs... gargs) : 
            m_fraction(fraction), m_jump(jump), m_low_count(low_count), m_steps(steps), m_min_count(min_count), m_min_kmer(min_kmer), m_max_kmer(max_kmer),
            m_max_kmer_paired(max_kmer_paired), m_maxkmercount(maxkmercount), m_ncores(ncores), m_average_count(0), m_raw_reads(raw_reads) {

            m_insert_size = 0;

            for(auto& reads : m_raw_reads) {
                m_raw_pairs.push_back({reads[0], CReadHolder(false)});
            }    
            m_connected_reads.resize(m_raw_reads.size(), {CReadHolder(false), CReadHolder(true)});

            double total_seq = 0;
            size_t total_reads = 0;
            size_t paired = 0;
            for(auto& reads : m_raw_reads) {
                if(forcesinglereads) {
                    for(CReadHolder::string_iterator is = reads[0].sbegin(); is != reads[0].send(); ++is)
                        reads[1].PushBack(is);
                    reads[0].Clear();
                }
                total_seq += reads[0].TotalSeq()+reads[1].TotalSeq();
                total_reads += reads[0].ReadNum()+reads[1].ReadNum();
                paired += reads[0].ReadNum();
            }
            bool usepairedends = paired > 0;
 
            //graph for minimal kmer
            double average_count = GetGraph(m_min_kmer, m_raw_reads, true, estimate_min_count ? total_seq : 0, gargs...);
            if(average_count == 0)
                throw runtime_error("Reads are too short for selected minimal kmer length");
            m_average_count = average_count;

            // estimate genome
            int read_len = total_seq/total_reads+0.5;
            cerr << endl << "Average read length: " << read_len << endl;
            size_t genome_size = m_graphs[m_min_kmer]->GenomeSize();
            cerr << "Genome size estimate: " << genome_size << endl << endl;

            {// first iteration
                if(!seeds.empty()) {
                    m_contigs.push_back(TContigSequenceList());
                    for(string& seed : seeds) {
                        m_contigs.back().emplace_back();  // empty contig
                        auto& contig = m_contigs.back().back();
                        contig.InsertNewChunk();
                        contig.InsertNewVariant(); // one empty list
                        for(char c : seed) {
                            string ambigs = FromAmbiguousIUPAC[c];
                            if(ambigs.size() == 1) {
                                contig.ExtendTopVariant(c);
                            } else {
                                contig.InsertNewChunk();
                                for(char c : ambigs)
                                    contig.InsertNewVariant(c);
                                contig.InsertNewChunk();
                                contig.InsertNewVariant(); // one empty list
                            }
                        }
                    }
                    CombineSimilarContigs(m_contigs.back());
                    m_seeds = m_contigs.back();

                    cerr << "Seeds: " << m_contigs.back().size() << endl;

                    int num = 0;
                    for(auto& contig : m_contigs.back()) {
                        string first_variant;
                        for(auto& lst : contig)
                            first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                        cerr << ">Seed_" << ++num << endl << first_variant << endl;

                        int pos = 0;
                        for(unsigned chunk = 0; chunk < contig.size(); ++chunk) { //output variants
                            int chunk_len = contig[chunk].front().size();
                            if(contig.VariableChunk(chunk)) {
                                int left = 0;
                                if(chunk > 0)
                                    left = min(100,(int)contig[chunk-1].front().size());
                                int right = 0;
                                if(chunk < contig.size()-1)
                                    right = min(100,(int)contig[chunk+1].front().size());
                                int var = 0;
                                auto it = contig[chunk].begin();
                                for(++it; it != contig[chunk].end(); ++it) {
                                    auto& variant = *it;
                                    cerr << ">Variant_" << ++var << "_for_Seed_" << num << ":" << pos-left+1 << "_" << pos+chunk_len+right << "\n";
                                    if(chunk > 0) {
                                        for(int l = left ; l > 0; --l)
                                            cerr << *(contig[chunk-1].front().end()-l);
                                    }
                                    for(char c : variant)
                                        cerr << c;
                                    if(chunk < contig.size()-1) {
                                        for(int r = 0; r < right; ++r)
                                            cerr << contig[chunk+1].front()[r];
                                    }
                                    cerr << endl;
                                }
                            }
                            pos += chunk_len;
                        }
                    }
                }
                    
                ImproveContigs(m_min_kmer, false);
                if(m_contigs.back().empty())
                    throw runtime_error("Was not able to assemble anything");
            }

            //estimate max_kmer
            if(m_max_kmer == 0) {
                if(m_steps > 1 && average_count > m_maxkmercount) {
                    m_max_kmer = read_len+1-double(m_maxkmercount)/average_count*(read_len-min_kmer+1);
                    m_max_kmer = min(TKmer::MaxKmer(), m_max_kmer);
                    EstimateMaxKmer(read_len, gargs...);
                } else {
                    m_max_kmer = m_min_kmer;
                }
            }

            cerr << endl << "Average count: " << average_count << " Max kmer: " << m_max_kmer << endl;
            
            //estimate insert size
            if(steps > 1 || usepairedends) {
                if(m_max_kmer_paired == 0 && usepairedends) {
                    size_t mates = 0;
                    for(auto& rh : m_raw_reads)
                        mates += rh[0].ReadNum();
                    unsigned sample_size = 10000; // use 10000 reads for connecting to estimate insert size
                    unordered_set<size_t> selection;
                    if(mates/2 > 2*sample_size) {  // make random choice for reads
                        default_random_engine generator;
                        uniform_int_distribution<size_t> distribution(0,mates/2-1);
                        for(unsigned s = 0; s < sample_size; ) {
                            if(selection.insert(distribution(generator)).second) 
                                ++s;
                        }
                    } else if(mates/2 > 0) { // too few paired reads so using all : may be > sample_size but <= twice that size
                        for(size_t i = 0; i <= mates/2-1; ++i)
                            selection.insert(i);
                    }

                    if(!selection.empty()) {
                        CStopWatch timer;
                        timer.Restart();

                        list<array<CReadHolder,2>> mate_pairs;
                        size_t mp = 0;
                        int sub_sample = sample_size/m_ncores;
                        size_t num = 0;
                        for(auto& reads : m_raw_reads) {
                            for(CReadHolder::string_iterator is = reads[0].sbegin(); is != reads[0].send(); ++is, ++mp) {
                                if(selection.count(mp)) {
                                    if((num++)%sub_sample == 0)
                                        mate_pairs.push_back({CReadHolder(true), CReadHolder(false)});
                                    mate_pairs.back()[0].PushBack(is);
                                    mate_pairs.back()[0].PushBack(++is);
                                } else {
                                    ++is;
                                }
                            }
                        }
                
                        int long_insert_size = 2000; // we don't expect inserts to be longer than 2000 bp for this program
                        GraphDigger graph_digger(*m_graphs[min_kmer], m_fraction, m_jump, m_low_count);
                        list<array<CReadHolder,2>> connected_mate_pairs = graph_digger.ConnectPairs(mate_pairs, long_insert_size, m_ncores, false);
                        CReadHolder connected_mates(false);
                        for(auto& mp : connected_mate_pairs) {
                            for(CReadHolder::string_iterator is = mp[0].sbegin(); is != mp[0].send(); ++is)
                                connected_mates.PushBack(is);
                        }

 
                        m_max_kmer_paired = connected_mates.N50();
                        cerr << endl << "N50 for inserts: " << m_max_kmer_paired << endl << endl;

                    }
                }  
                m_max_kmer_paired = min(m_max_kmer_paired,TKmer::MaxKmer());
                m_insert_size = 3*m_max_kmer_paired; // we don't expect spread of histogram to go beyond three times expected insert

                CleanReads();               
            }
                
            //main iterations
            if(m_steps > 1) {
                if(m_max_kmer > 1.5*m_min_kmer) {
                    double alpha = double(m_max_kmer-m_min_kmer)/(steps-1); // find desired distance between consecutive kmers
                    for(int step = 1; step < m_steps; ++step) {
                        int kmer_len = min_kmer+step*alpha+0.5;             // round to integer
                        kmer_len -= 1-kmer_len%2;                           // get odd kmer
                        if(kmer_len <= m_graphs.rbegin()->first)
                            continue;
                        if(GetGraph(kmer_len, m_raw_reads, true, 0, gargs...) == 0) {
                            cerr << "Empty graph for kmer length: " << kmer_len << " skipping this and longer kmers" << endl;
                            break;
                        }
                        ImproveContigs(kmer_len, false);
                        CleanReads();
                    }
                } else {
                    cerr << "WARNING: iterations are disabled" << endl;
                }
            }
            
            // three additional iterations with kmers (usually) longer than read length and upto insert size
            if(usepairedends && m_insert_size > 0 && m_max_kmer_paired > 1.5*m_max_kmer) {
                ConnectPairsIteratively();

                array<int,3> long_kmers;
                long_kmers[0] = 1.25*m_max_kmer;
                long_kmers[2] = m_max_kmer_paired;
                long_kmers[1] = (long_kmers[0]+long_kmers[2])/2;
                    
                for(int kmer_len : long_kmers) {
                    kmer_len -= 1-kmer_len%2;
                    if(GetGraph(kmer_len, m_connected_reads, false, 0, gargs...) == 0) {
                        cerr << "Empty graph for kmer length: " << kmer_len << " skipping this and longer kmers" << endl;
                        break;
                    }
                    ImproveContigs(kmer_len, false);
                }
            } 

            if(allow_snps) { // snp discovery
                for(auto it = m_graphs.rbegin(); it != m_graphs.rend(); ++it) {
                    int kmer_len = it->first;
                    ImproveContigs (kmer_len, true);
                }
            }
        }        

        map<int,DBGraph*>& Graphs() { return m_graphs; }
        TContigSequenceList& Contigs() { return m_contigs.back(); }
        vector<TContigSequenceList>& AllIterations() { return m_contigs; }
        CReadHolder ConnectedReads() const {
            CReadHolder connected_reads(false);
            for(const auto& cr : m_connected_reads) {
                for(CReadHolder::string_iterator is = cr[0].sbegin(); is != cr[0].send(); ++is)
                    connected_reads.PushBack(is);
            }
            return connected_reads;
        }

        virtual ~CDBGAssembler() {
            for(auto& graph : m_graphs)
                delete graph.second;    
        }

    private:
        // connects paired reads using all constructed de Bruijn graphs 
        void ConnectPairsIteratively() {
            for(auto& gr : m_graphs) {
                int kmer_len = gr.first;
                cerr << endl << "Connecting mate pairs using kmer length: " << kmer_len << endl;
                GraphDigger graph_digger(*gr.second, m_fraction, m_jump, m_low_count);
                list<array<CReadHolder,2>> connected_reads_temp = graph_digger.ConnectPairs(m_raw_pairs, m_insert_size, m_ncores, true);
                list<array<CReadHolder,2>>::iterator pairedi = m_connected_reads.begin();
                list<array<CReadHolder,2>>::iterator rawi = m_raw_pairs.begin();
                for(auto& pr : connected_reads_temp) {
                    swap((*rawi)[0], pr[1]);                                                         // keep still not connected            
                    for(CReadHolder::string_iterator is = pr[0].sbegin(); is != pr[0].send(); ++is)  // add new connected reads             
                        (*pairedi)[0].PushBack(*is);
                    ++rawi;
                    ++pairedi;
                }
            }

            size_t connected = 0;
            for(auto& rh : m_connected_reads)
                connected += rh[0].ReadNum();        
            cerr << "Totally connected: " << connected << endl;

            size_t added = 0;
            list<array<CReadHolder,2>>::iterator pairedi = m_connected_reads.begin();
            for(auto& reads : m_raw_pairs) {
                for(CReadHolder::string_iterator is = reads[0].sbegin(); is != reads[0].send(); ++is) {
                    if((int)is.ReadLen() > m_max_kmer) {
                        (*pairedi)[0].PushBack(*is);
                        ++added;
                    }
                }
                ++pairedi;
            }
            cerr << "Added notconnected: " << added << endl;
        }

        // scans kmers for all assembled contigs and creates a map 
        // the key is the smaller of two possible kmer directions
        // the value is a tupe:
        //     int -     position on contig
        //     bool -    the same as the key or reverse complemented
        //     CContigSequence* - pointer to the contig
        typedef CKmerHashMap<tuple<int, bool, const CContigSequence*>, 8> TKmerToContig;
        //        typedef CKmerMap<tuple<int, bool, const CContigSequence*>> TKmerToContig;           
        TKmerToContig GetAssembledKmers() {
            int kmer_len = m_graphs.rbegin()->first;

            CKmerMap<int> seed_kmers(kmer_len);
            for(auto& seed : m_seeds) {
                if((int)seed.LenMin() < kmer_len)
                    continue;
                seed.RemoveShortUniqIntervals(kmer_len);

                for(int i = seed.size()-1; i >= 0; i -= 2) {
                    if(i == (int)seed.size()-1) {
                        if((int)seed.ChunkLenMax(i) >= kmer_len) { // last chunk could be short
                            CReadHolder rh(false);
                            rh.PushBack(seed.back().front());
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) {
                                TKmer kmer = *ik;
                                TKmer rkmer = revcomp(kmer, kmer_len);
                                ++seed_kmers[kmer < rkmer ? kmer : rkmer];
                            }
                        }
                    } else { // all uniq chunks in the middle >= kmer_len; first/last could be short
                        if((int)seed.ChunkLenMax(i) >= kmer_len) {
                            TVariation seq(seed[i].front().begin(), seed[i].front().end());
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) {
                                TKmer kmer = *ik;
                                TKmer rkmer = revcomp(kmer, kmer_len);
                                ++seed_kmers[kmer < rkmer ? kmer : rkmer];
                            }
                        }
                        for(auto& variant : seed[i+1]) {
                            TVariation seq;
                            if((int)seed.ChunkLenMax(i) >= kmer_len-1)
                                seq.insert(seq.end(), seed[i].front().end()-kmer_len+1, seed[i].front().end());
                            else
                                seq.insert(seq.end(), seed[i].front().begin(), seed[i].front().end());
                            seq.insert(seq.end(), variant.begin(), variant.end());
                            if((int)seed.ChunkLenMax(i+2) >= kmer_len-1)
                                seq.insert(seq.end(), seed[i+2].front().begin(), seed[i+2].front().begin()+kmer_len-1);
                            else
                                seq.insert(seq.end(), seed[i+2].front().begin(), seed[i+2].front().end());
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) {
                                TKmer kmer = *ik;
                                TKmer rkmer = revcomp(kmer, kmer_len);
                                ++seed_kmers[kmer < rkmer ? kmer : rkmer];
                            }
                        }
                    }
                }
            }

            cerr << "Seed kmers: " << seed_kmers.Size() << endl;

            int min_len = max(m_max_kmer_paired, m_max_kmer);
            size_t knum = 0;
            list<pair<CContigSequence*, SAtomic<int8_t>>> contigs;
            for(auto& contig : m_contigs.back()) {
                if((int)contig.LenMin() >= min_len && contig.size() == 1) {
                    contigs.emplace_back(&contig, 0);
                    knum += contig.LenMin()+2*(kmer_len-1);  // overestimation for reserve
                }
            }
            TKmerToContig assembled_kmers(kmer_len, knum);
            
            list<function<void()>> jobs;
            for(int thr = 0; thr < m_ncores; ++thr)
                jobs.push_back(bind(&CDBGAssembler::AssembledKmersJob, this, ref(contigs), ref(assembled_kmers), ref(seed_kmers)));
            RunThreads(m_ncores, jobs);

            return assembled_kmers;
        }

        void AssembledKmersJob(list<pair<CContigSequence*, SAtomic<int8_t>>>& contigs, TKmerToContig& assembled_kmers, CKmerMap<int>& seed_kmers) const {
            for(auto& pr : contigs) {
                if(!pr.second.Set(1))
                    continue;
                auto& contig = *pr.first;
                int kmer_len = m_graphs.rbegin()->first;
                auto& graphp = m_graphs.rbegin()->second;
                int pos = contig.ChunkLenMax(0)-kmer_len;
                CReadHolder rh(false);
                if(contig.m_circular) {
                    auto cc = contig[0].front();
                    cc.insert(cc.end(), contig[0].front().begin(), contig[0].front().begin()+kmer_len-1); // add kmer-1 bases to get all kmers
                    rh.PushBack(cc);
                    pos = contig.ChunkLenMax(0)-1;
                } else {
                    rh.PushBack(contig[0].front());
                }
                bool found_repeat = false;
                list<pair<TKmer, tuple<int, bool, const CContigSequence*>>> contig_kmers;
                for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) { // iteration from last kmer to first
                    TKmer kmer = *ik;
                    auto node = graphp->GetNode(kmer);
                    if(graphp->Abundance(node)*m_fraction > m_average_count)
                        continue;
                    if(node.isValid() && graphp->IsMultContig(node)) {
                        found_repeat = true;
                        break;
                    }
                    TKmer rkmer = revcomp(kmer, kmer_len);
                    TKmer* kmerp = &kmer;
                    bool direct = true;
                    if(rkmer < kmer) {
                        kmerp = &rkmer;
                        direct = false;
                    }
                    contig_kmers.emplace_back(*kmerp, make_tuple(pos, direct, &contig));
                }
                if(!found_repeat) {
                    for(auto& kmer :  contig_kmers) {
                        if(seed_kmers.Find(kmer.first) == nullptr)
                            *assembled_kmers.FindOrInsert(kmer.first) = kmer.second;
                    }
                }
            }
        }

        // finds if a read belongs to any of the contigs
        // return tuple:
        //   int -     position on the contig (-1 if not found)
        //   int -     +1 if in positive strand; -1 if in negative strand
        //   CContigSequence* - pointer to the contig 
        static tuple<int, int, const CContigSequence*> FindMatchForRead(const CReadHolder::string_iterator& is, TKmerToContig& assembled_kmers) {
            int rlen = is.ReadLen();
            int kmer_len = assembled_kmers.KmerLen();

            int plus = 1;
            tuple<int, bool, const CContigSequence*>* rsltp = nullptr;
            int knum = rlen-kmer_len+1;
            for(CReadHolder::kmer_iterator ik = is.KmersForRead(kmer_len); rsltp == nullptr && knum > 0; --knum, ++ik) {
                TKmer kmer = *ik;
                TKmer rkmer = revcomp(kmer, kmer_len);
                TKmer* kmerp = &kmer;
                plus = 1;
                if(rkmer < kmer) {
                    kmerp = &rkmer;
                    plus = -plus;
                }
                rsltp = assembled_kmers.Find(*kmerp);
                if(rsltp != nullptr && get<0>(*rsltp) < 0)
                    rsltp = nullptr;
            }

            int pos = -1; // position on contig of the 'outer' read end (aka insert end)    
            const CContigSequence* sp = nullptr;
            if(rsltp != nullptr) {
                sp = get<2>(*rsltp); // pointer to the contig
                if(!get<1>(*rsltp))
                    plus = -plus;
                if(plus > 0) {
                    pos = get<0>(*rsltp)-knum;
                    if(pos < 0 && sp->m_circular)
                        pos += sp->LenMax();
                } else {
                    pos = get<0>(*rsltp)+kmer_len-1+knum;
                    if(pos >= (int)sp->LenMax() && sp->m_circular)
                        pos -= sp->LenMax();
                }
            }
            
            return make_tuple(pos, plus, sp);
        }
    
        // removes reads if they belong to already assembled contigs
        // using contig sequence creates artificial connected pairs when both mates are placed
        //
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // raw_reads - reads
        // connected_reads - pointer to connected reads (nullp if not used)
        static void RemoveUsedReadsJob(TKmerToContig& assembled_kmers, int margin, int insert_size, array<CReadHolder,2>& raw_reads, CReadHolder* connected_reads) {
            int kmer_len = assembled_kmers.KmerLen();

            {
                CReadHolder cleaned_reads(true);
                if(raw_reads[0].ReadNum() > 0)
                    cleaned_reads.Reserve(raw_reads[0].TotalSeq(), raw_reads[0].ReadNum());

                CReadHolder::string_iterator is1 = raw_reads[0].sbegin();
                CReadHolder::string_iterator is2 = raw_reads[0].sbegin();
                ++is2;
                for( ; is2 != raw_reads[0].send(); ++is1, ++is1, ++is2, ++is2) {
                    if((int)min(is1.ReadLen(), is2.ReadLen()) < kmer_len) {
                        if(connected_reads) {    // keep short pairs for connection         
                            cleaned_reads.PushBack(is1);
                            cleaned_reads.PushBack(is2);                                 
                        } else {                 // give chance to be used as unpaired  
                            raw_reads[1].PushBack(is1);
                            raw_reads[1].PushBack(is2);
                        }
                        continue;
                    }

                    tuple<int, int, const CContigSequence*> rslt1 = FindMatchForRead(is1, assembled_kmers);
                    int pos1 = get<0>(rslt1);
                    int plus1 = get<1>(rslt1);
                    const CContigSequence* sp1 = get<2>(rslt1);
                    int clen1 = 0;
                    int left_flank1 = 0;
                    int right_flank1 = 0;
                    if(pos1 >= 0) {
                        left_flank1 = sp1->m_left_repeat;
                        right_flank1 = sp1->m_right_repeat;
                        clen1 = sp1->LenMax();
                        if(sp1->m_circular || (plus1 > 0 && pos1 >= margin+left_flank1 && pos1+insert_size-1 < clen1-margin-right_flank1) || 
                           (plus1 < 0 && pos1-insert_size+1 >= margin+left_flank1 && pos1 < clen1-margin-right_flank1))
                            continue;
                    }
                    // check for second mate in case first mate was of bad quality and not found in contigs 
                    tuple<int, int, const CContigSequence*> rslt2 = FindMatchForRead(is2, assembled_kmers);
                    int pos2 = get<0>(rslt2);
                    int plus2 = get<1>(rslt2);
                    const CContigSequence* sp2 = get<2>(rslt2);
                    if(pos2 >= 0) {
                        int left_flank2 = sp2->m_left_repeat;
                        int right_flank2 = sp2->m_right_repeat;
                        int clen2 = sp2->LenMax();
                        if(sp2->m_circular || (plus2 > 0 && pos2 >= margin+left_flank2 && pos2+insert_size-1 < clen2-margin-right_flank2) || 
                           (plus2 < 0 && pos2-insert_size+1 >= margin+left_flank2 && pos2 < clen2-margin-right_flank2))
                            continue;
                    }
                    if(pos1 >= 0 && pos2 >= 0 && sp1 == sp2 && plus1 != plus2) { // same contig, different strands      
                        if((plus1 > 0 && pos1 >= margin+left_flank1 && pos2 < clen1-margin-right_flank1) || 
                           (plus1 < 0 && pos2 >= margin+left_flank1 && pos1 < clen1-margin-right_flank1)) { // deep inside     
                            continue; 
                        } else if(connected_reads) {
                            if((plus1 > 0 && pos1 >= 0 && pos2 < clen1) || (plus1 < 0 && pos2 >= 0 && pos1 < clen1)) {                     // inside but not deep     
                                int a = min(pos1,pos2);
                                int b = max(pos1,pos2);
                                if(b < (int)sp1->ChunkLenMax(0)) {              // in first uniq chunk
                                    TVariation seq(sp1->front().front().begin()+a, sp1->front().front().begin()+b+1); 
                                    connected_reads->PushBack(seq);
                                    continue;
                                } else if(clen1-a <= (int)sp1->ChunkLenMax(sp1->size()-1)) { // in last uniq chunk
                                    TVariation seq(sp1->back().front().end()-clen1+a, sp1->back().front().end()-clen1+b+1);
                                    connected_reads->PushBack(seq);
                                    continue;
                                }
                            }
                        }
                    }
                    cleaned_reads.PushBack(is1);
                    cleaned_reads.PushBack(is2);                                 
                }
                cleaned_reads.Swap(raw_reads[0]);
            }

            if(!connected_reads) {          
                CReadHolder cleaned_reads(false);
                if(raw_reads[1].ReadNum() > 0)
                    cleaned_reads.Reserve(raw_reads[1].TotalSeq(), raw_reads[1].ReadNum());
                for(CReadHolder::string_iterator is = raw_reads[1].sbegin() ;is != raw_reads[1].send(); ++is) {
                    int rlen = is.ReadLen();
                    if(rlen < kmer_len)
                        continue;        
            
                    tuple<int, int, const CContigSequence*> rslt = FindMatchForRead(is, assembled_kmers);
                    int pos = get<0>(rslt);
                    int plus = get<1>(rslt);
                    const CContigSequence* sp = get<2>(rslt);
                    if(pos >= 0) {
                        int left_flank = sp->m_left_repeat;
                        int right_flank = sp->m_right_repeat;
                        int clen = sp->LenMax();
                        if(sp->m_circular || (plus > 0 && pos >= margin+left_flank && pos+rlen-1 < clen-margin-right_flank) || 
                           (plus < 0 && pos-rlen+1 >= margin+left_flank && pos < clen-margin-right_flank))
                            continue;
                    }

                    cleaned_reads.PushBack(is);
                }            
                cleaned_reads.Swap(raw_reads[1]);
            }
        }

        // removes used reads from the read set used for de Bruijn graphs
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // ncores - number of threads
        // raw_reads - reads
        static void RemoveUsedReads(TKmerToContig& assembled_kmers, int margin, int insert_size, int ncores, list<array<CReadHolder,2>>& raw_reads) {
            list<function<void()>> jobs;
            for(auto& job_input : raw_reads) {
                jobs.push_back(bind(RemoveUsedReadsJob, ref(assembled_kmers), margin, insert_size, ref(job_input), (CReadHolder*)0));                
            }
            RunThreads(ncores, jobs);
        }

        // removes used reads from the read set used for pair connection and from already connected (by contig sequence) reads
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // ncores - number of threads
        // raw_reads - reads
        // connected_reads - already connected by contig sequence reads
        static void RemoveUsedPairs(TKmerToContig& assembled_kmers, int margin, int insert_size, int ncores, list<array<CReadHolder,2>>& raw_reads, list<array<CReadHolder,2>>& connected_reads) {
            list<function<void()>> jobs;
            auto icr = connected_reads.begin();
            for(auto& job_input : raw_reads) {
                jobs.push_back(bind(RemoveUsedReadsJob, ref(assembled_kmers), margin, insert_size, ref(job_input), &(*icr++)[1]));                
            }
            RunThreads(ncores, jobs);
        }

        // removes used reads from the read set used for de Bruijn graphs and from the read set used for pair connection 
        // removes paired reads not needed as they are already connected by contig sequence reads
        // creates new set of reads to use
        void CleanReads() {
            CStopWatch timer;
            timer.Restart();
            TKmerToContig assembled_kmers = GetAssembledKmers();

            if(assembled_kmers.TableSize() > 0) {
                int jump = 50; //TODO reconsile with what used in filterneighbors   
                RemoveUsedReads(assembled_kmers, m_max_kmer+jump, m_insert_size, m_ncores, m_raw_reads);
                RemoveUsedReads(assembled_kmers, jump, m_insert_size, m_ncores, m_connected_reads);
                RemoveUsedPairs(assembled_kmers, jump, m_insert_size, m_ncores, m_raw_pairs, m_connected_reads);
            }

            size_t reads = 0;
            for(auto& rh : m_raw_reads)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Cleaned reads: " << reads << endl;
            reads = 0;
            for(auto& rh : m_raw_pairs)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Reads for connection: " << reads << endl;
            reads = 0;
            for(auto& rh : m_connected_reads)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Internal reads: " << reads << endl;
            cerr << "Reads cleaned in " << timer.Elapsed();                                               
        }

        // improves previously assembled contigs using a longer kmer
        void ImproveContigs (int kmer_len, bool allow_snps) {
            DBGraph& graph = *m_graphs[kmer_len];
            int jump = m_jump;
            if(allow_snps)
                jump += kmer_len;
            GraphDigger graph_digger(graph, m_fraction, jump, m_low_count, allow_snps);
            cerr << "Kmer: " << kmer_len << " Graph size: " << graph.GraphSize() << " Contigs in: " << (m_contigs.empty() ? 0 : m_contigs.back().size()) << endl;
            cerr << "Valley: " << graph_digger.HistMin() << endl; 
            CStopWatch total;
            total.Restart();
            
            CStopWatch timer;
            timer.Restart();
            //convert strings to SContig and mark visited kmers 
            if(allow_snps)
                graph.ClearAllVisited();  
              
            TContigList<DBGraph> scontigs = ConverToSContigAndMarkVisited(graph_digger); 
            cerr << endl << "Mark used kmers in " << timer.Elapsed();

            if(allow_snps)
                graph_digger.CheckRepeats(scontigs);

            size_t singl = 0;
            size_t multipl = 0;
            for(auto it = graph.Begin(); it != graph.End(); ++it) {
                if(graph.IsMultContig(it))
                    ++multipl;
                else if(graph.IsVisited(it))
                    ++singl;
            }
            cerr << "Kmers in multiple/single contigs: " << multipl << " " << singl << endl;            


            // connect overlapping contigs if we had seeds 
            if(!m_seeds.empty() && !allow_snps) {
                timer.Restart();
                graph_digger.CheckRepeats(scontigs);
                cerr << "Check repeats in " << timer.Elapsed();
                timer.Restart();
                graph_digger.ConnectOverlappingContigs(scontigs);
                cerr << "Connect overlapping contigs in " << timer.Elapsed();
            }

            timer.Restart();
            //create new contigs using not yet included kmers
            GraphDigger graph_digger_no_jump(graph, m_fraction, 0, m_low_count);
            unsigned min_len_for_new_seeds = 3*kmer_len;      // short ones are likely to be noise

            GraphDigger test_graphdigger(*m_graphs[m_min_kmer], m_fraction, 0, m_low_count);
            GraphDigger* test_graphdiggerp = nullptr;
            if(kmer_len != m_min_kmer)
                test_graphdiggerp = &test_graphdigger;

            TContigList<DBGraph> new_seeds = graph_digger_no_jump.GenerateNewSeeds(min_len_for_new_seeds, m_ncores, test_graphdiggerp);
            cerr << "New seeds: " << new_seeds.size() << endl;
            //add new seeds
            scontigs.splice(scontigs.end(), new_seeds);            
            cerr << "New seeds in " << timer.Elapsed();

            timer.Restart();
            graph_digger.ConnectAndExtendContigs(scontigs, m_ncores);             

            // convert back to CContigSequence 
            m_contigs.push_back(TContigSequenceList());
            for(auto& contig : scontigs) {
                m_contigs.back().push_back(contig.m_seq);
            }
            m_contigs.back().sort();

            vector<size_t> contigs_len;
            size_t genome_len = 0;
            for(auto& contig : m_contigs.back()) {
                contigs_len.push_back(contig.LenMax());
                genome_len += contigs_len.back();                            
            }
            sort(contigs_len.begin(), contigs_len.end());
            size_t n50 = 0;
            int l50 = 0;
            size_t len = 0;
            for(int j = (int)contigs_len.size()-1; j >= 0 && len < 0.5*genome_len; --j) {
                ++l50;
                n50 = contigs_len[j];
                len += contigs_len[j];
            }
            cerr << "Connections and extensions in " << timer.Elapsed();

            cerr << "Contigs out: " << contigs_len.size() << " Genome: " << genome_len << " N50: " << n50 << " L50: " << l50 << endl; 
            cerr << "Assembled in " << total.Elapsed() << endl; 
        }


        // converts contigs from the previous iteration into SContig and marks visited the nodes in the graph
        TContigList<DBGraph> ConverToSContigAndMarkVisited(GraphDigger& graph_digger) {
            if(m_contigs.empty())
                return TContigList<DBGraph>();

            int kmer_len = graph_digger.Graph().KmerLen();

            for(auto& contig : m_contigs.back()) {

                //remove short snps
                if(!contig.m_circular) {
                    if(contig.size() > 1 && (int)contig.ChunkLenMax(0) < kmer_len) {
                        contig.m_left_repeat = 0;
                        contig.erase(contig.begin(), contig.begin()+2);
                    }
                    if(contig.size() > 1 && (int)contig.ChunkLenMax(contig.size()-1) < kmer_len) {
                        contig.m_right_repeat = 0;
                        contig.pop_back();
                        contig.pop_back();
                    }
                }
            }

            TContigList<DBGraph> scontigs;
            vector<pair<const CContigSequence*, SAtomic<uint8_t>>> contig_is_taken;
            for(const auto& contig : m_contigs.back())
                contig_is_taken.push_back(make_pair(&contig,SAtomic<uint8_t>(0)));
            vector<TContigList<DBGraph>> scontigs_for_threads(m_ncores);
            list<function<void()>> jobs;
            for(auto& sc : scontigs_for_threads) 
                jobs.push_back(bind(&CDBGAssembler::ConverToSContigAndMarkVisitedJob, this, ref(contig_is_taken), ref(sc), ref(graph_digger)));
            RunThreads(m_ncores, jobs);

            for(auto& sc : scontigs_for_threads)
                scontigs.splice(scontigs.end(), sc);
    
            return scontigs;
        }

        // one-thread worker for ConverToSContigAndMarkVisited()
        void ConverToSContigAndMarkVisitedJob(vector<pair<const CContigSequence*, SAtomic<uint8_t>>>& contig_is_taken, TContigList<DBGraph>& scontigs, GraphDigger& graph_digger) {
            DBGraph& graph = graph_digger.Graph();
            int kmer_len = graph.KmerLen();
            for(auto& cnt : contig_is_taken) {
                if(!cnt.second.Set(1))
                    continue;
                const CContigSequence& contig = *cnt.first;
                int contig_len = contig.LenMin();
                if(contig_len >= kmer_len)
                    scontigs.push_back(SContig<DBGraph>(contig, graph)); // constructor sets visited in graph                    
            }
        }

        // estimates available memory
        int64_t AvailableMemory(int memory) const {
            int64_t GB = 1000000000;
            int64_t mem_available = GB*memory;
            int64_t mem_used = 0;
            for(const auto& reads : m_raw_reads)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(const auto& reads : m_raw_pairs)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(const auto& reads : m_connected_reads)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(auto& graph : m_graphs)
                mem_used += graph.second->MemoryFootprint();
            for(auto& lst : m_contigs) {
                for(auto& contig : lst)
                    mem_used += contig.MemoryFootprint()+2*sizeof(CContigSequence*);  // contig and 2 list pointers
            }

            return mem_available-mem_used;

        }

        template<typename... GraphArgs>
        void EstimateMaxKmer(int read_len, GraphArgs... gargs) {
            static_assert(sizeof(DBGraph) != sizeof(DBGraph), "Unknown specialization of CDBGAssembler");
        }

        // counts kmers and build a de Bruijn graph; returns average count of kmers in the graph
        // kmer_len - the size of the kmer
        // reads - reads from input or connected internally
        // is_stranded - whether or not stranded information is meaningful
        template<typename... GraphArgs>
        double GetGraph(int kmer_len, const list<array<CReadHolder,2>>& reads, bool is_stranded, double total_seq, GraphArgs... gargs) {
            static_assert(sizeof(DBGraph) != sizeof(DBGraph), "Unknown specialization of CDBGAssembler");
            return 0;
        }


        double m_fraction;                                   // Maximal noise to signal ratio of counts acceptable for extension
        int m_jump;                                          // minimal length of accepted dead ends
        int m_low_count;                                     // minimal kmer count to be included in a contig
        int m_steps;                                         // number of main steps
        int m_min_count;                                     // minimal kmer count to be included in a de Bruijn graph
        int m_min_kmer;                                      // the minimal kmer size for the main steps
        int m_max_kmer;                                      // maximal kmer size for the main steps
        int m_max_kmer_paired;                               // insert size
        int m_insert_size;                                   // upper bound for the insert size
        int m_maxkmercount;                                  // the minimal average count for estimating the maximal kmer
        int m_ncores;                                        // number of threads

        double m_average_count;                              // average count for minimal kmers

        list<array<CReadHolder,2>>& m_raw_reads;             // original reads - will be reduced gradually
        list<array<CReadHolder,2>> m_raw_pairs;              // paired original reads for connection - will be reduced gradually
        list<array<CReadHolder,2>> m_connected_reads;        // connected pairs (long reads)
        map<int, DBGraph*> m_graphs;                         // De Bruijn graphs for mutiple kmers
        vector<TContigSequenceList> m_contigs;               // assembled contigs for each iteration
        TContigSequenceList m_seeds;
    };

    template<> template<> // one for graph, the other for args
    void CDBGAssembler<CDBGraph>::EstimateMaxKmer(int read_len, int memory) {
        while(m_max_kmer > m_min_kmer) {
            m_max_kmer -= 1-m_max_kmer%2;           // odd kmers desired
            CKmerCounter kmer_counter(m_raw_reads, m_max_kmer, m_min_count, true, AvailableMemory(memory), m_ncores);
            if(kmer_counter.Kmers().Size() < 100) { // find a kmer length with at least 100 distinct kmers at that length
                m_max_kmer -= read_len/25;          // reduce maximal kmer length by a small amount based on read length
                continue;
            }
            double average_count_for_max_kmer = kmer_counter.AverageCount();
            if(average_count_for_max_kmer >= m_maxkmercount)
                break;
            else 
                m_max_kmer -= read_len/25;                                    
        }
        m_max_kmer = max(m_max_kmer, m_min_kmer);
    }

    template<> template<> // one for graph, the other for args
    void CDBGAssembler<CDBHashGraph>::EstimateMaxKmer(int read_len, int estimated_kmer_num, bool skip_bloom_filter) {
        int64_t M = 1000000;
        while(m_max_kmer > m_min_kmer) {
            m_max_kmer -= 1-m_max_kmer%2;           // odd kmers desired
            CKmerHashCounter  kmer_counter(m_raw_reads, m_max_kmer, m_min_count, M*estimated_kmer_num, true, m_ncores, skip_bloom_filter);
            if(kmer_counter.KmerNum() < 100) {      // find a kmer length with at least 100 distinct kmers at that length
                m_max_kmer -= read_len/25;          // reduce maximal kmer length by a small amount based on read length
                continue;
            }
            double average_count_for_max_kmer = GetAverageCount(kmer_counter.Kmers().GetBins());           
            if(average_count_for_max_kmer >= m_maxkmercount)
                break;
            else 
                m_max_kmer -= read_len/25;                                    
        }
        m_max_kmer = max(m_max_kmer, m_min_kmer);
    }

    template<> template<> // one for graph, the other for args
    double CDBGAssembler<CDBGraph>::GetGraph(int kmer_len, const list<array<CReadHolder,2>>& reads, bool is_stranded, double total_seq, int memory) {
        CKmerCounter kmer_counter(reads, kmer_len, m_min_count, is_stranded, AvailableMemory(memory), m_ncores);
        if(kmer_counter.Kmers().Size() == 0)
            return 0;

        TKmerCount& sorted_kmers =  kmer_counter.Kmers();

        if(total_seq > 0) {
            map<int,size_t> hist;
            for(size_t index = 0; index < sorted_kmers.Size(); ++index) {
                ++hist[sorted_kmers.GetCount(index)];                  // count clipped to integer automatically
            }
            TBins bins(hist.begin(), hist.end());
            int genome_size = CalculateGenomeSize(bins);
            if(genome_size > 0) {
                int new_min_count = total_seq/genome_size/50+0.5;
                if(new_min_count > m_min_count) {
                    int new_maxkmercount = max(10, int(total_seq/genome_size/10+0.5));
                    cerr << "WARNING: --min_count changed from " << m_min_count << " to " << new_min_count << " because of high coverage for genome size " << genome_size << endl;
                    cerr << "WARNING: --max_kmer_count " << m_maxkmercount << " to " << new_maxkmercount << " because of high coverage for genome size " << genome_size << endl;
                    m_min_count = new_min_count;                
                    m_low_count = m_min_count;
                    m_maxkmercount = new_maxkmercount;
                    sorted_kmers.RemoveLowCountKmers(m_min_count);
                }
            }
        }            
        if(kmer_counter.Kmers().Size() == 0)
            return 0;
                
        double average_count = kmer_counter.AverageCount();
        kmer_counter.GetBranches();

        map<int,size_t> hist;
        for(size_t index = 0; index < sorted_kmers.Size(); ++index) {
            ++hist[sorted_kmers.GetCount(index)];                  // count clipped to integer automatically
        }
        TBins bins(hist.begin(), hist.end());
        m_graphs[kmer_len] = new CDBGraph(move(sorted_kmers), move(bins), is_stranded);

        return average_count;
    }

    template<> template<> // one for graph, the other for args
    double CDBGAssembler<CDBHashGraph>::GetGraph(int kmer_len, const list<array<CReadHolder,2>>& reads, bool is_stranded, double total_seq, int estimated_kmer_num, bool skip_bloom_filter) {
        int64_t M = 1000000;
        CKmerHashCounter  kmer_counter(reads, kmer_len, m_min_count, M*estimated_kmer_num, is_stranded, m_ncores, skip_bloom_filter);
        if(kmer_counter.KmerNum() == 0)
            return 0;

        if(total_seq > 0) {
            TBins bins = kmer_counter.Kmers().GetBins();
            int genome_size = CalculateGenomeSize(bins);
            if(genome_size > 0) {
                int new_min_count = min(255.,total_seq/genome_size/50+0.5);
                if(new_min_count > m_min_count) {
                    int new_maxkmercount = max(10, int(total_seq/genome_size/10+0.5));
                    cerr << "WARNING: --min_count changed from " << m_min_count << " to " << new_min_count << " because of high coverage for genome size " << genome_size << endl;
                    cerr << "WARNING: --max_kmer_count changed from " << m_maxkmercount << " to " << new_maxkmercount << " because of high coverage for genome size " << genome_size << endl;
                    m_min_count = new_min_count;                
                    m_low_count = m_min_count;
                    m_maxkmercount = new_maxkmercount;
                    kmer_counter.RemoveLowCountKmers(m_min_count);
                }
            }
        }            
        if(kmer_counter.KmerNum() == 0)
            return 0;
                
        kmer_counter.GetBranches();
        m_graphs[kmer_len] = new CDBHashGraph(move(kmer_counter.Kmers()), is_stranded);

        return m_graphs[kmer_len]->AverageCount();
    }


} // namespace
#endif /* _DBGAssembler_ */

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

#ifndef _GuidedAssembler_
#define _GuidedAssembler_

#include <smmintrin.h>

using namespace std;
namespace DeBruijn {

class CGuidedAssembler {
public:
    CGuidedAssembler(int kmer_len, int secondary_kmer_len, bool extend_ends, bool protect_ends, bool keep_subgraphs, int min_count, double fraction, int word_size, int gap_open, int gap_extend, 
                     int drop_off, int kmer_complexity, double max_fork_density, int buf_length, int ncores, list<array<CReadHolder,2>>& raw_reads, 
                     int estimated_kmer_num, bool skip_bloom_filter, int not_aligned_len, int not_aligned_count, int aligned_count, int maxp, double target_coverage, int min_hit_len, 
                     bool no_reads, bool no_pairs,  int secondary_kmer_threshold, bool remove_homopolymer_indels, int homopolymer_len, double homopolymer_ratio) :
        m_kmer_len(kmer_len), m_secondary_kmer_len(secondary_kmer_len),
        m_extend_ends(extend_ends), m_protect_ends(protect_ends), m_keep_subgraphs(keep_subgraphs), m_min_count(min_count), m_fraction(fraction), m_word_size(word_size), 
        m_gap_open(gap_open), m_gap_extend(gap_extend), m_drop_off(drop_off), m_kmer_complexity(kmer_complexity), m_buf_length(buf_length), m_ncores(ncores), 
        m_not_aligned_len(not_aligned_len), m_not_aligned_count(not_aligned_count), m_aligned_count(aligned_count), m_maxp(maxp), m_target_coverage(target_coverage), m_min_hit_len(min_hit_len), 
        m_graph_uniformity(0), m_no_reads(no_reads), m_no_pairs(no_pairs), m_raw_reads(raw_reads),
        m_secondary_kmer_threshold(secondary_kmer_threshold), m_remove_homopolymer_indels(remove_homopolymer_indels), m_homopolymer_len(homopolymer_len), m_homopolymer_ratio(homopolymer_ratio) {

        if(m_kmer_complexity == 0)
            m_kmer_complexity = numeric_limits<int>::max();
        m_min_bases_per_fork = 0;
        if(max_fork_density > 0)
            m_min_bases_per_fork = 1./max_fork_density;

        m_graphp.reset(CreateGraph(estimated_kmer_num, skip_bloom_filter, m_kmer_len, m_min_count));
        if(m_secondary_kmer_len != m_kmer_len)
            m_secondary_graphp.reset(CreateGraph(estimated_kmer_num, skip_bloom_filter, m_secondary_kmer_len, m_min_count));
        else
            m_secondary_graphp = m_graphp;
    }
    virtual ~CGuidedAssembler() = default;

    void PrintRslt(ofstream& gfa_out, ofstream& all_variants, int maxv, ofstream& selected_variants, bool use_ambiguous) {

        /*
        for(auto& gfa_graph : m_gfa_collection) {
            gfa_graph.CutToChunks();
            gfa_graph.GenerateKmersAndScores(*m_secondary_graphp);
        }                             
        */        

        if(use_ambiguous) {
            SVarNum variants = 0;
            for(auto& gfa_graph : m_gfa_collection) {
                gfa_graph.SnpsToAmbig();
                variants += gfa_graph.NumberOfVariants();
            }
            cerr << "Variants after SNP collaps: " << variants.ToString() << endl;
        }

        for(auto& gfa_graph : m_gfa_collection)
            gfa_graph.PrintGFA(gfa_out);                

        gfa_out.close();
        if(!gfa_out) {
            cerr << "Can't write to gfa output" << endl;
            exit(1);
        }

        if(all_variants.is_open()) {
            set<string> seqs;
            for(auto& gfa_graph : m_gfa_collection) {
                if(!m_keep_subgraphs)
                    seqs.clear();
                gfa_graph.PrintAllVariants(all_variants, maxv, seqs);
            }
            all_variants.close();
            if(!all_variants) {
                cerr << "Can't write to all variants output" << endl;
                exit(1);
            }
        }

        if(selected_variants.is_open()) {
            for(auto& gfa_graph : m_gfa_collection)
                gfa_graph.PrintSelectedVariants(selected_variants);
            selected_variants.close();
            if(!selected_variants) {
                cerr << "Can't write to selected variants output" << endl;
                exit(1);
            }
        }
    }

    DBGraph& Graph() { return *m_graphp; }

protected:
    virtual void ReadTargets(ifstream& targets_in) = 0;
    virtual void AssemblerJob(TGFACollection& rslts) = 0;

    void AssembleGraphs() { // assemble
        vector<TGFACollection> thread_rslts(m_ncores);
            
        list<function<void()>> jobs;
        for(int thr = 0; thr < m_ncores; ++thr) {
            jobs.push_back(bind(&CGuidedAssembler::AssemblerJob, this, ref(thread_rslts[thr])));
        }
        RunThreads(m_ncores, jobs);

        for(auto& trslt : thread_rslts)
            m_gfa_collection.splice(m_gfa_collection.end(), trslt);
    }

    void Assemble(ifstream& targets_in, bool clip_to_codons) {
        int jump = 50;  //not really used
        int low_count = m_min_count;
        m_graphdiggerp.reset(new GraphDigger(*m_graphp, m_fraction, jump, low_count));                
        m_secondary_graphdiggerp.reset(new GraphDigger(*m_secondary_graphp, m_fraction, jump, low_count));                
        
        //read fasta
        ReadTargets(targets_in);

        CStopWatch timer;

        timer.Restart();
        CreateGraphHash();
        cerr << "Graph hash in " << timer.Elapsed();

        timer.Restart();
        AssembleGraphs();
        cerr << "Assembling in " << timer.Elapsed();

        timer.Restart();
        ClearGraphHash();
        cerr << "Graph hash clear in " << timer.Elapsed();
                                             
        if(m_keep_subgraphs) {
            SortCollection(m_gfa_collection);
            EnumerateCollection(m_gfa_collection);
        } else {
            timer.Restart();
            RemoveRedundantGraphs(m_gfa_collection);
            cerr << "Remove redundant in " << timer.Elapsed();         
        }

        SVarNum total_variants;
        for(auto& gfa_graph : m_gfa_collection) {
            SVarNum  vars = gfa_graph.NumberOfVariants();
            total_variants += vars;
            auto& acc = gfa_graph.Target();
            int group = gfa_graph.front().m_group;
            cerr << acc << ":" << group << " variants before " << vars.ToString() << endl; 
        }
        cerr << "Variants: " << total_variants.ToString() << endl;
              
        if((m_no_reads && m_no_pairs) || m_gfa_collection.empty())
            return; 

        GraphCleaner<TGFACollection> graph_cleaner(m_gfa_collection, &m_targets, *m_graphp, m_fraction, 0.51, m_not_aligned_len, m_not_aligned_count, m_aligned_count, m_maxp, 
                                   m_no_reads, m_no_pairs, m_raw_reads, m_ncores);
        
        {
            TGFACollection rslts;
            for(auto& gfa_graph : m_gfa_collection) {
                gfa_graph.RemoveHair(*m_secondary_graphdiggerp, m_fraction);
                gfa_graph.CalculateChainLength();
                int min_len = get<1>(m_targets[gfa_graph.Target()]);
                gfa_graph.TrimGroups(m_graph_uniformity, min_len);   
                gfa_graph.MergeRedundantLinks();
                if(clip_to_codons)
                    gfa_graph.ClipToCodons();
                gfa_graph.AssignGroupNumber();
                auto splitted = gfa_graph.SplitGroups();
                for(auto& graph : splitted) {
                    if(m_min_count == 1) {
                        graph.RemoveSinglReadSnps(*m_secondary_graphp);
                    }
                    if(graph.CollapsFreeEnds() > 0)
                        graph.MergeRedundantLinks();
                    if(m_extend_ends) {
                        graph.GenerateKmers(*m_graphp);                // kmers for calculating extension
                        graph.ExtendToFirstFork(*m_graphdiggerp);  
                    } 
                    if(m_remove_homopolymer_indels)
                        graph.RemoveHomopolymerIndels(*m_secondary_graphp, m_homopolymer_ratio, m_homopolymer_len);
                    graph.GenerateKmersAndScores(*m_secondary_graphp); // kmers for scores
                    graph.ScoreGraph(get<2>(m_targets[graph.Target()]), m_word_size);
                }
                rslts.splice(rslts.end(), splitted);
            }
            swap(m_gfa_collection, rslts);
        }
                       
        if(m_keep_subgraphs) {
            SortCollection(m_gfa_collection);
            EnumerateCollection(m_gfa_collection);
        } else {
            timer.Restart();
            RemoveRedundantGraphs(m_gfa_collection);
            cerr << "Remove redundant in " << timer.Elapsed(); 
        }

        total_variants = 0;
        for(auto& gfa_graph : m_gfa_collection) {
            SVarNum  vars = gfa_graph.NumberOfVariants();
            total_variants += vars;
            auto& acc = gfa_graph.Target();
            int group = gfa_graph.front().m_group;
            cerr << acc << ":" << group << " variants after " << vars.ToString() << endl; 
        }
        cerr << "Variants: " << total_variants.ToString() << endl;
    }

    CDBHashGraph* CreateGraph(int estimated_kmer_num, bool skip_bloom_filter, int kmer_len, int min_count) {
        int64_t M = 1000000;
        CKmerHashCounter  kmer_counter(m_raw_reads, kmer_len, min_count, M*estimated_kmer_num, true, m_ncores, skip_bloom_filter);
        if(kmer_counter.KmerNum() == 0)
            throw runtime_error("Not enough quality reads which are at least "+to_string(kmer_len)+"bp long");
        kmer_counter.GetBranches();
        return new CDBHashGraph(move(kmer_counter.Kmers()), true);
    }

    typedef vector<vector<pair<Node,int>>> TGraphHash;   // pair is node,count

    void CreateGraphHash() { // hash for kmers
        m_kmer_hash.resize(1ULL << 2*m_word_size);
        vector<SAtomic<uint8_t>> centinel(m_kmer_hash.size());
        {
            vector<typename DBGraph::Iterator> chunks = m_graphp->Chunks(m_ncores);
            list<function<void()>> jobs;
            for(int thr = 0; thr < (int)chunks.size()-1; ++thr)
                jobs.push_back(bind(&CGuidedAssembler::GraphHashJob, this, chunks[thr], chunks[thr+1], ref(centinel)));        
            RunThreads(m_ncores, jobs);
        }
        {
            list<function<void()>> jobs;
            for(int thr = 0; thr < m_ncores; ++thr)
                jobs.push_back(bind(&CGuidedAssembler::GraphHashSortJob, this, ref(centinel)));
            RunThreads(m_ncores, jobs);
        }
    }

    void GraphHashJob(typename DBGraph::Iterator from, typename DBGraph::Iterator to, vector<SAtomic<uint8_t>>& centinel) {
        uint64_t mask = (1ULL << 2*m_word_size)-1;
        for(auto it = from; it != to; ++it) { 
            int abundance = m_graphp->Abundance(it);
            if(abundance <= 1) // for min_count 1
                continue;
            TKmer kmer = m_graphp->GetNodeKmer(it);
            uint64_t word = *kmer.getPointer()&mask;             // last symbols of forward kmer
            while(!centinel[word].Set(1, 0));
            m_kmer_hash[word].emplace_back(it, abundance);                                
            centinel[word] = 0;
            TKmer rkmer = revcomp(kmer, m_kmer_len);
            word = *rkmer.getPointer()&mask;                     // last symbols of reversed kmer
            while(!centinel[word].Set(1, 0));
            m_kmer_hash[word].emplace_back(m_graphp->ReverseComplement(it), abundance);
            centinel[word] = 0;
        }
    }

    void GraphHashSortJob(vector<SAtomic<uint8_t>>& centinel) {
        for(size_t i = 0; i < centinel.size(); ++i) {
            if(!centinel[i].Set(1,0))
                continue;
            sort(m_kmer_hash[i].begin(), m_kmer_hash[i].end(), [](const pair<Node,int>& a, const pair<Node, int>& b) { return a.second > b.second; });
        }
    }

    void ClearGraphHash() {
        int step = m_kmer_hash.size()/m_ncores+1;            
        list<function<void()>> jobs;
        for(int thr = 0; thr < m_ncores; ++thr)
            jobs.push_back(bind(&CGuidedAssembler::GraphHashClearJob, this, step*thr, step));
        RunThreads(m_ncores, jobs);

        m_kmer_hash.clear();
    }

    void GraphHashClearJob(size_t from, size_t step) {
        for(size_t i = from; i < min(from+step, m_kmer_hash.size()); ++i)
            m_kmer_hash[i].clear();
    }

    /*
    int BlclLim(int block) {
        int num = pow(4, block);
        vector<double> transision(num+1);
        for(int n = 0; n <= num; ++n)
            transision[n] = double(n)/num;
        vector<double> previous(num+1);
        previous[0] = 1;
        for(int m = 0; m < m_kmer_len/block; ++m) {
            vector<double> next(num+1);
            for(int n = 0; n <= num; ++n) {
                double p = transision[n];
                next[n] += previous[n]*p;
                if(n < num)
                    next[n+1] += previous[n]*(1-p);
            }
            swap(previous, next);
        }

        int lim = 0;
        double eps = 0;
        for( ; lim < num && eps+previous[lim+1] <= m_lc_level; ++lim, eps += previous[lim]);

        return lim;
    }

    bool SimplicityCheck(const TKmer& kmer, int block, int lim) {
        uint64_t present = 0;                        // each bit represents different symbol combination
        uint64_t mask = (1ULL << 2*block) - 1;
        const uint64_t* data = kmer.getPointer();
        for(int s = 0; s <= m_kmer_len-block; s += block) {
            int p = s/32;
            int shift = s%32;
            int remain = 32-shift;
            uint64_t w = (data[p] >> 2*shift);
            if(remain < block)                        // need symbols from next word
                w |= (data[p+1] << 2*(block-remain));
            w &= mask;
            present |= (1ULL << w);
            if(_mm_popcnt_u64(present) > lim)
                return false;
        }
        return true;
    }
    bool SimplKmer(const TKmer& kmer) {
        return SimplicityCheck(kmer, 2, m_blc2) || SimplicityCheck(kmer, 3, m_blc3);
    }
    int BCounter(const TKmer& kmer, int block) {
        uint64_t present = 0;                        // each bit represents different symbol combination    
        uint64_t mask = (1ULL << 2*block) - 1;
        const uint64_t* data = kmer.getPointer();
        for(int s = 0; s <= m_kmer_len-block; s += block) {
            int p = s/32;
            int shift = s%32;
            int remain = 32-shift;
            uint64_t w = (data[p] >> 2*shift);
            if(remain < block)                        // need symbols from next word    
                w |= (data[p+1] << 2*(block-remain));
            w &= mask;
            present |= (1ULL << w);
        }
        return _mm_popcnt_u64(present);
    }
    */

    bool CheckSeed(Node seed, int left_target_len, int right_target_len) {
        int margin = numeric_limits<int>::min();
        if(m_protect_ends)
            margin = max(100, m_kmer_len);
        for(int dir = 0; dir < 2; ++dir) {
            set<Node> kmers;
            kmers.insert(seed);
            for(int step = 0; step < m_kmer_len; ++step) {
                set<Node> new_kmers;
                for(auto& node : kmers) {
                    int remaining_len = right_target_len-step;
                    int back_step_len = left_target_len+step+1; // 1 because we will make a 1bp step
                    bool check_forward = margin < remaining_len;
                    bool check_backward = margin < back_step_len;
                    vector<Successor> neighbors = m_graphdiggerp->GetReversibleNodeSuccessorsF(node, nullptr, check_forward, check_backward);
                    for(auto& suc : neighbors)
                        new_kmers.insert(suc.m_node);
                }
                if(new_kmers.empty())
                    return false;
                swap(kmers, new_kmers);
            }
            seed = seed.ReverseComplement();
            swap(left_target_len, right_target_len);
        }
        
        return true;
    }

    TGraphHash m_kmer_hash;

    TGFACollection m_gfa_collection;
    TTargets m_targets;                         // [accession], seq, min_len, words, centinel
    unique_ptr<GraphDigger> m_graphdiggerp;
    unique_ptr<GraphDigger> m_secondary_graphdiggerp;
    shared_ptr<DBGraph> m_graphp;
    shared_ptr<DBGraph> m_secondary_graphp;
    int m_kmer_len;
    int m_secondary_kmer_len;
    bool m_extend_ends;
    bool m_protect_ends;
    bool m_keep_subgraphs;
    int m_min_count;
    double m_fraction;

    int m_word_size;
    int m_gap_open;
    int m_gap_extend;
    int m_drop_off;

    unsigned m_kmer_complexity;
    double m_min_bases_per_fork;
    int m_buf_length;

    int m_ncores;
    int m_not_aligned_len;
    int m_not_aligned_count;
    int m_aligned_count;
    int m_maxp;
    double m_target_coverage;
    int m_min_hit_len;
    double m_graph_uniformity;
    bool m_no_reads;
    bool m_no_pairs;

    list<array<CReadHolder,2>>& m_raw_reads;

    int m_secondary_kmer_threshold;
    bool m_remove_homopolymer_indels;
    int m_homopolymer_len;
    double m_homopolymer_ratio;

    mutex m_out_mutex;
};

       
class CGuidedAssemblerNA : public CGuidedAssembler {
public:

    CGuidedAssemblerNA(int kmer_len, int secondary_kmer_len, bool extend_ends, bool protect_ends, bool keep_subgraphs, int min_count, double fraction, int seed_prec, int word_size, int match, int mismatch, int gap_open, int gap_extend, 
                       int drop_off, int kmer_complexity, double max_fork_density, int buf_length, int ncores, list<array<CReadHolder,2>>& raw_reads, ifstream& targets_in, 
                       int estimated_kmer_num, bool skip_bloom_filter, int not_aligned_len, int not_aligned_count, int aligned_count, int maxp, double target_coverage, int min_hit_len, 
                       bool no_reads, bool no_pairs, int secondary_kmer_threshold, bool remove_homopolymer_indels, int homopolymer_len, double homopolymer_ratio) : 
        CGuidedAssembler(kmer_len, secondary_kmer_len, extend_ends, protect_ends, keep_subgraphs, min_count, fraction, word_size, gap_open, gap_extend, 
                         drop_off, kmer_complexity, max_fork_density, buf_length, ncores, raw_reads, 
                         estimated_kmer_num, skip_bloom_filter, not_aligned_len, not_aligned_count, aligned_count, maxp, target_coverage, min_hit_len, 
                         no_reads, no_pairs, secondary_kmer_threshold, remove_homopolymer_indels, homopolymer_len, homopolymer_ratio),  
                         m_seed_prec(seed_prec), m_match(match), m_mismatch(mismatch), m_delta(match, mismatch) {

        Assemble(targets_in, false);
    }

private:
    void ReadTargets(ifstream& targets_in) {
        char c;
        if(!(targets_in >> c) || c != '>')
            throw runtime_error("Invalid fasta file format for targets");
        string record;
        while(getline(targets_in, record, '>')) {
            while(!targets_in.eof() && record.back() != '\n') {
                string part;
                getline(targets_in, part, '>');
                record += '>'+part;
            }
            size_t first_ret = min(record.size(),record.find('\n'));
            if(first_ret == string::npos)
                throw runtime_error("Invalid fasta file format for targets");
            string acc = record.substr(0, first_ret);
            acc = acc.substr(0, acc.find_first_of(" \t"));
            string target = record.substr(first_ret+1);
            target.erase(remove(target.begin(),target.end(),'\n'),target.end());
            for(char& c : target) c = toupper(c);
            if(target.find_first_not_of("ACGTYRWSKMDVHBXN") != string::npos)
                throw runtime_error("Invalid sequence in fasta file for targets");
            get<0>(m_targets[acc]) = target;
            get<1>(m_targets[acc]) = min(target.size()*m_target_coverage, (double)m_min_hit_len); 
        }
        if(targets_in.bad())
            throw runtime_error("Error in reading targets");
    }

    void AssemblerJob(TGFACollection& rslts) {

        for(auto& item : m_targets) {
            if(!get<3>(item.second).Set(1,0))
                continue;

            string& target = get<0>(item.second);
            const string& acc = item.first;
            double anchor_frac = 0.25;
            int min_len = get<1>(item.second);
            int target_len = target.size();
        
            if((int)target.size() < m_kmer_len) {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Skipped short target: " << acc << " " << target.size() << endl;
                continue;
            }

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Started assembling: " << acc << endl;
            }

            string tcopy;
            for(char c : target)
                tcopy.push_back(toupper(c));

            for(int p = 0; p <= (int)tcopy.size()-m_word_size; ++p) {
                string seed = tcopy.substr(p, m_word_size);
                if(seed.find_first_not_of("ACGT") != string::npos)
                    continue;
                uint32_t word = 0;
                for(char c : seed) {
                    word = word << 2;
                    word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                }
                get<2>(item.second).insert(word);
            }

            size_t collapsedL = 0;
            size_t collapsedR = 0;
   
            CStopWatch timer;
            timer.Restart();

            double seed_prec = min((double)m_kmer_len-1, (double)m_mismatch/(m_match+m_mismatch)*m_kmer_len+m_seed_prec);

            vector<set<Node>> kmers_for_pos(target_len-(m_kmer_len-1));
            vector<int> max_counts_for_pos(target_len-(m_kmer_len-1));

            for(int p = 0; p < target_len-m_word_size; ++p) {
                string seed = tcopy.substr(p, m_word_size); // last symbols
                if(seed.find_first_not_of("ACGT") != string::npos)
                    continue;

                uint32_t word = 0;
                for(char c : seed) {
                    word = word << 2;
                    word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                }

                array<int,4> pos = {p+m_word_size-m_kmer_len, p+m_word_size-m_kmer_len, p, p};   // before word direct, before word reversed, after word direct, after word reversed
                array<unique_ptr<TKmer>,4> tkmerp;
                if(pos[0] >= 0) {
                    string tkmer = tcopy.substr(pos[0], m_kmer_len);        // kmer before word
                    if(tkmer.find_first_not_of("ACGT") == string::npos) {
                        tkmerp[0].reset(new TKmer(tkmer)); 
                        ReverseComplementSeq(tkmer.begin(), tkmer.end());
                        tkmerp[1].reset(new TKmer(tkmer)); 
                    }
                }
                if(pos[2] <= target_len-m_kmer_len) {
                    string tkmer = tcopy.substr(p, m_kmer_len);             // kmer after word
                    if(tkmer.find_first_not_of("ACGT") == string::npos) {
                        tkmerp[2].reset(new TKmer(tkmer));
                        ReverseComplementSeq(tkmer.begin(), tkmer.end());
                        tkmerp[3].reset(new TKmer(tkmer));
                    }
                }
                if(!tkmerp[0] && !tkmerp[2])
                    continue;
                 
                auto CountMatches = [](uint64_t a, uint64_t b) {
                    uint64_t w = ~(a^b);       // each match will produce 11                                           
                    w &= (w << 1);             // upper bit is 1 only if both bits are 1 (only for matches)
                    w &= 0xAAAAAAAAAAAAAAAAUL; // each match has 10; all mismatches 00
                    return _mm_popcnt_u64(w);  // count number set of bits == matches 
                };

                for(int wdir = 0; wdir < 2; ++wdir) {
                    for(auto& node_count : m_kmer_hash[word]) {
                        Node node = node_count.first;
                        int kdir = node.isMinus();
                        if(wdir > 0)
                            node = node.ReverseComplement();

                        int templ;
                        if(wdir == 0 && kdir == 0)
                            templ = 0;
                        else if(wdir == 0 && kdir > 0)
                            templ = 1;
                        else if(wdir > 0 && kdir == 0)
                            templ = 3;
                        else
                            templ = 2;                        

                        if(tkmerp[templ]) { 
                            int leftp = pos[templ];
                            int count = node_count.second;
                            const uint64_t* targetp = tkmerp[templ]->getPointer();
                            const uint64_t* graphp = m_graphp->getPointer(node);
                            int prec = tkmerp[templ]->getSize()/64;// number of 8-byte blocks   
                            int matches = -(prec*32-m_kmer_len);   // empty positions will match    
                            for(int i = 0; i < prec; ++i)
                                matches += CountMatches(*(targetp+i), *(graphp+i));
                            if(matches > seed_prec) {
                                kmers_for_pos[leftp].insert(node);                                    
                                max_counts_for_pos[leftp] = max(max_counts_for_pos[leftp], count);
                                if(kmers_for_pos[leftp].size() >= m_kmer_complexity)
                                    break;
                            }
                        }
                    }
                    
                    if(wdir == 0) {
                        //reverse
                        word = ((word & 0x33333333) << 2)  | ((word >> 2)  & 0x33333333); // swap adjacent pairs
                        word = ((word & 0x0F0F0F0F) << 4)  | ((word >> 4)  & 0x0F0F0F0F); // swap nibbles
                        word = ((word & 0x00FF00FF) << 8)  | ((word >> 8)  & 0x00FF00FF); // swap bytes
                        word = ((word & 0x0000FFFF) << 16) | ((word >> 16) & 0x0000FFFF); // swap 16 bit chunks
                        //complement
                        word ^= 0xAAAAAAAA;
                        //shift
                        word >>= 2*(16-m_word_size);
                    }
                }
            }

            /*            
            {
                lock_guard<mutex> guard(m_out_mutex);
                for(int l = 0; l < kmers_for_pos.size(); ++l)
                    cerr << "Word hits: " << acc << " " << l << " " << kmers_for_pos[l].size() << endl;
            }
            */            

            {
                int ksize = kmers_for_pos.size();
                vector<unsigned> kcount(ksize);
                for(int p = 0; p < ksize; ++p)
                    kcount[p] = kmers_for_pos[p].size();
                for(int p = 0; p < ksize; ++p) {
                    if(kcount[p] >= m_kmer_complexity) {
                        for(int l = p; l < p+m_kmer_len; ++l)
                            target[l] = 'N';
                        for(int l = max(0,p-m_kmer_len+1); l < min(p+m_kmer_len,ksize); ++l)
                            kmers_for_pos[l].clear();                        
                    }
                }
            }
            int simpl_kmers = 0;
            for(char c : target) {
                if(c == 'N')
                    ++simpl_kmers;
            }

            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                auto& kfp = kmers_for_pos[p];
                int mcount = max_counts_for_pos[p];
                for(auto it_loop = kfp.begin(); it_loop != kfp.end(); ) {
                    auto it = it_loop++;
                    auto abundance = m_graphp->Abundance(*it);
                    if(abundance == 1 || abundance < m_fraction*mcount) 
                        kfp.erase(it);
                }
            }

            unordered_map<Node, list<int>, Node::Hash> positions;
            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                for(Node node : kmers_for_pos[p]) 
                    positions[node].push_back(p);
            }
            set<int> periodic_points;
            for(auto pos : positions) {
                if(pos.second.size() > 1) {
                    for(int p : pos.second)
                        periodic_points.insert(p);
                }
            }
                                
            unordered_map<Node, int, Node::Hash> target_kmers; // [node] position on target
            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                if(periodic_points.count(p))
                    continue;
                for(Node node : kmers_for_pos[p])
                    target_kmers.emplace(node, p);
            }

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Periodic points for : " << acc << " " << periodic_points.size() << " " << target_len << "/" << simpl_kmers << "\n";
                cerr << "Seed kmers for " << acc << " " << target_kmers.size() << " in " << timer.Elapsed();
                timer.Restart();
            }                                             

            CGuidedGraph gugraph(m_kmer_len, m_secondary_kmer_len);

            while(!target_kmers.empty()) {
                unordered_map<Node, int, Node::Hash>::iterator first_kmerp;
                typedef unordered_map<Node, int, Node::Hash>::value_type elem_t;
                first_kmerp = min_element(target_kmers.begin(), target_kmers.end(), 
                                          [this](const elem_t& a, const elem_t& b) 
                                          { 
                                              if(a.second != b.second)
                                                  return a.second < b.second; 
                                              else if(m_graphp->Abundance(a.first) != m_graphp->Abundance(b.first))
                                                  return m_graphp->Abundance(a.first) > m_graphp->Abundance(b.first);
                                              else
                                                  return m_graphp->GetNodeSeq(a.first) < m_graphp->GetNodeSeq(b.first);
                                          });
                    
                Node initial_node = first_kmerp->first;
                string kmer_seq = m_graphp->GetNodeSeq(initial_node);
                int first_matching_kmer = first_kmerp->second; // position on target
                target_kmers.erase(first_kmerp);
                auto inode = m_graphp->GetNodeKmer(initial_node);  

                if(!CheckSeed(initial_node, first_matching_kmer, target_len-first_matching_kmer-m_kmer_len))
                    continue;

                int left_not_aligned = 0;
                int right_not_aligned = 0;
                int left_penalty = 0;
                int right_penalty = 0;
                {
                    int lminscore = numeric_limits<int>::max();
                    int lminpos = -1;
                    int score = 0;
                    for(int p = 0; p < m_kmer_len; ++p) {
                        if(kmer_seq[p] == target[first_matching_kmer+p])
                            score += m_match;
                        else
                            score -= m_mismatch;
                        if(score <= lminscore) {
                            lminscore = score;
                            lminpos = p;
                        }
                    }
                    int rminscore = numeric_limits<int>::max();
                    int rminpos = -1;
                    score = 0;
                    for(int p = 0; p < m_kmer_len; ++p) {
                        if(kmer_seq[m_kmer_len-1-p] == target[first_matching_kmer+m_kmer_len-1-p])
                            score += m_match;
                        else
                            score -= m_mismatch;
                        if(score <= rminscore) {
                            rminscore = score;
                            rminpos = p;
                        }
                    }
                    if(lminscore <= 0) {
                        left_not_aligned = lminpos+1;
                        left_penalty = -lminscore;
                    }
                    if(rminscore <= 0) {
                        right_not_aligned = rminpos+1;
                        right_penalty = -rminscore;
                    }
                }
                if(left_not_aligned+right_not_aligned > 0.5*m_kmer_len)
                    continue;

                gugraph.StartNewAssembly(kmer_seq, left_not_aligned, right_not_aligned);
                size_t assembled_len = 0;

                if(first_matching_kmer > 0) {                  // needs left extension
                    string left_target = target.substr(0, first_matching_kmer);
                    ReverseComplementSeq(left_target.begin(), left_target.end());
                    CGuidedPathNA extender(m_graphp->ReverseComplement(initial_node), left_penalty, left_not_aligned, m_protect_ends, left_target, target_len-first_matching_kmer-m_kmer_len, 
                                           *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, m_gap_open, m_gap_extend, m_drop_off, anchor_frac, m_secondary_kmer_threshold);
                    map<typename CGuidedGraph::TAnchor, int> lanchors;
                    map<Node, int> lnotaligned;
                    int total = 0;
                    deque<int> forks(m_buf_length+1);
                    unordered_set<Node, Node::Hash> dead_ends;
                    while(extender.ProcessNextEdge()) {
                        if(extender.PathEnd() || (extender.NotAligned() >= m_kmer_len && dead_ends.count(extender.LastStepNode()))) {
                            if(extender.NotAligned() >= m_kmer_len) {
                                int nkmer = extender.NotAligned()-m_kmer_len+1;
                                for(int p = extender.AssembledSeqLength()-nkmer; p < extender.AssembledSeqLength(); ++p) {
                                    auto& node = extender.AssembledSeq()[p].m_node;
                                    if(node.isValid())
                                        dead_ends.insert(node);
                                }
                            }
                            SPathChunk chunk = extender.GetLastSegment();
                            total += chunk.m_seq.size();
                            gugraph.AddLeftSegment(chunk, lanchors, lnotaligned, gugraph.End()); 
                            int check_len = extender.AssembledSeqLength();
                            check_len = min(check_len, check_len-extender.NotAligned()+m_kmer_len-1);
                            for(int p = extender.StartingShift(); p < check_len; ++p)
                                target_kmers.erase(m_graphp->ReverseComplement(extender.AssembledSeq()[p].m_node));
                            lanchors.clear();
                            lnotaligned.clear();
                            if(extender.PathEnd())
                                extender.DeleteNotAlignedForks(extender.AssembledSeqLength()-extender.NotAligned()+m_kmer_len);
                            else
                                extender.DeleteLastBranch();
                            gugraph.RewindLeftBranch(extender.StartingShift(), anchor_frac);
                        } else if(extender.LastStepNode().isValid()) {
                            typename CGuidedGraph::TAnchor anchor(extender.LastStepNode(), first_matching_kmer-1-extender.GetMaxPos());  // left end of alignment on target
                            CGuidedGraph::TSegmentP hook(gugraph.End());
                            int chunk_len = extender.AssembledSeqLength()-extender.StartingShift();

                            if(extender.SolidKmer()) {
                                hook = gugraph.KnownLeftAnchor(anchor);
                                lanchors[anchor] = chunk_len-1;
                            } else if(chunk_len > m_kmer_len) {
                                lnotaligned[extender.LastStepNode()] = chunk_len-1;
                                hook = gugraph.KnownLeftNotAligned(extender.LastStepNode(), chunk_len, extender.StartingShift(), extender.AssembledSeq());
                                if(hook != gugraph.End())
                                    ++collapsedL;
                            }                                                                                                                                             

                            if(hook != gugraph.End()) {
                                SPathChunk chunk = extender.GetLastSegment();
                                total += chunk.m_seq.size();
                                gugraph.AddLeftSegment(chunk, lanchors, lnotaligned, hook);
                                for(int p = extender.StartingShift(); p < extender.AssembledSeqLength(); ++p)
                                    target_kmers.erase(m_graphp->ReverseComplement(extender.AssembledSeq()[p].m_node));
                                lanchors.clear();
                                lnotaligned.clear();
                                extender.DeleteLastBranch();
                                gugraph.RewindLeftBranch(extender.StartingShift(), anchor_frac);
                            }
                        }
                        
                        forks.pop_front();
                        forks.push_back(extender.ForkCount());
                        if(total > 10*target_len || (total > 0 && m_min_bases_per_fork*(forks.back()-forks.front()) > m_buf_length)) {
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted assemblyA " << acc << " kmer " << kmer_seq << " " << first_matching_kmer << endl;
                            break;                        
                        }                                                
                    }
                    assembled_len += total;
                }

                gugraph.CleanBranch();
                if(first_matching_kmer < target_len-m_kmer_len) {                        
                    string right_target =  target.substr(first_matching_kmer+m_kmer_len);  
                    CGuidedPathNA extender(initial_node, right_penalty, right_not_aligned, m_protect_ends, target.substr(first_matching_kmer+m_kmer_len), first_matching_kmer, 
                                           *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, m_gap_open, m_gap_extend, m_drop_off, anchor_frac, m_secondary_kmer_threshold);
                    map<typename CGuidedGraph::TAnchor, int> ranchors;
                    map<Node, int> rnotaligned;
                    int total = 0;
                    deque<int> forks(m_buf_length+1);
                    unordered_set<Node, Node::Hash> dead_ends;  
                    while(extender.ProcessNextEdge()) {
                        if(extender.PathEnd() || (extender.NotAligned() >= m_kmer_len && dead_ends.count(extender.LastStepNode()))) {
                            if(extender.NotAligned() >= m_kmer_len) {
                                int nkmer = extender.NotAligned()-m_kmer_len+1;
                                for(int p = extender.AssembledSeqLength()-nkmer; p < extender.AssembledSeqLength(); ++p) {
                                    auto& node = extender.AssembledSeq()[p].m_node;
                                    if(node.isValid())
                                        dead_ends.insert(node);
                                }
                            }
                            SPathChunk chunk = extender.GetLastSegment();
                            total += chunk.m_seq.size();
                            gugraph.AddRightSegment(chunk, ranchors, rnotaligned, gugraph.End()); 
                            int check_len = extender.AssembledSeqLength();
                            check_len = min(check_len, check_len-extender.NotAligned()+m_kmer_len-1);
                            for(int p = extender.StartingShift(); p < check_len; ++p)
                                target_kmers.erase(extender.AssembledSeq()[p].m_node);
                            ranchors.clear();
                            rnotaligned.clear();
                            if(extender.PathEnd())
                                extender.DeleteNotAlignedForks(extender.AssembledSeqLength()-extender.NotAligned()+m_kmer_len);
                            else
                                extender.DeleteLastBranch();
                            gugraph.RewindRightBranch(extender.StartingShift(), anchor_frac);
                        } else if(extender.LastStepNode().isValid()) {
                            CGuidedGraph::TSegmentP hook(gugraph.End());
                            int chunk_len = extender.AssembledSeqLength()-extender.StartingShift();

                            if(extender.SolidKmer()) {
                                typename CGuidedGraph::TAnchor anchor(extender.LastStepNode(), first_matching_kmer+m_kmer_len+extender.GetMaxPos());  // right end of alignment on target
                                hook = gugraph.KnownRightAnchor(anchor);
                                ranchors[anchor] = chunk_len-1;                                   
                            } else if(chunk_len > m_kmer_len) {
                                rnotaligned[extender.LastStepNode()] = chunk_len-1;
                                hook = gugraph.KnownRightNotAligned(extender.LastStepNode(), chunk_len, extender.StartingShift(), extender.AssembledSeq());
                                if(hook != gugraph.End())
                                    ++collapsedR;
                            }                                                                                                                                            

                            if(hook != gugraph.End()) {
                                SPathChunk chunk = extender.GetLastSegment();
                                total += chunk.m_seq.size();
                                gugraph.AddRightSegment(chunk, ranchors, rnotaligned, hook);
                                for(int p = extender.StartingShift(); p < extender.AssembledSeqLength(); ++p)
                                    target_kmers.erase(extender.AssembledSeq()[p].m_node);
                                ranchors.clear();
                                rnotaligned.clear();
                                extender.DeleteLastBranch();
                                gugraph.RewindRightBranch(extender.StartingShift(), anchor_frac);
                            }                                                                                                
                        } 
                                               
                        forks.pop_front();
                        forks.push_back(extender.ForkCount());
                        if(total > 10*target_len || (total > 0 && m_min_bases_per_fork*(forks.back()-forks.front()) > m_buf_length)) {
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted assemblyB " << acc << " kmer " << kmer_seq << " " << first_matching_kmer << endl;
                            break;                        
                        }                                               
                    }
                    assembled_len += total;
                }
                gugraph.RemoveNotAlignedSegments(anchor_frac);
            }

            GFAGraph gfa_graph(acc, m_kmer_len);
            gugraph.GetGFAGraph(gfa_graph);
            gfa_graph.CalculateChainLength();                        
            //            gfa_graph.RemoveShortChains(min_len);
            gfa_graph.AssignGroupNumber();
            gfa_graph.TrimGroups(m_graph_uniformity, min_len);  
            gfa_graph.MergeRedundantLinks();
            gfa_graph.AssignGroupNumber(); 
            size_t gsize = gfa_graph.Size();

            TGFACollection splitted = gfa_graph.SplitGroups();
            for(auto& graph : splitted) {
                if(m_min_count == 1) {
                    if(graph.RemoveHair(*m_secondary_graphdiggerp, m_fraction))
                        graph.MergeForks(); 
                    graph.RemoveSinglReadSnps(*m_secondary_graphp);
                }
                if(m_remove_homopolymer_indels)
                    graph.RemoveHomopolymerIndels(*m_secondary_graphp, m_homopolymer_ratio, m_homopolymer_len);
                if(m_no_reads && m_no_pairs) {
                    if(graph.CollapsFreeEnds() > 0)
                        graph.MergeRedundantLinks();
                    if(m_extend_ends) {
                        graph.GenerateKmers(*m_graphp);                
                        graph.ExtendToFirstFork(*m_graphdiggerp);                
                    }
                    graph.GenerateKmersAndScores(*m_secondary_graphp);
                } else {
                    graph.GenerateKmersAndScores(*m_graphp);
                }
                graph.ScoreGraph(get<2>(item.second), m_word_size);
            }
            rslts.splice(rslts.end(), splitted);

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Finished assembling: " << acc << " " << collapsedL << " " << collapsedR << " " << gsize << " in " << timer.Elapsed();
            }
        } 
    }
                           
    int m_seed_prec;
    int m_match;
    int m_mismatch;
    SMatrix m_delta;
};


class CGuidedAssemblerAA : public CGuidedAssembler {
public:

    CGuidedAssemblerAA(int kmer_len, int secondary_kmer_len, bool extend_ends, bool protect_ends, bool keep_subgraphs, int min_count, double fraction, int word_size, int gap_open, int gap_extend, 
                       int drop_off, int kmer_complexity, double max_fork_density, int buf_length, int ncores, list<array<CReadHolder,2>>& raw_reads, ifstream& targets_in, 
                       int estimated_kmer_num, bool skip_bloom_filter, int not_aligned_len, int not_aligned_count, int aligned_count, int maxp, double target_coverage, int min_hit_len, 
                       bool no_reads, bool no_pairs, int genetic_code, bool translate_na, int fs_open, bool allow_frameshifts, 
                       int secondary_kmer_threshold, bool remove_homopolymer_indels, int homopolymer_len, double homopolymer_ratio) : 
        CGuidedAssembler(kmer_len, secondary_kmer_len, extend_ends, protect_ends, keep_subgraphs, min_count, fraction, word_size, gap_open, gap_extend, 
                         drop_off, kmer_complexity, max_fork_density, buf_length, ncores, raw_reads, 
                         estimated_kmer_num, skip_bloom_filter, not_aligned_len, not_aligned_count, aligned_count, maxp, target_coverage, min_hit_len, 
                         no_reads, no_pairs, secondary_kmer_threshold, remove_homopolymer_indels, homopolymer_len, homopolymer_ratio), m_genetic_code(genetic_code), m_translate_na(translate_na), 
                         m_fs_open(fs_open), m_allow_frameshifts(allow_frameshifts) {

        Assemble(targets_in, !m_allow_frameshifts);
    }
    void TranslateToAA() {
        for(auto& gfa_graph : m_gfa_collection)
            gfa_graph.TranslateToAA(m_genetic_code, *m_graphp);

    }

private:
    void ReadTargets(ifstream& targets_in) {
        string accepted_symbols = m_translate_na ? "ACGTYRWSKMDVHBXN" : "UARNDCQEGHILKMFPSTWYVBJZX*";
        char c;
        if(!(targets_in >> c) || c != '>')
            throw runtime_error("Invalid fasta file format for targets");
        string record;
        while(getline(targets_in, record, '>')) {
            while(!targets_in.eof() && record.back() != '\n') {
                string part;
                getline(targets_in, part, '>');
                record += '>'+part;
            }
            size_t first_ret = min(record.size(),record.find('\n'));
            if(first_ret == string::npos)
                throw runtime_error("Invalid fasta file format for targets");
            string acc = record.substr(0, first_ret);
            acc = acc.substr(0, acc.find_first_of(" \t"));
            string target = record.substr(first_ret+1);
            target.erase(remove(target.begin(),target.end(),'\n'),target.end());
            for(char& c : target) c = toupper(c);
            if(target.find_first_not_of(accepted_symbols) != string::npos)
                throw runtime_error("Invalid sequence in fasta file for targets");
            if(m_translate_na) {
                target = m_genetic_code.Translate(target, true);
            } else {
                replace(target.begin(),target.end(),'U', 'X');
                replace(target.begin(),target.end(),'J', 'X');
            }
                    
            get<0>(m_targets[acc]) = target;
            get<1>(m_targets[acc]) = min(3*target.size()*m_target_coverage, (double)m_min_hit_len);
        }
        if(targets_in.bad())
            throw runtime_error("Error in reading targets");
    }

    void AssemblerJob(TGFACollection& rslts) {
        
        for(auto& item : m_targets) {
            if(!get<3>(item.second).Set(1,0))
                continue;

            string& target = get<0>(item.second);
            const string& acc = item.first;
            double anchor_frac = 0.25;
            int min_len = get<1>(item.second);
            int target_len = target.size();

            if(target_len*3 < m_kmer_len) {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Skipped short target: " << acc << " " << target_len << endl;
                continue;
            }

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Started assembling: " << acc << endl;
            }

            string tcopy;
            for(char c : target)
                tcopy.push_back(toupper(c));
            bool tstart = (tcopy[0] == 'M');

            for(int p = 0; p <= (int)tcopy.size()-m_word_size/3; ++p) {
                string aa_seed = tcopy.substr(p, m_word_size/3);
                if(aa_seed.find_first_not_of("ARNDCQEGHILKMFPSTWYV") != string::npos) // B and Z removed
                    continue;
                list<string> seeds;
                if(p == 0 && tstart)
                    seeds = m_genetic_code.Starts();
                else
                    seeds = m_genetic_code.Codons(aa_seed[0]);
                for(int i = 1; i < m_word_size/3; ++i) {
                    list<string> codons = m_genetic_code.Codons(aa_seed[i]);
                    for(auto iseed = seeds.begin(); iseed != seeds.end(); ++iseed) {
                        for(auto it = next(codons.begin()); it != codons.end(); ++it)
                            seeds.push_front(*iseed+*it);
                        *iseed += codons.front();
                    }
                }
                for(auto& seed : seeds) {
                    uint32_t word = 0;
                    for(char c : seed) {
                        word = word << 2;
                        word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                    }
                    get<2>(item.second).insert(word);
                }
            }
   
            CStopWatch timer;
            timer.Restart();

            vector<set<Node>> kmers_for_pos(target_len-(m_kmer_len/3-1));
            vector<int> max_counts_for_pos(target_len-(m_kmer_len/3-1));

            for(int p = 0; p <= target_len-m_word_size/3; ++p) {
                string aa_seed = tcopy.substr(p, m_word_size/3); // last symbols
                if(aa_seed.find_first_not_of("ARNDCQEGHILKMFPSTWYV") != string::npos) // B and Z removed
                    continue;

                string target_protl;
                int selfscorel = 0;
                if(p+(m_word_size-m_kmer_len)/3 >= 0) {
                    target_protl = tcopy.substr(p+(m_word_size-m_kmer_len)/3, m_kmer_len/3);
                    if(target_protl.find_first_not_of("ARNDCQEGHILKMFPSTWYVBZX*")  != string::npos) {
                        target_protl.clear();
                    } else {
                        for(int p = 0; p < m_kmer_len/3; ++p)
                            selfscorel += m_delta.matrix[(int)target_protl[p]][(int)target_protl[p]];
                    }
                }
                string target_protr;
                int selfscorer = 0;
                if(p <= target_len-m_kmer_len/3) {
                    target_protr = tcopy.substr(p, m_kmer_len/3);
                    if(target_protr.find_first_not_of("ARNDCQEGHILKMFPSTWYVBZX*")  != string::npos) {
                        target_protr.clear();
                    } else {
                        for(int p = 0; p < m_kmer_len/3; ++p)
                            selfscorer += m_delta.matrix[(int)target_protr[p]][(int)target_protr[p]];
                    }
                }
                if(target_protl.empty() && target_protr.empty())
                    continue;

                list<string> seeds;
                if(p == 0 && tstart)
                    seeds = m_genetic_code.Starts();
                else
                    seeds = m_genetic_code.Codons(aa_seed[0]);
                for(int i = 1; i < m_word_size/3; ++i) {
                    list<string> codons = m_genetic_code.Codons(aa_seed[i]);
                    for(auto iseed = seeds.begin(); iseed != seeds.end(); ++iseed) {
                        for(auto it = next(codons.begin()); it != codons.end(); ++it)
                            seeds.push_front(*iseed+*it);
                        *iseed += codons.front();
                    }
                }

                for(auto& seed : seeds) {
                    uint32_t word = 0;
                    for(char c : seed) {
                        word = word << 2;
                        word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                    }

                    double fr = 0.75;
                    for(int wdir = 0; wdir < 2; ++wdir) {
                        for(auto& node_count : m_kmer_hash[word]) {
                            Node node = node_count.first;
                            int count = node_count.second;
                            if(wdir == 0) {
                                //before word
                                if(target_protl.empty())
                                    break;
                                int leftp = p+(m_word_size-m_kmer_len)/3;
                                bool alts = (p == (m_kmer_len-m_word_size)/3 && tstart);
                                string kmer_prot = m_genetic_code.Translate(m_graphp->GetNodeKmer(node), m_kmer_len, alts);
                                if(kmer_prot.find('*') != string::npos)
                                    continue;
                                int score = 0;
                                for(int p = 0; p < m_kmer_len/3; ++p)
                                    score += m_delta.matrix[(int)kmer_prot[p]][(int)target_protl[p]];
                                if(score > fr*selfscorel) {
                                    kmers_for_pos[leftp].insert(node);
                                    max_counts_for_pos[leftp] = max(max_counts_for_pos[leftp], count);
                                    if(kmers_for_pos[leftp].size() >= m_kmer_complexity)
                                        break;
                                }
                            } else {
                                //after word
                                if(target_protr.empty())
                                    break;
                                int leftp = p;
                                node = node.ReverseComplement();
                                bool alts = (p == 0 && tstart);
                                string kmer_prot = m_genetic_code.Translate(m_graphp->GetNodeKmer(node), m_kmer_len, alts);
                                if(kmer_prot.find('*') != string::npos)
                                    continue;
                                int score = 0;
                                for(int p = 0; p < m_kmer_len/3; ++p)
                                    score += m_delta.matrix[(int)kmer_prot[p]][(int)target_protr[p]];
                                if(score > fr*selfscorer) {
                                    kmers_for_pos[leftp].insert(node);
                                    max_counts_for_pos[leftp] = max(max_counts_for_pos[leftp], count);
                                    if(kmers_for_pos[leftp].size() >= m_kmer_complexity)
                                        break;
                                }
                            }
                        }

                        if(wdir == 0) {
                            //reverse
                            word = ((word & 0x33333333) << 2)  | ((word >> 2)  & 0x33333333); // swap adjacent pairs
                            word = ((word & 0x0F0F0F0F) << 4)  | ((word >> 4)  & 0x0F0F0F0F); // swap nibbles
                            word = ((word & 0x00FF00FF) << 8)  | ((word >> 8)  & 0x00FF00FF); // swap bytes
                            word = ((word & 0x0000FFFF) << 16) | ((word >> 16) & 0x0000FFFF); // swap 16 bit chunks
                            //complement    
                            word ^= 0xAAAAAAAA;
                            //shift 
                            word >>= 2*(16-m_word_size);
                        }
                    }
                }
            }

            /*
            {
                lock_guard<mutex> guard(m_out_mutex);
                for(int l = 0; l < kmers_for_pos.size(); ++l)
                    cerr << "Word hits: " << acc << " " << l << " " << kmers_for_pos[l].size() << endl;
            }
            */

            {
                int ksize = kmers_for_pos.size();
                vector<unsigned> kcount(ksize);
                for(int p = 0; p < ksize; ++p)
                    kcount[p] = kmers_for_pos[p].size();
                for(int p = 0; p < ksize; ++p) {
                    if(kcount[p] >= m_kmer_complexity) {
                        for(int l = p; l < p+m_kmer_len/3; ++l)
                            target[l] = 'X';
                        for(int l = max(0,p-m_kmer_len/3+1); l < min(p+m_kmer_len/3,ksize); ++l)
                            kmers_for_pos[l].clear(); 
                    }
                }
            }
            int simpl_kmers = 0;
            for(char c : target) {
                if(c == 'X')
                    ++simpl_kmers;
            }

            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                auto& kfp = kmers_for_pos[p];
                int mcount = max_counts_for_pos[p];
                for(auto it_loop = kfp.begin(); it_loop != kfp.end(); ) {
                    auto it = it_loop++;
                    auto abundance = m_graphp->Abundance(*it);
                    if(abundance == 1 || abundance < m_fraction*mcount) 
                        kfp.erase(it);
                }
            }

            unordered_map<Node, list<int>, Node::Hash> positions;
            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                for(Node node : kmers_for_pos[p]) 
                    positions[node].push_back(p);
            }
            set<int> periodic_points;
            for(auto pos : positions) {
                if(pos.second.size() > 1) {
                    for(int p : pos.second)
                        periodic_points.insert(p);
                }
            }

            int seed_positions = 0;
            unordered_map<Node, int, Node::Hash> target_kmers; // [node] position on target
            for(int p = 0; p < (int)kmers_for_pos.size(); ++p) {
                if(periodic_points.count(p) || kmers_for_pos[p].empty())
                    continue;

                ++seed_positions;
                for(Node node : kmers_for_pos[p])
                    target_kmers.emplace(node, p);                
            }

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Periodic points for : " << acc << " " << periodic_points.size() << " " << seed_positions << " " << target_len << "/" << simpl_kmers << "\n";
                cerr << "Seed kmers for " << acc << " " << target_kmers.size() << " in " << timer.Elapsed();
                timer.Restart();
            } 

            CGuidedGraph gugraph(m_kmer_len, m_secondary_kmer_len);

            int total_forks = 0;
            int aligned_forks = 0;
            int total_length = 0;
            int aligned_length = 0;

            while(!target_kmers.empty()) {
                unordered_map<Node, int, Node::Hash>::iterator first_kmerp;
                typedef unordered_map<Node, int, Node::Hash>::value_type elem_t;
                first_kmerp = min_element(target_kmers.begin(), target_kmers.end(), 
                                          [this,target_len](const elem_t& a, const elem_t& b) 
                                          { 
                                              int dista = abs((target_len-m_kmer_len)/2-a.second);
                                              int distb = abs((target_len-m_kmer_len)/2-b.second);
                                              if(dista != distb)
                                                  return dista < distb;
                                              else if(m_graphp->Abundance(a.first) != m_graphp->Abundance(b.first))
                                                  return m_graphp->Abundance(a.first) > m_graphp->Abundance(b.first);
                                              else
                                                  return m_graphp->GetNodeSeq(a.first) < m_graphp->GetNodeSeq(b.first);
                                          });
                    
                Node initial_node = first_kmerp->first;
                string kmer_seq = m_graphp->GetNodeSeq(initial_node);
                int first_matching_kmer = first_kmerp->second; // position on target
                target_kmers.erase(first_kmerp);

                if(!CheckSeed(initial_node, 3*first_matching_kmer, 3*target_len-3*first_matching_kmer-m_kmer_len))
                    continue;

                bool alts = (first_matching_kmer == 0 && tstart);
                string kmer_prot = m_genetic_code.Translate(m_graphp->GetNodeKmer(initial_node), m_kmer_len, alts);
                int left_not_aligned = 0;
                int right_not_aligned = 0;
                int left_penalty = 0;
                int right_penalty = 0;
                {
                    int lminscore = numeric_limits<int>::max();
                    int lminpos = -1;
                    int score = 0;
                    for(int p = 0; p < m_kmer_len/3; ++p) {
                        score += m_delta.matrix[(int)kmer_prot[p]][(int)target[first_matching_kmer+p]];
                        if(score <= lminscore) {
                            lminscore = score;
                            lminpos = p;
                        }
                    }
                    int rminscore = numeric_limits<int>::max();
                    int rminpos = -1;
                    score = 0;
                    for(int p = 0; p < m_kmer_len/3; ++p) {
                        score += m_delta.matrix[(int)kmer_prot[m_kmer_len/3-1-p]][(int)target[first_matching_kmer+m_kmer_len/3-1-p]];
                        if(score <= rminscore) {
                            rminscore = score;
                            rminpos = p;
                        }
                    }
                    if(lminscore <= 0) {
                        left_not_aligned = 3*(lminpos+1);
                        left_penalty = -lminscore;
                    }
                    if(rminscore <= 0) {
                        right_not_aligned = 3*(rminpos+1);
                        right_penalty = -rminscore;
                    }
                }
                if(left_not_aligned+right_not_aligned > 0.5*m_kmer_len)
                    continue;

                gugraph.StartNewAssembly(kmer_seq, left_not_aligned, right_not_aligned);
                total_length += m_kmer_len;
                aligned_length += m_kmer_len;

                size_t assembled_len = 0;
                if(first_matching_kmer > 0) {                  // needs left extension
                    string left_target = target.substr(0, first_matching_kmer);
                    reverse(left_target.begin(), left_target.end());
                    unique_ptr<CGuidedPathBase> extenderp;
                    if(m_allow_frameshifts)
                        extenderp.reset(new CGuidedPathFS(m_graphp->ReverseComplement(initial_node), left_penalty, left_not_aligned, m_protect_ends, left_target, 
                                                          target_len-first_matching_kmer-m_kmer_len/3, *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, 
                                                          m_gap_open, m_gap_extend, m_fs_open, m_drop_off, anchor_frac, m_genetic_code, false, m_secondary_kmer_threshold));
                    else
                        extenderp.reset(new CGuidedPathAA(m_graphp->ReverseComplement(initial_node), left_penalty, left_not_aligned, m_protect_ends, left_target, 
                                                          target_len-first_matching_kmer-m_kmer_len/3, *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, 
                                                          m_gap_open, m_gap_extend, m_drop_off, anchor_frac, m_genetic_code, false, m_secondary_kmer_threshold));
                    map<typename CGuidedGraph::TAnchor, int> lanchors;
                    map<Node, int> lnotaligned; // not used
                    int total = 0;
                    deque<int> forks(m_buf_length+1);
                    while(extenderp->ProcessNextEdge()) {
                        if(extenderp->PathEnd()) {
                            SPathChunk chunk = extenderp->GetLastSegment();
                            total += chunk.m_seq.size();
                            total_length += chunk.m_seq.size();
                            int al = chunk.m_seq.size()-extenderp->NotAligned();
                            if(al > 0)
                                aligned_length += al;
                            gugraph.AddLeftSegment(chunk, lanchors, lnotaligned, gugraph.End()); 
                            for(int p = extenderp->StartingShift(); p < extenderp->AssembledSeqLength()-extenderp->NotAligned(); ++p)
                                target_kmers.erase(m_graphp->ReverseComplement(extenderp->AssembledSeq()[p].m_node));
                            lanchors.clear();
                            extenderp->DeleteNotAlignedForks(extenderp->AssembledSeqLength()-extenderp->NotAligned()+m_kmer_len);
                            gugraph.RewindLeftBranch(extenderp->StartingShift(), anchor_frac);
                        } else if(extenderp->NotAligned()%3 == 0 && extenderp->LastStepNode().isValid() && extenderp->SolidKmer()) {
                            typename CGuidedGraph::TAnchor anchor(extenderp->LastStepNode(), first_matching_kmer-1-extenderp->GetMaxPos());  // left end of alignment on target
                            int chunk_len = extenderp->AssembledSeqLength()-extenderp->StartingShift();
                            lanchors[anchor] = chunk_len-1;
                            CGuidedGraph::TSegmentP hook = gugraph.KnownLeftAnchor(anchor);
                            if(hook != gugraph.End()) {
                                SPathChunk chunk = extenderp->GetLastSegment();
                                total += chunk.m_seq.size();
                                total_length += chunk.m_seq.size();
                                aligned_length += chunk.m_seq.size();
                                gugraph.AddLeftSegment(chunk, lanchors, lnotaligned, hook);
                                for(int p = extenderp->StartingShift(); p < extenderp->AssembledSeqLength()-extenderp->NotAligned(); ++p) 
                                    target_kmers.erase(m_graphp->ReverseComplement(extenderp->AssembledSeq()[p].m_node));
                                lanchors.clear();
                                extenderp->DeleteLastBranch();
                                gugraph.RewindLeftBranch(extenderp->StartingShift(), anchor_frac);
                            }
                        }                                                                                                                                                          

                        forks.pop_front();
                        forks.push_back(extenderp->ForkCount());
                        if(total > 10*target_len*3 || (total > 0 && m_min_bases_per_fork*(forks.back()-forks.front()) > m_buf_length)) {                            
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted assemblyA " << acc << " kmer " << kmer_seq << " " << first_matching_kmer << endl;                            
                            break;                        
                        }                                                
                    }
                    assembled_len += total;
                    total_forks += extenderp->ForkCount();
                    aligned_forks += extenderp->AlignedForkCount();
                }

                gugraph.CleanBranch();
                if(first_matching_kmer < target_len-m_kmer_len/3) {                        
                    string right_target =  target.substr(first_matching_kmer+m_kmer_len/3);
                    unique_ptr<CGuidedPathBase> extenderp;
                    if(m_allow_frameshifts)
                        extenderp.reset(new CGuidedPathFS(initial_node, right_penalty, right_not_aligned, m_protect_ends, right_target, 
                                                          first_matching_kmer, *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, m_gap_open, 
                                                          m_gap_extend, m_fs_open, m_drop_off, anchor_frac, m_genetic_code, true, m_secondary_kmer_threshold));
                    else
                        extenderp.reset(new CGuidedPathAA(initial_node, right_penalty, right_not_aligned, m_protect_ends, right_target, 
                                                          first_matching_kmer, *m_graphdiggerp, *m_secondary_graphdiggerp, m_delta, m_gap_open, 
                                                          m_gap_extend, m_drop_off, anchor_frac, m_genetic_code, true, m_secondary_kmer_threshold));
                    map<typename CGuidedGraph::TAnchor, int> ranchors;
                    map<Node, int> rnotaligned;  // not used
                    int total = 0;
                    deque<int> forks(m_buf_length+1);
                    while(extenderp->ProcessNextEdge()) {
                        if(extenderp->PathEnd()) {
                            SPathChunk chunk = extenderp->GetLastSegment();
                            total += chunk.m_seq.size();
                            total_length += chunk.m_seq.size();
                            int al = chunk.m_seq.size()-extenderp->NotAligned();
                            if(al > 0)
                                aligned_length += al;
                            gugraph.AddRightSegment(chunk, ranchors, rnotaligned, gugraph.End());
                            for(int p = extenderp->StartingShift(); p < extenderp->AssembledSeqLength()-extenderp->NotAligned(); ++p)
                                target_kmers.erase(extenderp->AssembledSeq()[p].m_node);                            
                            ranchors.clear();
                            extenderp->DeleteNotAlignedForks(extenderp->AssembledSeqLength()-extenderp->NotAligned()+m_kmer_len);
                            gugraph.RewindRightBranch(extenderp->StartingShift(), anchor_frac);
                        } else if(extenderp->NotAligned()%3 == 0 && extenderp->LastStepNode().isValid() && extenderp->SolidKmer()) {
                            typename CGuidedGraph::TAnchor anchor(extenderp->LastStepNode(), first_matching_kmer+m_kmer_len/3+extenderp->GetMaxPos());  // right end of alignment on target
                            int chunk_len = extenderp->AssembledSeqLength()-extenderp->StartingShift();
                            ranchors[anchor] = chunk_len-1; 
                            CGuidedGraph::TSegmentP hook = gugraph.KnownRightAnchor(anchor);
                            if(hook != gugraph.End()) {
                                SPathChunk chunk = extenderp->GetLastSegment();
                                total += chunk.m_seq.size();
                                total_length += chunk.m_seq.size();
                                aligned_length += chunk.m_seq.size();
                                gugraph.AddRightSegment(chunk, ranchors, rnotaligned, hook);
                                for(int p = extenderp->StartingShift(); p < extenderp->AssembledSeqLength()-extenderp->NotAligned(); ++p) 
                                    target_kmers.erase(extenderp->AssembledSeq()[p].m_node);                                
                                ranchors.clear();
                                extenderp->DeleteLastBranch();
                                gugraph.RewindRightBranch(extenderp->StartingShift(), anchor_frac);
                            }                                                                                                
                        }                                                
                        forks.pop_front();
                        forks.push_back(extenderp->ForkCount());
                        if(total > 10*target_len*3 || (total > 0 && m_min_bases_per_fork*(forks.back()-forks.front()) > m_buf_length)) {                            
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted assemblyB " << acc << " kmer " << kmer_seq << " " << first_matching_kmer << endl;                            
                            break;                        
                        }                                               
                    }
                    assembled_len += total;
                    total_forks += extenderp->ForkCount();
                    aligned_forks += extenderp->AlignedForkCount();
                }
                gugraph.RemoveNotAlignedSegments(anchor_frac);
            }

            GFAGraph gfa_graph(acc, m_kmer_len);
            gugraph.GetGFAGraph(gfa_graph);
            gfa_graph.CalculateChainLength();
            if(!m_allow_frameshifts) {
                for(auto& seg : gfa_graph)
                    seg.m_frame = seg.m_left_len%3;
            }
            gfa_graph.AssignGroupNumber();
            gfa_graph.TrimGroups(m_graph_uniformity, min_len);  
            gfa_graph.MergeRedundantLinks();
            gfa_graph.AssignGroupNumber();            

            TGFACollection splitted = gfa_graph.SplitGroups();
            for(auto& graph : splitted) {
                if(m_min_count == 1) {
                    if(graph.RemoveHair(*m_secondary_graphdiggerp, m_fraction))
                        graph.MergeForks();                        
                    graph.RemoveSinglReadSnps(*m_secondary_graphp);
                }
                if(m_remove_homopolymer_indels)
                    graph.RemoveHomopolymerIndels(*m_secondary_graphp, m_homopolymer_ratio, m_homopolymer_len);
                if(m_no_reads && m_no_pairs) {
                    if(graph.CollapsFreeEnds() > 0)
                        graph.MergeRedundantLinks();
                    if(m_extend_ends) {
                        graph.GenerateKmers(*m_graphp);                
                        graph.ExtendToFirstFork(*m_graphdiggerp);                
                    }
                    graph.GenerateKmersAndScores(*m_secondary_graphp);
                } else {
                    graph.GenerateKmersAndScores(*m_graphp);
                }
                graph.ScoreGraph(get<2>(item.second), m_word_size);
            }
            rslts.splice(rslts.end(), splitted);

            {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Finished assembling: " << acc << " Forks: " << aligned_forks << "/" << total_forks << " Length: " << aligned_length << "/" << total_length << " in " << timer.Elapsed();
            }
        } 
    }
                           
    SMatrix m_delta;  //blosum62
    GeneticCode m_genetic_code;
    bool m_translate_na;
    int m_fs_open;
    bool m_allow_frameshifts;
};
    
} // namespace
#endif /* _GuidedAssembler_ */

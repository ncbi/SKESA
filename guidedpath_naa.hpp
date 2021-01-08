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

#ifndef _GuidedPathNAA_
#define _GuidedPathNAA_

#include "graphdigger.hpp"
#include "glb_align.hpp"
#include "genetic_code.hpp"
#include <stack> 

using namespace std;
namespace DeBruijn {

    struct SBranch {
        vector<int> m_sm;    // best scores in previs a-row
        vector<int> m_gapb;  // best score with b-gap
        Node m_node;         // last node
        int m_na;
        int m_maxscore;
        int m_maxposa;
        int m_maxposb;
        int m_jmin;          // b-interval evaluated using prevous row results
        int m_jmax;
        int m_isfork;
    };

    struct SBranchFS {
        array<vector<int>,4> m_s;    // best scores in 4 last raws
        array<vector<int>,4> m_gapb; // best gaps (any length) in b (insertions in a)
        Node m_node;                 // last node
        int m_na;
        int m_maxscore;
        int m_maxposa;
        int m_maxposb;
        int m_isfork;
        void Rotate() {
            rotate(m_s.begin(), m_s.end()-1, m_s.end());
            rotate(m_gapb.begin(), m_gapb.end()-1, m_gapb.end());
        }
    };

    struct SPathBase {
        SPathBase(char c = 0) : m_nt(c) {}
        int m_fork = eNoFork;              // eRightFork - base - 1 is right fork
                                           // eLeftFork  - base - (kmer-1) is left fork
        char m_nt;
    };
    struct SPathBaseKmer : public SPathBase {
        Node m_node;
    };
    struct SPathChunk {
        SPathChunk() : m_score(numeric_limits<int>::min()), m_tlen(0), m_starting_shift(0) {}
        int m_score;
        int m_tlen;
        int m_starting_shift;
        int m_not_aligned;
        vector<SPathBase> m_seq;
    };

    class CGuidedPathBase {
    public:
        virtual ~CGuidedPathBase() = default; 
        virtual bool ProcessNextEdge() = 0;
        virtual void DeleteLastBranch() = 0;
        virtual void DeleteNotAlignedForks(int starting_shift) = 0;
        virtual SPathChunk GetBestPart() const = 0;
        virtual SPathChunk GetLastSegment() const = 0;
        virtual bool SolidKmer() = 0;
        virtual int StartingShift() const = 0;
        virtual int GetMaxPos() const = 0;
        virtual const Node& LastStepNode() const = 0;
        virtual int AssembledSeqLength() const = 0;
        virtual int NotAligned() const = 0;
        virtual const vector<SPathBaseKmer>& AssembledSeq() const = 0;
        virtual bool PathEnd() const = 0;
        virtual int Num() const = 0;
        virtual int ForkCount() const = 0;
        virtual int AlignedForkCount() const = 0;

    protected:
        virtual bool AddOneBase(const Successor& c) = 0;
    };

    template <class Branch>
    class CGuidedPath : public CGuidedPathBase {
    public:
        CGuidedPath(const Node& initial_node, int initial_penalty, int initial_not_aligned, bool protect_ends, const string& target_extension, int position_on_target,  GraphDigger& graph_digger, 
                    GraphDigger& secondary_graph_digger, SMatrix& delta, int gapopen, int gapextend, int dropoff, double anchor_frac, int secondary_kmer_threshold) : 
            m_b(target_extension), m_position_on_target(position_on_target), m_fork_count(0), m_aligned_fork_count(0), m_path_end(false), m_num(0), m_starting_shift(0), m_initial_penalty(initial_penalty), 
            m_initial_not_aligned(initial_not_aligned), m_negative_limit(0), m_protect_ends(protect_ends), m_anchor_frac(anchor_frac), m_graph_digger(graph_digger), m_secondary_graph_digger(secondary_graph_digger),
            m_delta(delta), m_rho(gapopen), m_sigma(gapextend), m_dropoff(dropoff), m_is_prot(false), m_initial_node(m_graph_digger.Graph().GetNodeSeq(initial_node)), m_secondary_kmer_threshold(secondary_kmer_threshold) {}

        bool ProcessNextEdge() {
            if(m_edges.empty())
                return false;            
        
            ++m_num;

            if(m_path_end) 
                DeleteLastBranch();

            m_branch.m_node = m_edges.top().m_node;
            m_path_end = !AddOneBase(m_edges.top()) || (m_branch.m_maxposb == (int)m_b.size()-1);
            m_edges.pop();

            m_branch.m_isfork = eNoFork;

            int kmer_len = m_graph_digger.Graph().KmerLen();
            int remaining_len = m_b.size()-GetMaxPos()-1;
            if(m_is_prot)
                remaining_len *= 3;
            remaining_len -= NotAligned();
            int target_len = m_b.size()+m_position_on_target;
            if(m_is_prot)
                target_len *= 3;
            target_len += kmer_len; // target size in bp
            int back_step_len = target_len-(remaining_len+kmer_len-1); // 1 because we will make a 1bp step
            int margin = max(100, kmer_len);
            bool check_forward = (!m_protect_ends || margin < remaining_len);
            bool check_backward = (!m_protect_ends || margin < back_step_len);

            int fork_info;
            vector<Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessorsF(m_branch.m_node, &fork_info, check_forward, check_backward); 

            bool need_secondary = neighbors.empty();
            if(neighbors.size() == 1 && (m_graph_digger.LowCount() == 1 || m_secondary_kmer_threshold > 1))
                need_secondary = m_graph_digger.Graph().Abundance(neighbors.front().m_node) <= m_secondary_kmer_threshold;

            int secondary_kmer_len = m_secondary_graph_digger.Graph().KmerLen();
            int assembled_len = m_a.size();

            if(need_secondary && secondary_kmer_len < kmer_len && NotAligned() < kmer_len && assembled_len >= kmer_len) {
                string long_kmer;
                if(assembled_len < kmer_len)
                    long_kmer = m_initial_node.substr(assembled_len);
                for(int k = max(0,assembled_len-kmer_len); k < assembled_len; ++k)
                    long_kmer.push_back(m_a[k].m_nt);

                bool good_extension = true;
                deque<Node> nodes;
                CReadHolder rh(false);
                rh.PushBack(long_kmer);
                for(CReadHolder::kmer_iterator ik = rh.kbegin(secondary_kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first    
                    nodes.push_front(m_secondary_graph_digger.Graph().GetNode(*ik));
                for(int i = 0; i < (int)nodes.size()-1 && good_extension; ++i) {
                    int rl = remaining_len+kmer_len-secondary_kmer_len+i;
                    int bcl = target_len-(rl+secondary_kmer_len-1); // 1 because we will make a 1bp step    
                    bool check_frwd = (!m_protect_ends || margin < rl);
                    bool check_bckwd = (!m_protect_ends || margin < bcl);
                    vector<Successor> nbrs = m_secondary_graph_digger.GetReversibleNodeSuccessorsF(nodes[i], nullptr, check_frwd, check_bckwd);
                    good_extension = false;
                    for(auto& nbr : nbrs) {
                        if(nbr.m_node == nodes[i+1]) {
                            good_extension = true;
                            break;
                        }
                    }
                } 
                if(good_extension) {
                    string secondary_kmer;
                    if(assembled_len < secondary_kmer_len)
                        secondary_kmer = m_initial_node.substr(kmer_len-(secondary_kmer_len-assembled_len));
                    for(int k = max(0,assembled_len-secondary_kmer_len); k < assembled_len; ++k)
                        secondary_kmer.push_back(m_a[k].m_nt);
                    Node secondary_node = m_secondary_graph_digger.Graph().GetNode(secondary_kmer);
                    int rl = remaining_len;
                    int bcl = target_len-(rl+secondary_kmer_len-1); // 1 because we will make a 1bp step
                    bool check_frwd = (!m_protect_ends || margin < rl);
                    bool check_bckwd = (!m_protect_ends || margin < bcl);
                    vector<Successor>  secondary_neighbors = m_secondary_graph_digger.GetReversibleNodeSuccessorsF(secondary_node, &fork_info, check_frwd, check_bckwd);
                    if(!secondary_neighbors.empty()) {
                        neighbors = secondary_neighbors;
                        if(fork_info & eLeftFork)
                            fork_info |= eSecondaryKmer;
                        for(auto& neighbor : neighbors) {
                            string kmer = long_kmer.substr(1);
                            kmer.push_back(neighbor.m_nt);
                            neighbor.m_node = m_graph_digger.Graph().GetNode(kmer);
                        }                                                
                    }
                }
            }
           
            if(neighbors.empty())
                m_path_end = true;
            else
                m_branch.m_isfork |= fork_info;
                
            if(m_path_end)
                return true;

            for(int i = (int)neighbors.size()-1; i >= 0; --i)
                m_edges.push(neighbors[i]);
            if(neighbors.size() > 1) {
                int fcount = neighbors.size()-1;
                m_forks.push(make_pair(m_branch, fcount));
                ++m_fork_count;
                if(NotAligned() < kmer_len)
                    ++m_aligned_fork_count;
            }            

            return true;
        }
        void DeleteNotAlignedForks(int starting_shift) {
            while(!m_forks.empty() && m_forks.top().first.m_na > starting_shift) {
                int neigbors_num = m_forks.top().second;
                while(neigbors_num-- > 0)
                    m_edges.pop();    
                m_forks.pop();
            }

            if(!m_forks.empty()) {
                m_path_end = false;
                m_branch = m_forks.top().first;
                m_a.erase(m_a.begin()+m_branch.m_na, m_a.end());
                m_starting_shift = m_branch.m_na;
                if(--m_forks.top().second == 0)
                    m_forks.pop(); 
            } else {
                m_a.clear();
            }
        }
        void DeleteLastBranch() {
            if(m_edges.empty())
                return;    

            if(!m_path_end) {
                int neigbors_num = 1;
                if(!m_forks.empty() && m_forks.top().first.m_na == m_branch.m_na) {  // last base is a fork 
                    neigbors_num += m_forks.top().second;
                    m_forks.pop();
                }
                while(neigbors_num-- > 0)
                    m_edges.pop();    
            }

            if(!m_forks.empty()) {
                m_path_end = false;
                m_branch = m_forks.top().first;
                m_a.erase(m_a.begin()+m_branch.m_na, m_a.end());
                m_starting_shift = m_branch.m_na;
                if(--m_forks.top().second == 0)
                    m_forks.pop(); 
            } else {
                m_a.clear();
            }
        }
        SPathChunk GetBestPart() const { // whole sequence
            SPathChunk rslt;
            rslt.m_score = m_branch.m_maxscore;
            rslt.m_tlen = m_branch.m_maxposb+1;
            rslt.m_starting_shift = m_starting_shift;
            rslt.m_not_aligned = 0;
            if(m_branch.m_maxscore == 0)
                rslt.m_not_aligned = m_initial_not_aligned;
            rslt.m_seq.insert(rslt.m_seq.end(), m_a.begin(), m_a.begin()+m_branch.m_maxposa+1);

            return rslt;
        }
        SPathChunk GetLastSegment() const { // last segment after fork
            SPathChunk rslt;
            rslt.m_score = m_branch.m_maxscore;
            rslt.m_tlen = m_branch.m_maxposb+1;
            rslt.m_starting_shift = m_starting_shift;
            rslt.m_not_aligned = m_branch.m_na-1-m_branch.m_maxposa;
            if(m_branch.m_maxscore == 0)
                rslt.m_not_aligned += m_initial_not_aligned;
            rslt.m_seq.insert(rslt.m_seq.end(), m_a.begin()+m_starting_shift, m_a.end());

            return rslt;
        }
        bool SolidKmer() { 
            int kmer_len = m_graph_digger.Graph().KmerLen();
            if(!m_branch.m_node.isValid() || m_branch.m_na < kmer_len || m_branch.m_na-1-m_branch.m_maxposa >= m_anchor_frac*kmer_len)
                return false;

            return true;
        }
        int StartingShift() const { return m_starting_shift; }
        int GetMaxPos() const { return m_branch.m_maxposb; }
        const Node& LastStepNode() const { return m_branch.m_node; }
        int AssembledSeqLength() const { return m_a.size(); }
        int NotAligned() const {
            int not_aligned = m_branch.m_na-1-m_branch.m_maxposa;
            if(m_branch.m_maxscore == 0)
                not_aligned += m_initial_not_aligned;
            return not_aligned;
        }
        const vector<SPathBaseKmer>& AssembledSeq() const { return m_a; }

        bool PathEnd() const { return m_path_end; }
        int Num() const { return m_num; }
        int ForkCount() const { return m_fork_count; }
        int AlignedForkCount() const { return m_aligned_fork_count; }

    protected:

        Branch m_branch;
        vector<SPathBaseKmer> m_a;
        string m_b;
        int m_position_on_target;
        stack<pair<Branch,int> > m_forks;
        stack<Successor> m_edges;
        int m_fork_count;
        int m_aligned_fork_count;
        bool m_path_end;
        int m_num;
        int m_starting_shift;
        int m_initial_penalty;
        int m_initial_not_aligned;
        int m_negative_limit;
        bool m_protect_ends;

        double m_anchor_frac;
        GraphDigger m_graph_digger; 
        GraphDigger m_secondary_graph_digger; 
        SMatrix& m_delta;
        int m_rho;
        int m_sigma;
        int m_dropoff;

        bool m_is_prot;
        string m_initial_node;

        int m_secondary_kmer_threshold;
    };

    class CGuidedPathNA : public CGuidedPath<SBranch> {
    public:
        CGuidedPathNA(const Node& initial_node, int initial_penalty, int initial_not_aligned, bool protect_ends, const string& target_extension, int position_on_target, GraphDigger& graph_digger, GraphDigger& secondary_graph_digger, SMatrix& delta, 
                      int gapopen, int gapextend, int dropoff, double anchor_frac, int secondary_kmer_threshold) : 
            CGuidedPath(initial_node, initial_penalty, initial_not_aligned, protect_ends, target_extension, position_on_target, graph_digger, secondary_graph_digger, delta, gapopen, gapextend, dropoff, anchor_frac, secondary_kmer_threshold) {

            int kmer_len = m_graph_digger.Graph().KmerLen();
            m_negative_limit = min(m_delta.matrix['A']['C']*kmer_len, -m_rho-m_sigma*kmer_len);

            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            m_branch.m_sm.resize(nb+1,bignegative);
            m_branch.m_gapb.resize(nb+1,bignegative);

            m_branch.m_jmin = 0;
            m_branch.m_jmax = nb-1;
            m_branch.m_sm[0] = 0;
            m_branch.m_sm[1] = -m_rho-m_sigma;                                          // scores for --------------      
            for(int i = 2; i <= nb && m_branch.m_sm[i-1]-m_sigma > -2*m_dropoff; ++i) { //            BBBBBBBBBBBBBB      
                m_branch.m_sm[i] = m_branch.m_sm[i-1]-m_sigma;
                m_branch.m_jmax = min(nb-1,i);
            }    
                
            m_branch.m_maxscore = 0;
            m_branch.m_maxposa = -1;
            m_branch.m_maxposb = -1;
            m_branch.m_na = 0;

            m_branch.m_node = initial_node;
            m_branch.m_isfork = eNoFork;
            if(m_branch.m_node.isValid()) {
                int remaining_len = m_b.size();
                int back_step_len = m_position_on_target+1; // 1 because we will make a 1bp step
                int margin = max(100, kmer_len);
                bool check_forward = (!m_protect_ends || margin < remaining_len);
                bool check_backward = (!m_protect_ends || margin < back_step_len);
                int fork_info;
                vector<Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessorsF(m_branch.m_node, &fork_info, check_forward, check_backward);            
                if(!neighbors.empty())
                    m_branch.m_isfork |= fork_info;

                for(auto& neighbor : neighbors)
                    m_edges.push(neighbor);
                if(neighbors.size() > 1) {
                    int fcount = neighbors.size()-1;
                    m_forks.push(make_pair(m_branch, fcount));
                    ++m_fork_count;
                    ++m_aligned_fork_count;
                }
            }                    
            m_path_end = m_edges.empty();
        }

    private:
        bool AddOneBase(const Successor& c) {
            m_a.emplace_back();
            m_a.back().m_node = c.m_node;
            m_a.back().m_nt = c.m_nt;
            m_a.back().m_fork = m_branch.m_isfork;
            ++m_branch.m_na;
            return UpdateScore();
        }
        bool UpdateScore() {
            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            int rs = m_rho+m_sigma;
            int next_jmax = -1;
            int next_jmin = nb;

            vector<int> s(nb+1,bignegative);             // best scores in current a-raw    
            if(-m_rho-m_branch.m_na*m_sigma > m_branch.m_maxscore-2*m_dropoff) {
                next_jmin = 0;
                s[0] = -m_rho-m_branch.m_na*m_sigma;      // score for AAAAAAAAAAA  
            }                                             //           -----------  
            
            int gapa = bignegative;
            int ai = m_a[m_branch.m_na-1].m_nt;
            const char* matrix = (m_delta.matrix)[ai];
            int* sp = &s[m_branch.m_jmin];
            int smax = *sp;

            for(int j = m_branch.m_jmin; j <= m_branch.m_jmax; ) { // here j is 'real' position in b
                int ss = m_branch.m_sm[j]+matrix[(int)m_b[j]];

                gapa -= m_sigma;
                if(*sp-rs > gapa)
                    gapa = *sp-rs;
			
                int& gapbj = m_branch.m_gapb[++j]; // here j is one-shifted to account for extra element in vectors (nb+1)
                gapbj -= m_sigma;
                if(m_branch.m_sm[j]-rs > gapbj)
                    gapbj = m_branch.m_sm[j]-rs;

                if(gapa > gapbj) {
                    if(ss >= gapa) {
                        *(++sp) = ss;
                        if(ss-m_initial_penalty > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss-m_initial_penalty;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapa;
                    }
                } else {
                    if(ss >= gapbj) {
                        *(++sp) = ss;
                        if(ss-m_initial_penalty > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss-m_initial_penalty;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapbj;
                    }
                }
                
                if(*sp > m_branch.m_maxscore-2*m_dropoff) {
                    next_jmin = min(next_jmin, j-1);
                    next_jmax = min(nb-1,j);
                }                

                smax = max(smax, *sp);
            }
            swap(m_branch.m_sm,s);            
            // jmin never decreases
            m_branch.m_jmin = next_jmin;
            //right may decrease
            for(int l = next_jmax+1; l <= m_branch.m_jmax; ++l) {
                m_branch.m_gapb[l+1] = bignegative;
                m_branch.m_sm[l+1] = bignegative;
            }
            m_branch.m_jmax = next_jmax;            

            return smax-m_initial_penalty >= m_branch.m_maxscore-m_dropoff && m_branch.m_jmax >= m_branch.m_jmin && smax-m_initial_penalty > m_negative_limit;
        }
    };

    class CGuidedPathAA : public CGuidedPath<SBranch> {
    public:
        CGuidedPathAA(const Node& initial_node, int initial_penalty, int initial_not_aligned, bool protect_ends, const string& target_extension, int position_on_target, GraphDigger& graph_digger, GraphDigger& secondary_graph_digger, SMatrix& delta, 
                      int gapopen, int gapextend, int dropoff, double anchor_frac, const GeneticCode& genetic_code, bool forward, int secondary_kmer_threshold) : 
            CGuidedPath(initial_node, initial_penalty, initial_not_aligned, protect_ends, target_extension, position_on_target, graph_digger, secondary_graph_digger, delta, gapopen, gapextend, dropoff, anchor_frac, secondary_kmer_threshold), 
            m_genetic_code(genetic_code), m_forward(forward) { 

            m_is_prot = true;

            int kmer_len = m_graph_digger.Graph().KmerLen();            
            char min_element = numeric_limits<char>::max();
            for(int i = 0; i < 256; ++i) {
                for(int j = 0; j < i; ++j)
                    min_element = min(min_element, m_delta.matrix[i][j]);
            }
            m_negative_limit = min(min_element*kmer_len, -m_rho-m_sigma*kmer_len);

            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            m_branch.m_jmin = 0;
            m_branch.m_jmax = nb-1;
            m_branch.m_sm.resize(nb+1,bignegative);
            m_branch.m_gapb.resize(nb+1,bignegative);

            // scores for --------------
            //            BBBBBBBBBBBBBB
            m_branch.m_sm[0] = 0;
            m_branch.m_sm[1] = -m_rho-m_sigma;                  
            for(int i = 2; i <= nb; ++i)    
                m_branch.m_sm[i] = m_branch.m_sm[i-1]-m_sigma;
                
            m_branch.m_maxscore = 0;
            m_branch.m_maxposa = -1;
            m_branch.m_maxposb = -1;
            m_branch.m_na = 0;

            m_branch.m_node = initial_node;
            m_branch.m_isfork = eNoFork;
            if(m_branch.m_node.isValid()) {
                int remaining_len = 3*m_b.size();
                int back_step_len = 3*m_position_on_target+1; // 1 because we will make a 1bp step
                int margin = max(100, kmer_len);
                bool check_forward = (!m_protect_ends || margin < remaining_len);
                bool check_backward = (!m_protect_ends || margin < back_step_len);
                int fork_info;
                vector<Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessorsF(m_branch.m_node, &fork_info, check_forward, check_backward);            
                if(!neighbors.empty())
                    m_branch.m_isfork |= fork_info;

                for(auto& neighbor : neighbors)
                    m_edges.push(neighbor);
                if(neighbors.size() > 1) {
                    int fcount = neighbors.size()-1;
                    m_forks.push(make_pair(m_branch, fcount));
                    ++m_fork_count;
                    ++m_aligned_fork_count;
                }
            } 
            if(m_edges.empty())
                m_path_end = true;
        }

    private:
        bool AddOneBase(const Successor& c) {
            m_a.emplace_back();
            m_a.back().m_node = c.m_node;
            m_a.back().m_nt = c.m_nt;
            m_a.back().m_fork = m_branch.m_isfork;
            if((++m_branch.m_na)%3)
                return true;
            else
                return UpdateScore();
        }
        bool UpdateScore() {
            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            int rs = m_rho+m_sigma;

            vector<int> s(nb+1,bignegative);             // best scores in current a-raw 
            s[0] = -m_rho-m_branch.m_na*m_sigma;         // score for AAAAAAAAAAA
                                                         //           -----------  
           
            int gapa = bignegative;
            string codon(3, 0);
            if(m_forward) {
                codon[0] = m_a[m_branch.m_na-3].m_nt;
                codon[1] = m_a[m_branch.m_na-2].m_nt;
                codon[2] = m_a[m_branch.m_na-1].m_nt;
            } else {
                codon[0] = Complement(m_a[m_branch.m_na-1].m_nt);
                codon[1] = Complement(m_a[m_branch.m_na-2].m_nt);
                codon[2] = Complement(m_a[m_branch.m_na-3].m_nt);
            }
            int ai = m_genetic_code.AA(codon);
            const char* matrix = (m_delta.matrix)[ai];
            int* sp = &s[0];
            int smax = *sp;

            for(int j = 0; j < nb; ) { // here j is 'real' position in b
                int ss = m_branch.m_sm[j]+matrix[(int)m_b[j]];
                if(!m_forward && j == nb-1 && toupper(m_b[j]) == 'M' && m_genetic_code.IsStart(codon))
                    ss = m_branch.m_sm[j]+(m_delta.matrix)[(int)'M'][(int)'M'];

                gapa -= m_sigma;
                if(*sp-rs > gapa)
                    gapa = *sp-rs;
			
                int& gapbj = m_branch.m_gapb[++j]; // here j is one-shifted to account for extra element in vectors (nb+1)
                gapbj -= m_sigma;
                if(m_branch.m_sm[j]-rs > gapbj)
                    gapbj = m_branch.m_sm[j]-rs;

                if(gapa > gapbj) {
                    if(ss >= gapa) {
                        *(++sp) = ss;
                        if(ss-m_initial_penalty > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss-m_initial_penalty;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapa;
                    }
                } else {
                    if(ss >= gapbj) {
                        *(++sp) = ss;
                        if(ss-m_initial_penalty > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss-m_initial_penalty;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapbj;
                    }
                }

                smax = max(smax, *sp);
            }
            swap(m_branch.m_sm,s);  

            return smax-m_initial_penalty >= m_branch.m_maxscore-m_dropoff && smax-m_initial_penalty > m_negative_limit;
        }

        const GeneticCode& m_genetic_code;
        bool m_forward;
    };

    class CGuidedPathFS : public CGuidedPath<SBranchFS> {
    public:
        CGuidedPathFS(const Node& initial_node, int initial_penalty, int initial_not_aligned, bool protect_ends, const string& target_extension, int position_on_target, GraphDigger& graph_digger, GraphDigger& secondary_graph_digger, SMatrix& delta, 
                      int gapopen, int gapextend, int fsopen, int dropoff, double anchor_frac, const GeneticCode& genetic_code, bool forward, int secondary_kmer_threshold) : 
            CGuidedPath(initial_node, initial_penalty, initial_not_aligned, protect_ends, target_extension, position_on_target, graph_digger, secondary_graph_digger, delta, gapopen, gapextend, dropoff, anchor_frac, secondary_kmer_threshold), 
            m_rhofs(fsopen), m_genetic_code(genetic_code), m_forward(forward) { 

            m_is_prot = true;

            int kmer_len = m_graph_digger.Graph().KmerLen();            
            char min_element = numeric_limits<char>::max();
            for(int i = 0; i < 256; ++i) {
                for(int j = 0; j < i; ++j)
                    min_element = min(min_element, m_delta.matrix[i][j]);
            }
            m_negative_limit = min(min_element*kmer_len, -m_rho-m_sigma*kmer_len);

            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();

            // scores for na == 0,1,2
            for(int k = 0; k < 4; ++k) {      
                m_branch.m_s[k].resize(nb+1, bignegative);
                m_branch.m_gapb[k].resize(nb+1, bignegative);
            }
            m_branch.m_s[3][0] = -m_initial_penalty;                  // nothing aligned
            m_branch.m_s[3][1] = m_branch.m_s[3][0]-m_rho-m_sigma;    // 3*n deletion
            for(int i = 2; i <= nb; ++i)                              // extend 3*n deletion
                m_branch.m_s[3][i] = m_branch.m_s[3][i-1]-m_sigma;

            m_branch.m_s[2][0] = -m_initial_penalty-m_rhofs-m_sigma;  // first a base is insertion
            m_branch.m_s[2][1] = m_branch.m_s[2][0]-m_rho-m_sigma;    
            for(int i = 2; i <= nb; ++i)
                m_branch.m_s[2][i] = m_branch.m_s[2][i-1]-m_sigma;
            m_branch.m_gapb[2] = m_branch.m_s[2];

            m_branch.m_s[1][0] = -m_initial_penalty-m_rhofs-m_sigma;  // first two a bases is insertion
            m_branch.m_s[1][1] = m_branch.m_s[1][0]-m_rho-m_sigma;
            for(int i = 2; i <= nb; ++i)
                m_branch.m_s[1][i] = m_branch.m_s[1][i-1]-m_sigma;
            m_branch.m_gapb[1] = m_branch.m_s[1];

            m_branch.m_maxscore = 0;
            m_branch.m_maxposa = -1;
            m_branch.m_maxposb = -1;
            m_branch.m_na = 0;

            m_branch.m_node = initial_node;
            m_branch.m_isfork = eNoFork;
            if(m_branch.m_node.isValid()) {
                int remaining_len = 3*m_b.size();
                int back_step_len = 3*m_position_on_target+1; // 1 because we will make a 1bp step
                int margin = max(100, kmer_len);
                bool check_forward = (!m_protect_ends || margin < remaining_len);
                bool check_backward = (!m_protect_ends || margin < back_step_len);
                int fork_info;
                vector<Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessorsF(m_branch.m_node, &fork_info, check_forward, check_backward);            
                if(!neighbors.empty())
                    m_branch.m_isfork |= fork_info;

                for(auto& neighbor : neighbors)
                    m_edges.push(neighbor);
                if(neighbors.size() > 1) {
                    int fcount = neighbors.size()-1;
                    m_forks.push(make_pair(m_branch, fcount));
                    ++m_fork_count;
                    ++m_aligned_fork_count;
                }
            } 
            if(m_edges.empty())
                m_path_end = true;
        }
    
    private:
        bool AddOneBase(const Successor& c) {
            m_a.emplace_back();
            m_a.back().m_node = c.m_node;
            m_a.back().m_nt = c.m_nt;
            m_a.back().m_fork = m_branch.m_isfork;
            if(++m_branch.m_na < 3)
                return true;
            else
                return UpdateScore();
        }
        bool UpdateScore() {
            // bacause we don't need back trace and any gap extension is the same cost we may have one gap score as long as we penilize properly for the gap opening
            vector<int>& s = m_branch.m_s[0];         // best scores in current a-raw
            vector<int>& sm1 = m_branch.m_s[1];       // best scores in -1 a-raw
            vector<int>& sm2 = m_branch.m_s[2];       // best scores in -2 a-raw
            vector<int>& sm3 = m_branch.m_s[3];       // best scores in -3 a-raw
            vector<int>& gapb = m_branch.m_gapb[0];   // best score with b-gap (any insertion in a) in current a-raw
            vector<int>& gapbm3 = m_branch.m_gapb[3]; // best score with b-gap (any insertion in a) in -3 a-raw

            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            int rs = m_rho+m_sigma;
            int rsfs = m_rhofs+m_sigma;

            s[0] = (m_branch.m_na%3 ? -m_rhofs-m_sigma : -m_rho)-m_branch.m_na/3*m_sigma;  // score for AAAAAAAAAAA
                                                                                           //           -----------  
            string codon(3, 0);
            if(m_forward) {
                codon[0] = m_a[m_branch.m_na-3].m_nt;
                codon[1] = m_a[m_branch.m_na-2].m_nt;
                codon[2] = m_a[m_branch.m_na-1].m_nt;
            } else {
                codon[0] = Complement(m_a[m_branch.m_na-1].m_nt);
                codon[1] = Complement(m_a[m_branch.m_na-2].m_nt);
                codon[2] = Complement(m_a[m_branch.m_na-3].m_nt);
            }
            int ai = m_genetic_code.AA(codon);
            const char* matrix = (m_delta.matrix)[ai];
            int smax = s[0];

            int gapa = bignegative; // any deletion from a

            for(int j = 0; j < nb; ) {
                // b[j] current aa
                // any_score[j] score for previous aa
                // any_score[j+1] score for current aa
                int ss = sm3[j]+matrix[(int)m_b[j]];    // diagonal extension
                if(!m_forward && j == nb-1 && toupper(m_b[j]) == 'M' && m_genetic_code.IsStart(codon))
                    ss = sm3[j]+(m_delta.matrix)[(int)'M'][(int)'M'];

                gapa -= m_sigma;              // gapa extension
                if(s[j]-rs > gapa)
                    gapa = s[j]-rs;           // new 3bp deletion from a
                if(sm2[j]-rsfs > gapa)
                    gapa = sm2[j]-rsfs;       // new 1bp deletion from a
                if(sm1[j]-rsfs > gapa)
                    gapa = sm1[j]-rsfs;       // new 2bp deletion from a
                
                int& gapbj = gapb[++j];       // here j is one-shifted to account for extra element in vectors (nb+1)
                gapbj = gapbm3[j]-m_sigma;    // gapb extension
                if(sm3[j]-rs > gapbj)
                    gapbj = sm3[j]-rs;        // new 3bp insertion in a
                if(sm1[j]-rsfs > gapbj) 
                    gapbj = sm1[j]-rsfs;      // new 1bp insertion in a            
                if(sm2[j]-rsfs > gapbj)
                    gapbj = sm2[j]-rsfs;      // new 2bp insertion in a             

                s[j] = gapa;
                if(gapbj > s[j]) 
                    s[j] = gapbj;
                if(ss >= s[j]) {
                    s[j] = ss;
                    if(ss > m_branch.m_maxscore) {
                        m_branch.m_maxscore = ss;
                        m_branch.m_maxposa = m_branch.m_na-1;
                        m_branch.m_maxposb = j-1;
                    }
                }

                smax = max(smax, s[j]);
            }
            m_branch.Rotate();

            return smax >= m_branch.m_maxscore-m_dropoff && smax > m_negative_limit;
        }

        int m_rhofs;
        const GeneticCode& m_genetic_code;
        bool m_forward;
    };

} // namespace
#endif /* _GuidedPathNAA_ */

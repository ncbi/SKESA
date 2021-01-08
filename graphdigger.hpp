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

#ifndef _GraphDigger_
#define _GraphDigger_

#include <stack> 
#include <deque>
#include <forward_list>
#include <unordered_set>
#include "glb_align.hpp"
#include "DBGraph.hpp"

namespace DeBruijn {

/************************
General description

Class CDBGraphDigger, defined in this file, performs most of the actual assembling work. 

struct SContig is used to hold assembled sequences. Members are 
        deque<char> m_seq;               // sequence representing contig
        deque<CDBGraph::Node> m_kmers;   // all kmers of this contig sequence
        int m_kmer_len;                  // size of kmers used for building the contig
        CDBGraph::Node m_next_left;      // denied left kmer (connection possible but it is already owned)
        CDBGraph::Node m_next_right;     // denied right kmer (connection possible but it is already owned)
        SContig* m_left_link;            // if set points to 'left' contig
        SContig* m_right_link;           // if set points to 'right' contig
        int m_left_shift;                // shift+1 for m_next_left in this contig (positive for the right end)
        int m_right_shift;               // shift+1 for m_next_right in this contig (positive for the right end)
        int m_left_extend;               // number of newly assembled bases which could be clipped
        int m_right_extend;              // number of newly assembled bases which could be clipped
        SAtomic<uint8_t> m_is_taken;

There are three main scenarios how this structure may be created

1) From a previously assembled contig represented by a c++ string. 
   In this case, no members other than m_seq, m_kmers, and m_kmer_len change their default zero value. 

2) Assembled starting from one of the kmers which has not been used so far. Because assembly is done in multiple threads,
   two or more threads could start assembling the same contig from different starting kmers. At some point they will
   collide with each other and will try to obtain a kmer which has already been used by the other contig in a different
   thread. When this happens, the thread stops extending the contig and assigns the denied kmer to m_next_left or m_next_right.
   These partially assembled contigs (which internally are called fragments) could be connected to each other using
   m_next_left/m_next_right. It is done in ConnectFragments().

3) When we increase the kmer size, some previously assembled contigs could be extended or connected because
   the longer kmer could resolve some repeats. To achieve this, we assemble new contigs starting from each of the
   flank kmers. When these contigs are started, m_left_link/m_right_link are assigned to point to the
   parent contig and m_left_shift/m_right_shift are assigned to indicate the start position. Because the work is done
   in mutiple threads, the contigs could come in two fragments if they started from different contigs and come together.
   The connected fragments will have links on both sides. These contigs are 'connectors'. The rest of contigs are
   'extenders'. They are used by ConnectAndExtendContigs() to form a new contig set. There is an important corner case
   for connectors. A thread could start from contig A and finish assembling a connector all the way to contig B before
   some other threads starts dealing with B. In this case the sequence starting from B will not contain any real bases
   but will only have m_next_left/m_next_right and a link. Those are mentioned in the code as 'empty linkers' and should
   be treated as special cases. 

A note about multiprocessing in 2) and 3): In both cases, the threads should be able to decide which kmers or contigs are
still available for work. For communicating this information between the threads, the code uses lock-free c++ atomic
variables. For the kmers, this is stored in m_visited vector in CDBGraph. For the contigs, this is stored in m_is_taken.

The length of newly assembled sequence is stored in  m_left_extend/m_right_extend.

************************/


    mutex out_mutex;

    enum EForkType { eNoFork = 0, eLeftFork = 1, eRightFork = 2, eLeftBranch = 4, eRightBranch = 8, eSecondaryKmer = 16 };

    typedef vector<char> TVariation;
    struct SeqInterval {
        SeqInterval(TVariation::iterator b, TVariation::iterator e) : begin(b), end(e) {}
        bool operator<(const SeqInterval& other) const { return lexicographical_compare(begin, end, other.begin, other.end); }
        bool operator==(const SeqInterval& other) const { return equal(begin, end, other.begin); }
        
        TVariation::iterator begin;
        TVariation::iterator end;
    };
    typedef forward_list<TVariation> TLocalVariants;
    class CContigSequence : public vector<TLocalVariants> {
    public:
        int m_left_repeat = 0;     // number of bases which COULD be in repeat
        int m_right_repeat = 0;    // number of bases which COULD be in repeat 
        bool m_circular = false;

        size_t MemoryFootprint() const {
            size_t total = sizeof(CContigSequence)+sizeof(TLocalVariants)*capacity();
            for(auto& lst : *this) {
                for(auto& v : lst)
                    total += sizeof(TVariation)+sizeof(TVariation*)+v.capacity();
            }
            return total;
        }
        int VariantsNumber(int chunk) { return distance((*this)[chunk].begin(), (*this)[chunk].end()); }
        bool UniqueChunk(int chunk) const {
            auto it = (*this)[chunk].begin();
            return (it != (*this)[chunk].end() && ++it == (*this)[chunk].end());
        }
        bool VariableChunk(int chunk) const {
            auto it = (*this)[chunk].begin();
            return (it != (*this)[chunk].end() && ++it != (*this)[chunk].end());
        }
        size_t ChunkLenMax(int chunk) const {
            size_t mx = 0;
            for(auto& seq : (*this)[chunk])
                mx = max(mx, seq.size());
            return mx;
        }
        size_t ChunkLenMin(int chunk) const {
            size_t mn = numeric_limits<size_t>::max();
            for(auto& seq : (*this)[chunk])
                mn = min(mn, seq.size());            
            return mn;
        }
        size_t LenMax() const { 
            size_t len = 0;
            for(unsigned chunk = 0; chunk < size(); ++chunk)
                len += ChunkLenMax(chunk);
            return len; 
        }
        size_t LenMin() const { 
            size_t len = 0;
            for(unsigned chunk = 0; chunk < size(); ++chunk)
                len += ChunkLenMin(chunk);
            return len; 
        }

        void InsertNewVariant() { back().emplace_front(); }
        void InsertNewVariant(char c) { back().emplace_front(1, c); }
        template <typename ForwardIterator>
        void InsertNewVariant(ForwardIterator b, ForwardIterator e) { back().emplace_front(b, e); }

        void ExtendTopVariant(char c) { back().front().push_back(c); }
        template <typename ForwardIterator>
        void ExtendTopVariant(ForwardIterator b, ForwardIterator e) { back().front().insert(back().front().end(), b, e); }

        void InsertNewChunk() { emplace_back(); }
        template <typename ForwardIterator>
        void InsertNewChunk(ForwardIterator b, ForwardIterator e) {
            InsertNewChunk();
            InsertNewVariant(b, e);
        }

        void StabilizeVariantsOrder() {
            for(auto& chunk : *this) 
                chunk.sort();
        }
        void ReverseComplement() {
            std::swap(m_left_repeat, m_right_repeat);
            reverse(begin(), end());
            for(auto& chunk : *this) {
                for(auto& seq : chunk)
                    ReverseComplementSeq(seq.begin(), seq.end());
            }
            StabilizeVariantsOrder();
        }
        bool RemoveShortUniqIntervals(int min_uniq_len) {
            if(size() >= 5) {
                for(unsigned i = 2; i < size()-2; ) {
                    if((int)ChunkLenMax(i) < min_uniq_len) {
                        auto& new_chunk = *insert(begin()+i+2, TLocalVariants()); // empty chunk
                        for(auto& var1 : (*this)[i-1]) {
                            for(auto& var2 : (*this)[i+1]) {
                                new_chunk.push_front(var1);
                                auto& seq = new_chunk.front();
                                seq.insert(seq.end(), (*this)[i].front().begin(), (*this)[i].front().end());
                                seq.insert(seq.end(), var2.begin(), var2.end());
                            }
                        }
                        //erase i-1, i, i+1
                        erase(begin()+i-1, begin()+i+2);
                    } else {
                        i += 2;
                    }
                }
            }

            if(m_circular && size() >= 5 && (int)(ChunkLenMax(0)+ChunkLenMax(size()-1)) < min_uniq_len) {
                // rotate contig so that short interval is in the middle
                ExtendTopVariant(front().front().begin(), front().front().end()); // add first chunk to the last
                erase(begin());
                rotate(begin(), begin()+1, end()); // move first variable chunk to end
                InsertNewChunk();
                auto& seq = front().front();
                ExtendTopVariant(seq[0]); // move first base to the end
                seq.erase(seq.begin());
                RemoveShortUniqIntervals(min_uniq_len);
                return true;
            }

            return false;
        }
        void ContractVariableIntervals() {
            if(size() > 2) {
                for(unsigned i = 1; i < size()-1; ++i) {
                    if(VariableChunk(i)) {
                        auto& fseq = (*this)[i].front();
                        unsigned len = 0;
                        while(true) {
                            bool all_same = true;
                            for(auto& seq : (*this)[i]) {
                                if(seq.size() == len || seq[len] != fseq[len]) {
                                    all_same = false;
                                    break;
                                }
                            }
                            if(all_same)
                                ++len;
                            else
                                break;
                        }
                        if(len > 0) {
                            (*this)[i-1].front().insert((*this)[i-1].front().end(), fseq.begin(), fseq.begin()+len);
                            for(auto& seq : (*this)[i])
                                seq.erase(seq.begin(), seq.begin()+len);
                        }
                        len = 0;
                        while(true) {
                            bool all_same = true;
                            for(auto& seq : (*this)[i]) {
                                if(seq.size() == len || *(seq.rbegin()+len) != *(fseq.rbegin()+len)) {
                                    all_same = false;
                                    break;
                                }
                            }
                            if(all_same)
                                ++len;
                            else
                                break;
                        }
                        if(len > 0) {
                            (*this)[i+1].front().insert((*this)[i+1].front().begin(), fseq.end()-len, fseq.end());
                            for(auto& seq : (*this)[i])
                                seq.erase(seq.end()-len, seq.end());
                        }
                    }
                }
            }
        }

        bool AllSameL(int chunk, int shift) const {
            if(!VariableChunk(chunk))
                return false;
            
            auto Symbol_i = [](const TVariation& seq, const TVariation& next, unsigned i) {
                if(i < seq.size())
                    return seq[i];
                else if(i < seq.size()+next.size()-1 )  // -1 - we don't want to absorb all next
                    return next[i-seq.size()];
                else
                    return char(0);
            };

            char symb = Symbol_i((*this)[chunk].front(), (*this)[chunk+1].front(), shift);
            if(!symb)
                return false;

            auto it = (*this)[chunk].begin();
            for(++it; it != (*this)[chunk].end(); ++it) {
                if(Symbol_i(*it, (*this)[chunk+1].front(), shift) != symb)
                    return false;
            }
           
            return true;
        }

        bool AllSameR(int chunk, int shift) const {
            if(!VariableChunk(chunk))
                return false;
            
            auto Symbol_i = [](const TVariation& seq, const TVariation& prev, unsigned i) {
                if(i < seq.size())
                    return seq[seq.size()-1-i];
                else if(i < seq.size()+prev.size()-1)  // -1 - we don't want to absorb all prev
                    return prev[prev.size()+seq.size()-1-i];
                else
                    return char(0);
            };
            
            char symb = Symbol_i((*this)[chunk].front(), (*this)[chunk-1].front(), shift);
            if(!symb)
                return false;

            auto it = (*this)[chunk].begin();
            for(++it; it != (*this)[chunk].end(); ++it) {
                if(Symbol_i(*it, (*this)[chunk-1].front(), shift) != symb)
                    return false;
            }
           
            return true;
        }
        
        void IncludeRepeatsInVariableIntervals() {
            for(unsigned chunk = 1; chunk < size()-1; chunk += 2) {
                int min_len = ChunkLenMin(chunk);
                int len = 0;
                for(int shift = 0; AllSameL(chunk, shift); ++shift) {
                    if(shift >= min_len)
                        ++len;
                }
                if(len > 0) {
                    min_len += len;
                    auto& nseq = (*this)[chunk+1].front();
                    for(auto& seq : (*this)[chunk])
                        seq.insert(seq.end(), nseq.begin(), nseq.begin()+len);
                    nseq.erase(nseq.begin(), nseq.begin()+len);
                }
                len = 0;
                for(int shift = 0; AllSameR(chunk, shift); ++shift) {
                    if(shift >= min_len)
                        ++len;
                }
                if(len > 0) {
                    auto& pseq = (*this)[chunk-1].front();
                    for(auto& seq : (*this)[chunk])
                        seq.insert(seq.begin(), pseq.end()-len, pseq.end());
                    pseq.erase(pseq.end()-len, pseq.end());
                }
            }
        }               
    };

    typedef list<CContigSequence> TContigSequenceList;

    void CombineSimilarContigs(TContigSequenceList& contigs) {
        int match = 1;
        int mismatch = 2;
        int gap_open = 5;
        int gap_extend = 2;
        SMatrix delta(match, mismatch);

        list<string> all_variants;
        for(auto& contig : contigs) {
            list<string> variants;
            for(auto& seq : contig[0])
                variants.emplace_back(seq.begin(), seq.end());
            for(unsigned l = 1; l < contig.size(); ++l) {
                if(contig.UniqueChunk(l)) {
                    for(auto& seq : variants)
                        seq.insert(seq.end(), contig[l].front().begin(), contig[l].front().end());
                } else {
                    list<string> new_variants;
                    for(auto& seq : variants) {
                        for(auto& var : contig[l]) {
                            new_variants.push_back(seq);
                            new_variants.back().insert(new_variants.back().end(), var.begin(), var.end());
                        }
                    }
                    swap(variants, new_variants);
                }
            }
            all_variants.splice(all_variants.end(), variants);
        }
        all_variants.sort([](const string& a, const string& b) { return a.size() > b.size(); });

        list<list<TVariation>> all_groups;
        while(!all_variants.empty()) {
            string& query = all_variants.front();
            list<TVariation> group;
            auto it_loop = all_variants.begin();
            for(++it_loop; it_loop != all_variants.end(); ) {
                auto it = it_loop++;
                string& subject = *it;
                if(query.size()-subject.size() > 0.1*query.size())
                    continue;

                //                CCigar cigar = LclAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gap_open, gap_extend, delta.matrix);
                CCigar cigar = BandAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gap_open, gap_extend, delta.matrix, 0.1*query.size());
                if(cigar.QueryRange().first != 0 || cigar.QueryRange().second != (int)query.size()-1)
                    continue;
                if(cigar.SubjectRange().first != 0 || cigar.SubjectRange().second != (int)subject.size()-1)
                    continue;
                if(cigar.Matches(query.c_str(), subject.c_str()) < 0.9*query.size())
                    continue;

                TCharAlign align = cigar.ToAlign(query.c_str(), subject.c_str());
                if(group.empty()) {
                    group.emplace_back(align.first.begin(), align.first.end());
                    group.emplace_back(align.second.begin(), align.second.end());
                } else {
                    TVariation& master = group.front();
                    int mpos = 0;
                    TVariation new_member;
                    for(unsigned i = 0; i < align.first.size(); ++i) {
                        if(align.first[i] == master[mpos]) {
                            new_member.push_back(align.second[i]);
                            ++mpos;
                        } else if(master[mpos] == '-') {
                            while(master[mpos] == '-') {
                                new_member.push_back('-');
                                ++mpos;
                            }
                            new_member.push_back(align.second[i]);
                            ++mpos;                                                       
                        } else {  // align.first[i] == '-'
                            for(TVariation& seq : group)
                                seq.insert(seq.begin()+mpos, '-');
                            new_member.push_back(align.second[i]);
                            ++mpos;
                        }
                    }
                    group.push_back(new_member);
                }
                all_variants.erase(it);                
            }
            if(group.empty())
                group.emplace_back(query.begin(), query.end());
            all_groups.push_back(move(group));
            all_variants.pop_front();
        }

        TContigSequenceList new_contigs;
        for(auto& group : all_groups) {
            if(group.size() == 1) {
                new_contigs.push_back(CContigSequence());
                new_contigs.back().InsertNewChunk(group.front().begin(), group.front().end());
                continue;
            } 

            auto NextMismatch = [&](unsigned pos) {
                for( ; pos < group.front().size(); ++pos) {
                    for(auto& seq : group) {
                        if(seq[pos] != group.front()[pos])
                            return pos;
                    }
                }
                return pos;
            };
            
            CContigSequence combined_seq;
            int min_uniq_len = 21;
            for(unsigned mism = NextMismatch(0); mism < group.front().size(); mism = NextMismatch(0)) {
                if(mism > 0) {
                    combined_seq.InsertNewChunk(group.front().begin(), group.front().begin()+mism);
                    for(auto& seq : group)
                        seq.erase(seq.begin(), seq.begin()+mism);
                }
                for(unsigned len = 1; len <= group.front().size(); ) {
                    unsigned next_mism = NextMismatch(len);
                    if(next_mism >= len+min_uniq_len || next_mism == group.front().size()) {
                        map<SeqInterval,set<SeqInterval>> varmap;
                        for(auto& seq : group)
                            varmap[SeqInterval(seq.begin(), seq.begin()+len)].emplace(seq.begin()+len, seq.end());
                        bool all_same = true;
                        for(auto it = varmap.begin(); all_same && ++it != varmap.end(); ) {
                            if(varmap.begin()->second != it->second)
                                all_same = false;
                        }
                        if(all_same) {
                            combined_seq.InsertNewChunk(varmap.begin()->first.begin, varmap.begin()->first.end);
                            for(auto it = varmap.begin(); ++it != varmap.end(); )
                                combined_seq.InsertNewVariant(it->first.begin, it->first.end);
                            for(auto& seq : group)
                                seq.erase(seq.begin(), seq.begin()+len);
                            group.sort();
                            group.erase(unique(group.begin(),group.end()), group.end());
                            break;
                        }
                    }
                    len = next_mism+1;
                }
            }
            if(!group.front().empty())
                combined_seq.InsertNewChunk(group.front().begin(), group.front().end());

            for(auto& chunk : combined_seq) {
                for(auto& seq : chunk)
                    seq.erase(remove(seq.begin(),seq.end(),'-'), seq.end());
            }
            
            new_contigs.push_back(move(combined_seq));
        }

        swap(contigs, new_contigs);
    }    

    typedef list<string> TStrList;
    template<class DBGraph> using TBases = deque<typename DBGraph::Successor>;
    
    template<class DBGraph> struct SContig;
    template<class DBGraph> using TContigList = list<SContig<DBGraph>>;


    template<class DBGraph>
    struct SContig {
        typedef typename DBGraph::Node Node;
        typedef forward_list<Node> TNodeList;
        SContig(DBGraph& graph) :m_graph(graph), m_kmer_len(graph.KmerLen()) {}
        SContig(const CContigSequence& contig, DBGraph& graph) : m_seq(contig), m_graph(graph), m_kmer_len(graph.KmerLen()) { GenerateKmersAndCleanSNPs(); }        

        void GenerateKmersAndCleanSNPs() {
            if(m_seq.RemoveShortUniqIntervals(m_kmer_len))
                RotateCircularToMinKmer();

            int rotation = 0;
            bool extended = false;
            auto& first_chunk =  m_seq.front().front();
            auto& last_chunk = m_seq.back().front();
            if(m_seq.m_circular && (int)(last_chunk.size()+first_chunk.size()) >= m_kmer_len-1) {
                extended = true;
                if((int)first_chunk.size() < m_kmer_len-1) {
                    rotation = m_kmer_len-1-first_chunk.size();
                    first_chunk.insert(first_chunk.begin(), last_chunk.end()-rotation, last_chunk.end());
                    last_chunk.erase(last_chunk.end()-rotation, last_chunk.end());
                }
                last_chunk.insert(last_chunk.end(), first_chunk.begin(), first_chunk.begin()+m_kmer_len-1);
            }

            for(int i = m_seq.size()-1; i >= 0; ) {
                if(i == (int)m_seq.size()-1) {
                    if((int)m_seq.ChunkLenMax(i) >= m_kmer_len) {  // last chunk size >= kmer_len
                        CReadHolder rh(false);
                        rh.PushBack(m_seq.back().front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik) {
                            Node node = m_graph.GetNode(*ik);
                            if(node.isValid() && !m_graph.SetVisited(node))
                                m_graph.SetMultContig(node);
                        }
                    }
                    
                    --i;
                } else { // all uniq chunks >= kmer_len-1
                    if((int)m_seq.ChunkLenMax(i-1) >= m_kmer_len) {
                        CReadHolder rh(false);
                        rh.PushBack(m_seq[i-1].front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len); ik != rh.kend(); ++ik) {
                            Node node = m_graph.GetNode(*ik);
                            if(node.isValid() && !m_graph.SetVisited(node))
                                m_graph.SetMultContig(node);
                        }
                    }

                    unordered_set<Node, typename Node::Hash> kmers;
                    list<deque<Node>> failed_nodes;
                    list<TLocalVariants::iterator> failed_variants;
                    for(auto prev = m_seq[i].before_begin(); ;++prev) {
                        auto current = prev;
                        if(++current == m_seq[i].end())
                            break;

                        auto& variant = *current;
                        int left = min(m_kmer_len-1, (int)m_seq.ChunkLenMax(i-1));
                        TVariation var_seq(m_seq[i-1].front().end()-left, m_seq[i-1].front().end());
                        var_seq.insert(var_seq.end(), variant.begin(), variant.end());
                        int right = min(m_kmer_len-1, (int)m_seq.ChunkLenMax(i+1));
                        var_seq.insert(var_seq.end(), m_seq[i+1].front().begin(), m_seq[i+1].front().begin()+right);
                        if(i == 1 && m_seq.m_circular && !extended) { // can happen only if there is one long varianle chunk and short ends
                            if(m_seq.size() != 3)
                                throw runtime_error("Error in circular extension");
                            var_seq.insert(var_seq.end(), var_seq.begin(), var_seq.begin()+m_kmer_len-1);
                        }
                        CReadHolder rh(false);
                        rh.PushBack(var_seq);
                        deque<Node> var_nodes;
                        bool failed = false;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik)  {
                            var_nodes.emplace_back(m_graph.GetNode(*ik));
                            if(!var_nodes.back().isValid())
                                failed = true;
                        }
                        if(failed) {
                            failed_nodes.push_back(move(var_nodes));
                            failed_variants.push_front(prev);           // reverse order for deleting
                        } else {
                            kmers.insert(var_nodes.begin(), var_nodes.end()); 
                        }
                    }

                    if((int)failed_variants.size() == m_seq.VariantsNumber(i)) { // all failed
                        for(auto& nodes : failed_nodes)
                            kmers.insert(nodes.begin(), nodes.end());
                    } else { // some are good
                        for(auto prev : failed_variants)
                            m_seq[i].erase_after(prev);
                        if(m_seq.UniqueChunk(i)) { // only one left
                            m_seq[i-1].front().insert(m_seq[i-1].front().end(), m_seq[i].front().begin(), m_seq[i].front().end());
                            m_seq[i-1].front().insert(m_seq[i-1].front().end(), m_seq[i+1].front().begin(), m_seq[i+1].front().end());
                            m_seq.erase(m_seq.begin()+i, m_seq.begin()+i+2);
                        }
                    }

                    for(auto& node : kmers) {
                        if(node.isValid() && !m_graph.SetVisited(node))
                            m_graph.SetMultContig(node);                        
                    }
                     
                    i -= 2;   
                }
            }
            if(m_seq.m_circular && extended) {
                last_chunk.erase(last_chunk.end()-m_kmer_len+1, last_chunk.end());
                if(rotation > 0) {
                    last_chunk.insert(last_chunk.end(), first_chunk.begin(), first_chunk.begin()+rotation);
                    first_chunk.erase(first_chunk.begin(), first_chunk.begin()+rotation);
                }
            }
            m_seq.ContractVariableIntervals();
            m_seq.IncludeRepeatsInVariableIntervals();
            if(m_seq.RemoveShortUniqIntervals(m_kmer_len))
                RotateCircularToMinKmer();
            m_seq.StabilizeVariantsOrder();
        }

        SContig(const SContig& to_left, const SContig& to_right, const Node& initial_node, const Node& lnode, const Node& rnode, DBGraph& graph) :  
        m_next_left(lnode), m_next_right(rnode), m_graph(graph),  m_kmer_len(graph.KmerLen()) {                                                                                                                                                   
//          initial_node - the starting kmer
//          to_left - left extension of the starting kmer
//          to_right - right extension of the starting kmer
//          lnode - left denied node
//          rnode - right denied node
//          graph - de Bruijn graph

                                                                       
            // take parts of the assembled sequence and put them together in SContig

            if(!to_left.m_seq.empty()) {
                m_seq = to_left.m_seq;
                ReverseComplement();
            }

            // could be changed by ReverseComplement
            m_next_left = lnode;
            m_next_right = rnode;

            string ikmer = graph.GetNodeSeq(initial_node);
            if(m_seq.empty() || m_seq.VariableChunk(m_seq.size()-1))  // empty or variant
                m_seq.InsertNewChunk(ikmer.begin(), ikmer.end());
            else
                m_seq.ExtendTopVariant(ikmer.begin(), ikmer.end());

            if(!to_right.m_seq.empty()) {
                if(to_right.m_seq.UniqueChunk(0)) {
                    m_seq.ExtendTopVariant(to_right.m_seq.front().front().begin(), to_right.m_seq.front().front().end());
                    m_seq.insert(m_seq.end(), to_right.m_seq.begin()+1, to_right.m_seq.end());
                } else {
                    m_seq.insert(m_seq.end(), to_right.m_seq.begin(), to_right.m_seq.end());
                }
            }
            m_seq.StabilizeVariantsOrder();
                
            m_left_extend = m_right_extend = LenMax();
        }
        SContig(SContig* link, int shift, const Node& takeoff_node, const SContig& extension, const Node& rnode, DBGraph& graph) :
            m_next_left(takeoff_node), m_next_right(rnode), m_left_link(link), m_left_shift(shift),  m_graph(graph), m_kmer_len(graph.KmerLen()) {

            string kmer = graph.GetNodeSeq(takeoff_node);
            m_seq.InsertNewChunk(kmer.begin()+1, kmer.end()); // don't include first base
            if(!extension.m_seq.empty()) {
                if(extension.m_seq.UniqueChunk(0)) {
                    m_seq.ExtendTopVariant(extension.m_seq.front().front().begin(), extension.m_seq.front().front().end());
                    m_seq.insert(m_seq.end(), extension.m_seq.begin()+1, extension.m_seq.end());
                } else {
                    m_seq.insert(m_seq.end(), extension.m_seq.begin(), extension.m_seq.end());
                }
            }
            m_seq.StabilizeVariantsOrder();

            m_left_extend = m_right_extend = LenMax();
        }

        Node FrontKmer() const {
            if(m_seq.VariableChunk(0) || (int)m_seq.ChunkLenMax(0) < m_kmer_len)
                return Node();

            string kmer_seq(m_seq.front().front().begin(), m_seq.front().front().begin()+m_kmer_len);
            TKmer kmer(kmer_seq); // front must be unambiguous
            return m_graph.GetNode(kmer);
        }
        Node BackKmer() const { 
            int last = m_seq.size()-1;
            if(m_seq.VariableChunk(last) || (int)m_seq.ChunkLenMax(last) < m_kmer_len)
                return Node();

            string kmer_seq(m_seq.back().front().end()-m_kmer_len, m_seq.back().front().end());
            TKmer kmer(kmer_seq);
            return m_graph.GetNode(kmer);
        }

        // don't 'own' any kmers
        bool EmptyLinker() const { return ((int)max(m_seq.ChunkLenMax(0), m_seq.ChunkLenMax(m_seq.size()-1)) < m_kmer_len && m_seq.size() <= 3); }

        bool RightSNP() const { return (m_seq.size() >= 3 && m_seq.UniqueChunk(m_seq.size()-1) && (int)m_seq.ChunkLenMax(m_seq.size()-1) < m_kmer_len); }
        bool LeftSNP() const { return (m_seq.size() >= 3 && m_seq.UniqueChunk(0) &&  (int)m_seq.ChunkLenMax(0) < m_kmer_len); }

        Node RightConnectingNode() const {
            int last_index = m_seq.size()-1;
            if((int)m_seq.ChunkLenMax(last_index) >= m_kmer_len) {   // normal end
                return BackKmer(); 
            } else if(m_seq.size() >= 3) {                           // snp
                if((int)m_seq.ChunkLenMax(last_index-2) >= m_kmer_len) {
                    string kmer_seq(m_seq[last_index-2].front().end()-m_kmer_len, m_seq[last_index-2].front().end());
                    TKmer kmer(kmer_seq);
                    return m_graph.GetNode(kmer);
                }
            }

            return m_next_left;                        // empty linker
        }
        Node LeftConnectingNode() const {
            if((int)m_seq.ChunkLenMax(0) >= m_kmer_len) {  // normal end
                return FrontKmer();
            } else if(m_seq.size() >= 3) {                 // snp
                if((int)m_seq.ChunkLenMax(2) >= m_kmer_len) {
                    string kmer_seq(m_seq[2].front().begin(), m_seq[2].front().begin()+m_kmer_len);
                    TKmer kmer(kmer_seq); // front must be unambiguous
                    return m_graph.GetNode(kmer);
                }
            }

            return m_next_right;                        // empty linker
        }
        
        void ReverseComplement() {  
            m_seq.ReverseComplement();
            swap(m_next_left, m_next_right);
            m_next_left = DBGraph::ReverseComplement(m_next_left);
            m_next_right = DBGraph::ReverseComplement(m_next_right);
            swap(m_left_link, m_right_link);
            swap(m_left_shift, m_right_shift);
            swap(m_left_extend, m_right_extend);
        }
        void AddToRight(const SContig& other) {
            m_seq.m_circular = false;

            m_next_right = other.m_next_right;
            m_right_link = other.m_right_link;
            m_right_shift = other.m_right_shift; 
            if(EmptyLinker() && other.EmptyLinker())
                return;            
 
            auto& last_chunk = m_seq.back().front();
            int last_chunk_len = last_chunk.size();
            int overlap = m_kmer_len-1;
            auto first_other_chunk_it = other.m_seq.begin(); 
            if(RightSNP() && other.LeftSNP()) {    // skip snp chunk
                overlap = last_chunk_len+other.m_seq.ChunkLenMax(1)+first_other_chunk_it->front().size();
                first_other_chunk_it += 2;
            }

            if(other.m_right_extend < (int)other.LenMax()) {
                m_right_extend = other.m_right_extend;
            } else {
                m_right_extend += other.m_right_extend-overlap;
                if(m_left_extend == (int)LenMax())
                    m_left_extend = m_right_extend;
            }

            auto& first_other_chunk = first_other_chunk_it->front();
            last_chunk.insert(last_chunk.end(), first_other_chunk.begin()+min(m_kmer_len-1,last_chunk_len), first_other_chunk.end());  // combine overlapping chunks
            m_seq.insert(m_seq.end(), first_other_chunk_it+1, other.m_seq.end());                                                       // insert remaining chunks
        }
        void AddToLeft(const SContig& other) {
            m_seq.m_circular = false;

            m_next_left = other.m_next_left;
            m_left_link = other.m_left_link;
            m_left_shift = other.m_left_shift;
            if(EmptyLinker() && other.EmptyLinker())
                return;            
 
            auto& first_chunk = m_seq.front().front();
            int first_chunk_len = first_chunk.size();
            int overlap = m_kmer_len-1;
            auto last_other_chunk_it = other.m_seq.end()-1; 
            if(LeftSNP() && other.RightSNP()) {    // skip snp chunk
                overlap = first_chunk_len+other.m_seq.ChunkLenMax(other.m_seq.size()-2)+last_other_chunk_it->front().size();
                last_other_chunk_it -= 2;
            }
                
            if(other.m_left_extend < (int)other.LenMax()) {
                m_left_extend = other.m_left_extend;
            } else {
                m_left_extend += other.m_left_extend-overlap; 
                if(m_right_extend == (int)LenMax())
                    m_right_extend = m_left_extend;
            }

            auto& last_other_chunk = last_other_chunk_it->front();
            first_chunk.insert(first_chunk.begin(),last_other_chunk.begin(), last_other_chunk.end()-min(m_kmer_len-1,first_chunk_len));  // combine overlapping chunks
            m_seq.insert(m_seq.begin(), other.m_seq.begin(), last_other_chunk_it);                                     // insert remaining chunks
        }
        
        void ClipRight(int clip) {
            if(clip <= 0)
                return;

            m_seq.m_circular = false;
            m_next_right = Node();
            m_right_link = nullptr;
            m_right_shift = 0;            

            while(!m_seq.empty() && (m_seq.VariableChunk(m_seq.size()-1) || (int)m_seq.ChunkLenMax(m_seq.size()-1) <= clip)) {
                int chunk_len = m_seq.ChunkLenMax(m_seq.size()-1);
                clip -= chunk_len;
                m_right_extend = max(0, m_right_extend-chunk_len);
                m_seq.pop_back();
            }
            if(clip > 0 && !m_seq.empty()) {
                m_right_extend = max(0, m_right_extend-clip);
                m_seq.back().front().erase(m_seq.back().front().end()-clip, m_seq.back().front().end());
            }

            if((int)LenMin() < m_kmer_len-1)
                m_seq.clear();
        }
        void ClipLeft(int clip) {
            if(clip <= 0)
                return;

            m_seq.m_circular = false;
            m_next_left = Node();
            m_left_link = nullptr;
            m_left_shift = 0; 

            while(!m_seq.empty() && (m_seq.VariableChunk(0) || (int)m_seq.ChunkLenMax(0) <= clip)) {
                int chunk_len = m_seq.ChunkLenMax(0);
                clip -= chunk_len;
                m_left_extend = max(0, m_left_extend-chunk_len);
                m_seq.erase(m_seq.begin());
            }
            if(clip > 0 && !m_seq.empty()) {
                m_left_extend = max(0, m_left_extend-clip);
                m_seq.front().front().erase(m_seq.front().front().begin(), m_seq.front().front().begin()+clip);
            }         

            if((int)LenMin() < m_kmer_len-1)
                m_seq.clear();
        }        

        size_t LenMax() const { return m_seq.LenMax(); }
        size_t LenMin() const { return m_seq.LenMin(); }

        tuple<int, int, int> MinKmerPosition() const {  //chunk, position in chunk, strand
            int kmer_len = min(21, m_kmer_len); // duplicated in RotateCircularToMinKmer()
            typedef LargeInt<1> large_t;
            unordered_map<large_t, tuple<int, int, int>, SKmerHash> kmers; // [kmer], chunk, position in chunk, strand/notvalid

            for(int i = m_seq.size()-1; i >= 0; i -= 2) {
                deque<forward_list<LargeInt<1>>> chunk_kmers;
                if(i == (int)m_seq.size()-1) {
                    if((int)m_seq.ChunkLenMax(i) >= kmer_len) { // last chunk could be short
                        chunk_kmers.resize(m_seq.ChunkLenMax(i)-kmer_len+1);
                        CReadHolder rh(false);
                        rh.PushBack(m_seq.back().front());
                        int pos = chunk_kmers.size();
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) // iteration from last kmer to first     
                            chunk_kmers[--pos].push_front(get<large_t>(TKmer::Type(*ik)));
                    }
                } else { // all uniq chunks in the middle >= kmer_len-1; first/last could be short
                    chunk_kmers.resize(m_seq.ChunkLenMax(i)+m_seq.ChunkLenMax(i+1));
                    if((int)m_seq.ChunkLenMax(i) >= kmer_len) {
                        TVariation seq(m_seq[i].front().begin(), m_seq[i].front().end());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = seq.size()-kmer_len+1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first 
                            chunk_kmers[--pos].push_front(get<large_t>(TKmer::Type(*ik)));
                    }
                    for(auto& variant : m_seq[i+1]) {
                        TVariation seq;
                        if((int)m_seq.ChunkLenMax(i) >= kmer_len-1)
                            seq.insert(seq.end(), m_seq[i].front().end()-kmer_len+1, m_seq[i].front().end());
                        else
                            seq.insert(seq.end(), m_seq[i].front().begin(), m_seq[i].front().end());
                        seq.insert(seq.end(), variant.begin(), variant.end());
                        if((int)m_seq.ChunkLenMax(i+2) >= kmer_len-1)
                            seq.insert(seq.end(), m_seq[i+2].front().begin(), m_seq[i+2].front().begin()+kmer_len-1);
                        else
                            seq.insert(seq.end(), m_seq[i+2].front().begin(), m_seq[i+2].front().end());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = seq.size()-kmer_len+1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first 
                            chunk_kmers[--pos].push_front(get<large_t>(TKmer::Type(*ik)));
                    }
                }
                
                for(unsigned pos = 0; pos < chunk_kmers.size(); ++pos) {
                    int k = pos;
                    int chunk = i;
                    if(pos >= m_seq.ChunkLenMax(i)) {
                        k = pos-m_seq.ChunkLenMax(i);
                        chunk = i+1;
                    }
                    for(auto& kmer : chunk_kmers[pos]) {
                        int strand = 1;
                        large_t* min_kmerp = &kmer;
                        large_t rkmer = revcomp(kmer, kmer_len);
                        if(rkmer < kmer) {
                            strand = -1;
                            min_kmerp = &rkmer;
                        }
                        auto rslt = kmers.insert(make_pair(*min_kmerp, make_tuple(chunk, k, strand)));
                        if(!rslt.second)
                            get<2>(rslt.first->second) = 0;                        
                    }
                }                
            }

            tuple<int, int, int> rslt(0, 0, 0);
            large_t min_kmer;            
            for(auto& elem : kmers) {
                if(get<2>(elem.second)) { // not a repeat                    
                    if(!get<2>(rslt) || elem.first < min_kmer) {
                        min_kmer = elem.first;
                        rslt = elem.second;
                    }
                }
            }
            
            return rslt;
        }

        // stabilize contig orientation using minimal kmer in the contig
        void SelectMinDirection() {
            int strand = get<2>(MinKmerPosition());
            if(strand < 0)
               ReverseComplement();  
        }

        // finds stable origin for circular contigs by placing minimal kmer at the beginning of the sequence
        void RotateCircularToMinKmer() { // assumes that the next extension of sequence would give the first kmer (m_next_right == m_kmers.front())

            int kmer_len = min(21, m_kmer_len);
            m_seq.back().front().erase(m_seq.back().front().end()-(m_kmer_len-kmer_len), m_seq.back().front().end()); // clip extra portion of overlap
            auto rslt = MinKmerPosition();
            if(get<2>(rslt) == 0)
                return;

            m_seq.back().front().erase(m_seq.back().front().end()-kmer_len+1, m_seq.back().front().end());  // clip remaining overlap
            size_t first_chunk = get<0>(rslt);
            size_t first_base = get<1>(rslt);

            if(get<2>(rslt) < 0) {
                first_base += min(21, m_kmer_len);
                while(first_base >= m_seq.ChunkLenMax(first_chunk)) {
                    first_base -= m_seq.ChunkLenMax(first_chunk);
                    first_chunk = (first_chunk+1)%m_seq.size();
                }
                if(m_seq.VariableChunk(first_chunk)) {               // ambiguous interval - we don't want to cut it
                    ++first_chunk;  // variable chunk cant be last
                    first_base = 1;
                } else if(first_chunk > 0 && first_base == 0) {       // we want some uniq intervals on both ends
                    first_base = 1;
                }
            } else {
                if(m_seq.VariableChunk(first_chunk)) {                // ambiguous interval - we don't want to cut it
                    --first_chunk;  // variable chunk cant be first
                    first_base = m_seq.ChunkLenMax(first_chunk)-1;     // leave one base
                } else if(first_chunk > 0 && first_base == 0) {       // we want some uniq intervals on both ends
                    first_chunk -= 2;
                    first_base = m_seq.ChunkLenMax(first_chunk)-1;     // leave one base
                }
            }

            if(m_seq.size() == 1) {
                rotate(m_seq.front().front().begin(), m_seq.front().front().begin()+first_base, m_seq.front().front().end()); 
            } else {
                if(first_chunk > 0) {
                    auto& last_seq = m_seq.back().front();
                    last_seq.insert(last_seq.end(), m_seq.front().front().begin(), m_seq.front().front().end());
                    m_seq.erase(m_seq.begin());
                    rotate(m_seq.begin(), m_seq.begin()+first_chunk-1, m_seq.end());
                }
                if(first_base > 0) {
                    if(m_seq.VariableChunk(m_seq.size()-1)) {
                        m_seq.InsertNewChunk();
                        m_seq.InsertNewVariant();
                    }
                    auto& last_seq = m_seq.back().front();
                    last_seq.insert(last_seq.end(), m_seq.front().front().begin(), m_seq.front().front().begin()+first_base);
                    m_seq.front().front().erase(m_seq.front().front().begin(), m_seq.front().front().begin()+first_base);
                }
            }
            
            //clean edges   
            m_next_left = Node();
            m_next_right = Node();
            m_left_link = nullptr;
            m_left_shift = 0;
            m_right_link = nullptr;
            m_right_shift = 0;
            m_left_extend = 0;   // prevents any further clipping    
            m_right_extend = 0;  // prevents any further clipping    
            m_seq.m_circular = true;
        }

        bool operator<(const SContig& other) const { return m_seq < other.m_seq; }

        // connects fragments created in different threads and combines doubled 'empty' linkers
        static TContigList<DBGraph> ConnectFragments(vector<TContigList<DBGraph>>& fragments, const DBGraph& graph) {

            int total = 0;
            size_t len = 0;
            for(auto& ns : fragments) {
                for(auto& seq : ns) {
                    /*
                    cerr << "Fragment: " << seq.m_left_link << " " << seq.m_right_link << endl;
                    cerr << "Lnode: ";
                    if(seq.m_next_left.isValid())
                        cerr << graph.GetNodeSeq(seq.m_next_left);
                    cerr << endl;
                    cerr << "Rnode: ";
                    if(seq.m_next_right.isValid())
                        cerr << graph.GetNodeSeq(seq.m_next_right);
                    cerr << endl;

                    for(auto& chunk : seq.m_seq) {
                        cerr << "Chunk:" << endl;
                        for(auto& var : chunk) {
                            cerr << "Variant    : ";
                            for(char c : var) cerr << c;
                            cerr << endl;
                        }
                    }                    
                    */
                    ++total;
                    len += seq.LenMax()-graph.KmerLen()+1;
                }
            }
            
            cerr << "Fragments before: " << total << " " <<  len << endl;

            TContigList<DBGraph> connected;        
            unordered_map<Node, typename TContigList<DBGraph>::iterator, typename Node::Hash> denied_left_nodes;
            unordered_map<Node, typename TContigList<DBGraph>::iterator, typename Node::Hash> denied_right_nodes;
            for(auto& ns : fragments) {
                for(auto iloop = ns.begin(); iloop != ns.end(); ) {
                    auto ic = iloop++;
                    connected.splice(connected.begin(), ns, ic);
                    SContig& contig = *connected.begin();                
                    if(contig.m_next_left > contig.m_next_right)  // need this to pair two identical empty links
                        contig.ReverseComplement();                

                    if(contig.m_next_left.isValid()) {
                        auto rslt = denied_left_nodes.insert(make_pair(contig.m_next_left, connected.begin()));
                        if(!rslt.second) {
                            typename TContigList<DBGraph>::iterator other = rslt.first->second;

                            if(contig.m_left_link && other->m_right_link) { // other started from end of contig and went all the way to another contig
                                other->m_left_link = contig.m_left_link;    // add left link to other
                                connected.pop_front();
                                continue;
                            } else if(other->m_left_link && contig.m_right_link) { // contig started from end of contig and went all the way to another contig
                                contig.m_left_link = other->m_left_link;           // add left link to contig
                                rslt.first->second = connected.begin();
                                if(other->m_next_right.isValid())
                                    denied_right_nodes.erase(other->m_next_right);
                                connected.erase(other);
                            
                            /*
                            if(contig.EmptyLinker() && contig.m_left_link && !other->m_left_link && contig.m_next_right == other->LeftConnectingNode()) {
                                other->AddToLeft(contig); // add left link to other
                                connected.pop_front();
                                continue;
                            }else if(other->EmptyLinker() && other->m_left_link && !contig.m_left_link && other->m_next_right == contig.LeftConnectingNode()) {
                                contig.AddToLeft(*other); // add left link to contig    
                                rslt.first->second = connected.begin();
                                denied_right_nodes.erase(other->m_next_right);
                                connected.erase(other);
                            */

                            } else {
                                cerr << "Unexpected left fork: " << graph.GetNodeSeq(contig.m_next_left) << endl;

                                cerr << "Contig: " << contig.m_left_link << " " << contig.m_right_link << " ";
                                if(contig.m_next_left.isValid())
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_left) << " ";
                                else
                                    cerr << "LC notvalid ";
                                if(contig.m_next_right.isValid())
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_right) << " ";
                                else
                                    cerr << "RC notvalid ";
                                cerr << endl;
                                for(auto& chunk : contig.m_seq) {
                                    cerr << "Chunk:" << endl;
                                    for(auto& var : chunk) {
                                        cerr << "Variant: ";
                                        for(char c : var) cerr << c;
                                        cerr << endl;
                                    }
                                }

                                auto& other = *rslt.first->second;

                                cerr << "Other: " << other.m_left_link << " " << other.m_right_link << " ";
                                if(other.m_next_left.isValid())
                                    cerr << other.m_graph.GetNodeSeq(other.m_next_left) << " ";
                                else
                                    cerr << "LC notvalid ";
                                if(other.m_next_right.isValid())
                                    cerr << other.m_graph.GetNodeSeq(other.m_next_right) << " ";
                                else
                                    cerr << "RC notvalid ";
                                cerr << endl;
                                for(auto& chunk : other.m_seq) {
                                    cerr << "Chunk:" << endl;
                                    for(auto& var : chunk) {
                                        cerr << "Variant: ";
                                        for(char c : var) cerr << c;
                                        cerr << endl;
                                    }
                                }
                                
                            }
                        }
                    }
                    if(contig.m_next_right.isValid()) {
                        auto rslt = denied_right_nodes.insert(make_pair(contig.m_next_right, connected.begin()));
                        if(!rslt.second) {
                            typename TContigList<DBGraph>::iterator other = rslt.first->second;

                            if(contig.m_right_link && other->m_left_link) { // other started from end of contig and went all the way to another contig
                                other->m_right_link = contig.m_right_link;  // add right link to other
                                denied_left_nodes.erase(contig.m_next_left);
                                connected.pop_front();
                            } else if(other->m_right_link && contig.m_left_link) { // contig started from end of contig and went all the way to another contig
                                contig.m_right_link = other->m_right_link;         // add right link to contig
                                rslt.first->second = connected.begin();
                                if(other->m_next_left.isValid())
                                    denied_left_nodes.erase(other->m_next_left);
                                connected.erase(other);

                                /*
                            if(contig.EmptyLinker() && contig.m_right_link && !other->m_right_link && contig.m_next_left == other->RightConnectingNode()) {
                                other->AddToRight(contig); // add right link to other
                                denied_left_nodes.erase(contig.m_next_left);
                                connected.pop_front();
                            } else if (other->EmptyLinker() && other->m_right_link &&  !contig.m_right_link && other->m_next_left == contig.RightConnectingNode()) {
                                contig.AddToRight(*other); // add right link to contig
                                rslt.first->second = connected.begin();
                                denied_left_nodes.erase(other->m_next_left);
                                connected.erase(other);
                                */

                            } else {
                                cerr << "Unexpected right fork: " << graph.GetNodeSeq(contig.m_next_right) << endl;                                

                                cerr << "Contig: " << contig.m_left_link << " " << contig.m_right_link << " ";
                                if(contig.m_next_left.isValid())
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_left) << " ";
                                else
                                    cerr << "LC notvalid ";
                                if(contig.m_next_right.isValid())
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_right) << " ";
                                else
                                    cerr << "RC notvalid ";
                                cerr << endl;
                                for(auto& chunk : contig.m_seq) {
                                    cerr << "Chunk:" << endl;
                                    for(auto& var : chunk) {
                                        cerr << "Variant: ";
                                        for(char c : var) cerr << c;
                                        cerr << endl;
                                    }
                                }

                                auto& other = *rslt.first->second;

                                cerr << "Other: " << other.m_left_link << " " << other.m_right_link << " ";
                                if(other.m_next_left.isValid())
                                    cerr << other.m_graph.GetNodeSeq(other.m_next_left) << " ";
                                else
                                    cerr << "LC notvalid ";
                                if(other.m_next_right.isValid())
                                    cerr << other.m_graph.GetNodeSeq(other.m_next_right) << " ";
                                else
                                    cerr << "RC notvalid ";
                                cerr << endl;
                                for(auto& chunk : other.m_seq) {
                                    cerr << "Chunk:" << endl;
                                    for(auto& var : chunk) {
                                        cerr << "Variant: ";
                                        for(char c : var) cerr << c;
                                        cerr << endl;
                                    }
                                }

                            }
                        }
                    }
                }
            }

            for(SContig& contig : connected) {
                if(contig.EmptyLinker())
                    continue; 

                if(contig.m_next_right.isValid())
                    denied_right_nodes.erase(contig.m_next_right);
                if(contig.m_next_left.isValid())
                    denied_left_nodes.erase(contig.m_next_left);
                bool keep_doing = true;
                while(keep_doing) {
                    keep_doing = false;
                    if(contig.m_next_right.isValid()) {
                        Node rnode = contig.RightConnectingNode();
                        auto rslt = denied_left_nodes.find(rnode);
                        if(rslt != denied_left_nodes.end()) {
                            keep_doing = true;
                            SContig& rcontig = *rslt->second;
                            if(rcontig.m_next_right.isValid()) 
                                denied_right_nodes.erase(rcontig.m_next_right);
                            contig.AddToRight(rcontig);
                            connected.erase(rslt->second);
                            denied_left_nodes.erase(rslt);
                        } else if((rslt = denied_right_nodes.find(DBGraph::ReverseComplement(rnode))) != denied_right_nodes.end()) {
                            keep_doing = true;
                            SContig& rcontig = *rslt->second;
                            if(rcontig.m_next_left.isValid()) 
                                denied_left_nodes.erase(rcontig.m_next_left);
                            rcontig.ReverseComplement();
                            contig.AddToRight(rcontig);
                            connected.erase(rslt->second);
                            denied_right_nodes.erase(rslt);
                        }
                    }               
                    if(contig.m_next_left.isValid()) {
                        Node lnode = contig.LeftConnectingNode();
                        auto rslt = denied_right_nodes.find(lnode);
                        if(rslt != denied_right_nodes.end()) {
                            keep_doing = true;
                            SContig& lcontig = *rslt->second;
                            if(lcontig.m_next_left.isValid()) 
                                denied_left_nodes.erase(lcontig.m_next_left);
                            contig.AddToLeft(lcontig);
                            connected.erase(rslt->second);
                            denied_right_nodes.erase(rslt);
                        } else if((rslt = denied_left_nodes.find(DBGraph::ReverseComplement(lnode))) != denied_left_nodes.end()) {
                            keep_doing = true;
                            SContig& lcontig = *rslt->second;
                            if(lcontig.m_next_right.isValid()) 
                                denied_right_nodes.erase(lcontig.m_next_right);
                            lcontig.ReverseComplement();
                            contig.AddToLeft(lcontig);
                            connected.erase(rslt->second);
                            denied_left_nodes.erase(rslt);
                        }
                    }                                
                }
            
                if(contig.m_next_right == contig.LeftConnectingNode() && (int)contig.LenMax() >= 2*graph.KmerLen()-1)  // circular and not very short  
                    contig.RotateCircularToMinKmer();                     
            }
            
            total = 0;
            len = 0;
            for(auto& seq : connected) {
                ++total;
                len += seq.LenMax()-graph.KmerLen()+1;
            }            
            
            cerr << "Fragments after: " << total << " " <<  len << endl;

            return connected;                
        }

        // connects and extends contigs from previous iteration using a longer kmer
        // scontigs - previous contigs
        // extensions - connectors and extenders produced by longer kmer
        static void ConnectAndExtendContigs(TContigList<DBGraph>& scontigs, TContigList<DBGraph>& extensions, int ncores) {
            if(scontigs.empty())
                return;

            int kmer_len = scontigs.front().m_kmer_len;
            int connectors = 0;
            int extenders = 0;
            //assign links to main contigs
            for(auto& ex : extensions) {
                if(ex.m_left_link && ex.m_right_link)
                    ++connectors;
                else
                    ++extenders;

                if(ex.m_left_link) {
                    auto& contig = *ex.m_left_link;
                    if(contig.RightConnectingNode() == ex.m_next_left) {
                        if(contig.m_right_link)
                            throw runtime_error("Multiple connection of contigs");        
                        contig.m_right_link = &ex;
                    } else if(contig.LeftConnectingNode() == DBGraph::ReverseComplement(ex.m_next_left)) {
                        if(contig.m_left_link)
                            throw runtime_error("Multiple connection of contigs");        
                        contig.m_left_link = &ex;
                    } else {

                        cerr << "Corrupted connection of contigs L" << endl;
                        cerr << "Contig: ";
                        if(contig.LeftConnectingNode().isValid())
                            cerr << contig.m_graph.GetNodeSeq(contig.LeftConnectingNode()) << " ";
                        else
                            cerr << "LC notvalid ";
                        if(contig.RightConnectingNode().isValid())
                            cerr << contig.m_graph.GetNodeSeq(contig.RightConnectingNode()) << " ";
                        else
                            cerr << "RC notvalid ";
                        cerr << endl;
                        for(auto& chunk : contig.m_seq) {
                            cerr << "Chunk:" << endl;
                            for(auto& var : chunk) {
                                cerr << "Variant: ";
                                for(char c : var) cerr << c;
                                cerr << endl;
                            }
                        }

                        cerr << "Extension: ";
                        if(ex.m_next_left.isValid())
                            cerr << contig.m_graph.GetNodeSeq(ex.m_next_left);
                        else
                            cerr << "Notvalid";
                        cerr << endl;
                        for(auto& chunk : ex.m_seq) {
                            cerr << "Chunk:" << endl;
                            for(auto& var : chunk) {
                                cerr << "Variant: ";
                                for(char c : var) cerr << c;
                                cerr << endl;
                            }
                        }

                        if(ex.m_next_right.isValid())
                            cerr << "NR: " << contig.m_graph.GetNodeSeq(ex.m_next_right) << endl;
                        if(ex.m_right_link) {
                            auto& contig = *ex.m_right_link;
                            cerr << "RContig: ";
                            if(contig.LeftConnectingNode().isValid())
                                cerr << contig.m_graph.GetNodeSeq(contig.LeftConnectingNode()) << " ";
                            else
                                cerr << "LC notvalid ";
                            if(contig.RightConnectingNode().isValid())
                                cerr << contig.m_graph.GetNodeSeq(contig.RightConnectingNode()) << " ";
                            else
                                cerr << "RC notvalid ";
                            cerr << endl;
                            for(auto& chunk : contig.m_seq) {
                                cerr << "Chunk:" << endl;
                                for(auto& var : chunk) {
                                    cerr << "Variant: ";
                                    for(char c : var) cerr << c;
                                    cerr << endl;
                                }
                            }
                        }

                        //                        throw runtime_error("Corrupted connection of contigs L");
                    }
                }
                if(ex.m_right_link) {
                    auto& contig = *ex.m_right_link;
                    if(contig.LeftConnectingNode() == ex.m_next_right) {
                        if(contig.m_left_link)
                            throw runtime_error("Multiple connection of contigs");        
                        contig.m_left_link = &ex;
                    } else if(contig.RightConnectingNode() == DBGraph::ReverseComplement(ex.m_next_right)) {
                        if(contig.m_right_link)
                            throw runtime_error("Multiple connection of contigs");        
                        contig.m_right_link = &ex;
                    } else {
                        cerr << "Corrupted connection of contigs R" << endl;
                        cerr << "Contig: ";
                        if(contig.LeftConnectingNode().isValid())
                            cerr << contig.m_graph.GetNodeSeq(contig.LeftConnectingNode()) << " ";
                        else
                            cerr << "LC notvalid ";
                        if(contig.RightConnectingNode().isValid())
                            cerr << contig.m_graph.GetNodeSeq(contig.RightConnectingNode()) << " ";
                        else
                            cerr << "RC notvalid ";
                        cerr << endl;
                        for(auto& chunk : contig.m_seq) {
                            cerr << "Chunk:" << endl;
                            for(auto& var : chunk) {
                                cerr << "Variant: ";
                                for(char c : var) cerr << c;
                                cerr << endl;
                            }
                        }

                        cerr << "Extension: ";
                        if(ex.m_next_right.isValid())
                            cerr << contig.m_graph.GetNodeSeq(ex.m_next_right);
                        else
                            cerr << "Notvalid";
                        cerr << endl;
                        for(auto& chunk : ex.m_seq) {
                            cerr << "Chunk:" << endl;
                            for(auto& var : chunk) {
                                cerr << "Variant: ";
                                for(char c : var) cerr << c;
                                cerr << endl;
                            }
                        }

                        if(ex.m_next_left.isValid())
                            cerr << "NL: " << contig.m_graph.GetNodeSeq(ex.m_next_left) << endl;
                        if(ex.m_left_link) {
                            auto& contig = *ex.m_left_link;
                            cerr << "LContig: ";
                            if(contig.LeftConnectingNode().isValid())
                                cerr << contig.m_graph.GetNodeSeq(contig.LeftConnectingNode()) << " ";
                            else
                                cerr << "LC notvalid ";
                            if(contig.RightConnectingNode().isValid())
                                cerr << contig.m_graph.GetNodeSeq(contig.RightConnectingNode()) << " ";
                            else
                                cerr << "RC notvalid ";
                            cerr << endl;
                            for(auto& chunk : contig.m_seq) {
                                cerr << "Chunk:" << endl;
                                for(auto& var : chunk) {
                                    cerr << "Variant: ";
                                    for(char c : var) cerr << c;
                                    cerr << endl;
                                }
                            }
                        }


                        //                        throw runtime_error("Corrupted connection of contigs R");
                    }
                }            
            }
            cerr << "Connectors: " << connectors << " Extenders: " << extenders << endl;


            for(auto& contig : scontigs)
                contig.m_is_taken = 0;

            //select starting points for chains
            for(auto& contig : scontigs) {
                if(contig.m_is_taken)
                    continue;
                if(contig.m_left_link == nullptr && contig.m_right_link == nullptr) {
                    contig.m_is_taken = 1;
                    continue;
                }

                //mark as taken all chain members except the starting point
                auto parent = &contig;
                bool circular = false;
                for(auto child = parent->m_right_link; child != nullptr && !circular; ) {
                    child->m_is_taken = 1;
                    if(child->m_left_link == parent) {
                        parent = child;
                        child = child->m_right_link;
                    } else {
                        parent = child;
                        child = child->m_left_link;
                    }
                    circular = (child == &contig);
                }
                if(circular)
                    continue;
                parent = &contig;
                for(auto child = parent->m_left_link; child != nullptr; ) {
                    child->m_is_taken = 1;
                    if(child->m_left_link == parent) {
                        parent = child;
                        child = child->m_right_link;
                    } else {
                        parent = child;
                        child = child->m_left_link;
                    }
                }
            }

            list<function<void()>> jobs;
            for(int thr = 0; thr < ncores; ++thr) {
                jobs.push_back(bind(ConnectContigsJob, ref(scontigs)));
            }
            RunThreads(ncores, jobs);

            //remove fragments
            for(auto iloop = scontigs.begin(); iloop != scontigs.end(); ) {
                auto ic = iloop++;
                if(ic->m_is_taken == 2 || (int)ic->LenMin() < kmer_len)
                    scontigs.erase(ic);
                else
                    ic->m_is_taken = 0; 
            }
        }

        static void ConnectContigsJob(TContigList<DBGraph>& scontigs) {
            int kmer_len = scontigs.front().m_kmer_len;
            for(auto& contig : scontigs) {
                if(!contig.m_is_taken.Set(1))  // grab contig
                    continue;

                int num = 0;
                bool circular = false;
                for(auto parent = &contig; parent->m_right_link != nullptr && !circular; ++num) {
                    auto child = parent->m_right_link;
                    if(child->m_left_link != parent)
                        child->ReverseComplement();

                    if(child->m_right_link == &contig) {     // circular
                        circular = true;
                        if(child->m_left_link == &contig) {  // special case of a single connector needs additional check of orientation
                            if(contig.RightConnectingNode() != child->m_next_left)
                                child->ReverseComplement();
                        }
                    }

                    contig.AddToRight(*child);
                    if(num%2) // child is contig, not connector/extender
                        contig.m_seq.m_right_repeat = child->m_seq.m_right_repeat;
                    child->m_is_taken = 2;    // will be removed

                    if(circular && (int)contig.LenMax() >= 2*kmer_len-1)  //stabilize circular contig            
                        contig.RotateCircularToMinKmer();                    

                    parent = child;
                }
                if(circular)
                    continue;

                num = 0;
                for(auto parent = &contig; parent->m_left_link != nullptr; ++num) {
                    auto child = parent->m_left_link;
                    if(child->m_right_link != parent)
                        child->ReverseComplement();
                    contig.AddToLeft(*child);
                    if(num%2) // child is contig, not connector/extender
                        contig.m_seq.m_left_repeat = child->m_seq.m_left_repeat;
                    child->m_is_taken = 2;    // will be removed

                    parent = child;
                }

                //clip flanks which are not 'double' checked 
                auto& graph = contig.m_graph;

                for(int low_abundance_clip = 10; low_abundance_clip > 0 && contig.m_left_extend > 0; --low_abundance_clip) {
                    auto kmer = contig.FrontKmer();
                    if(!kmer.isValid() || graph.Abundance(kmer) > 5)
                        break;
                    
                    contig.ClipLeft(1);
                }
                int left_clip = min(kmer_len,contig.m_left_extend);
                contig.ClipLeft(left_clip);
                if(contig.m_left_extend > 0)
                    contig.m_seq.m_left_repeat = min(kmer_len-1, contig.m_left_extend+contig.m_seq.m_left_repeat);

                for(int low_abundance_clip = 10; low_abundance_clip > 0 && contig.m_right_extend > 0; --low_abundance_clip) {
                    auto kmer = contig.BackKmer();
                    if(!kmer.isValid() || graph.Abundance(kmer) > 5)
                        break;
                    
                    contig.ClipRight(1);
                }
                int right_clip = min(kmer_len,contig.m_right_extend); 
                contig.ClipRight(right_clip);
                if(contig.m_right_extend > 0)
                    contig.m_seq.m_right_repeat = min(kmer_len-1, contig.m_right_extend+contig.m_seq.m_right_repeat);
            }
        }

        CContigSequence m_seq;               // sequence

        Node m_next_left;          // denied left kmer (connection possible but it is already owned)
        Node m_next_right;         // denied right kmer (connection possible but it is already owned)

        SContig* m_left_link = nullptr;      // if set points to 'left' contig
        int m_left_shift = 0;                // shift+1 for m_next_left in this contig (positive for the right end)
        SContig* m_right_link = nullptr;     // if set points to 'right' contig
        int m_right_shift = 0;               // shift+1 for m_next_right in this contig (positive for the right end)

        int m_left_extend = 0;               // number of newly assembled bases which could be clipped
        int m_right_extend = 0;              // number of newly assembled bases which could be clipped

        DBGraph& m_graph;
        int m_kmer_len;
        SAtomic<uint8_t> m_is_taken = 0;
    };


    // This is a very lightweight class holding a reference to de Bruijn graph and main assembling parameters
    // It provides function used in assembling
    template<class DBGraph>
    class CDBGraphDigger {
    public:
        typedef typename DBGraph::Node Node;
        typedef typename DBGraph::Successor Successor;

        CDBGraphDigger(DBGraph& graph, double fraction, int jump, int low_count, bool allow_snps = false) : m_graph(graph), m_fraction(fraction), m_jump(jump), m_hist_min(graph.HistogramMinimum()), m_low_count(low_count), m_allow_snps(allow_snps) { 
            m_max_branch = 200; // maximum number of paths explored before quitting
        }

    private:
        typedef tuple<TStrList::iterator,int> TContigEnd;

    public:

        // starting from a node, find an extension of len l with maximal abundance
        string MostLikelyExtension(Node node, unsigned len) const { //don't do FilterNeighbors because it is called in it     
            string s;
            while(s.size() < len) {
                vector<Successor> successors = m_graph.GetNodeSuccessors(node);
                if(successors.empty())
                    return s;            
                sort(successors.begin(), successors.end(), [&](const Successor& a, const Successor& b) {return m_graph.Abundance(a.m_node) > m_graph.Abundance(b.m_node);}); 
                node = successors[0].m_node;
                s.push_back(successors[0].m_nt);
            }
            return s;
        }            

        string MostLikelySeq(Successor base, unsigned len) const {
            string s(1, base.m_nt);
            return s+MostLikelyExtension(base.m_node, len-1);
        }

        // starting from a node, find an extension of len l without forks (simple path); returns true if hit dead end
        pair<string, bool> StringentExtension(Node node, unsigned len) const {
            string s;
            while(s.size() < len) {
                vector<Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterNeighbors(successors, false);
                if(successors.empty())
                    return make_pair(s, true);
                if(successors.size() != 1)
                    return make_pair(s, false);            
                node = successors[0].m_node;
                s.push_back(successors[0].m_nt);
            }
            return make_pair(s, false);
        }

        bool ExtendableSuccessor(const Successor& initial_suc, double factor) const {
            int kmer_len = m_graph.KmerLen();
            int total_len = max(100, kmer_len);

            unordered_map<Node, int, typename Node::Hash> node_len;
            node_len.emplace(initial_suc.m_node,0);

            stack<pair<Node,int>> active_nodes;
            active_nodes.emplace(initial_suc.m_node,0);

            while(!active_nodes.empty()) {

                Node node = active_nodes.top().first;
                int len = active_nodes.top().second;
                active_nodes.pop();
                
                if(len == kmer_len) {
                    vector<Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(node));
                    FilterLowAbundanceNeighbors(step_back, factor);
                    bool found = false;
                    for(auto& back : step_back) {
                        if(back.m_nt == Complement(initial_suc.m_nt)) {
                            found = true;
                            break;
                        }
                    }
                    if(!found)
                        continue;
                }                

                if(len == total_len)
                    return true;

                if(len > kmer_len) {
                    int& l = node_len[node];
                    if(len > l)
                        l = len;
                    else
                        continue;
                }
        
                vector<Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterLowAbundanceNeighbors(successors, factor);

                if(!successors.empty()) {
                    for(int i = successors.size()-1; i >= 0; --i)
                        active_nodes.emplace(successors[i].m_node, len+1);
                } 
            }

            return false;
        }
        
        vector<Successor> GetReversibleNodeSuccessorsF(const Node& node, int* fork_infop, bool check_forward_extension, bool check_backward_extension) const {
            if(fork_infop != nullptr)
                *fork_infop = eNoFork;
            vector<Successor> neighbors = m_graph.GetNodeSuccessors(node);             
            FilterNeighbors(neighbors, check_forward_extension, 1.);
            if(neighbors.size() > 1 && fork_infop != nullptr && (LowCount() > 1 || m_graph.Abundance(neighbors[1].m_node) > 1))
                *fork_infop |= eRightFork;            

            for(auto& neighbor : neighbors) {
                vector<Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(neighbor.m_node));                
                FilterNeighbors(step_back, check_backward_extension, 1.);
                if(step_back.size() > 1 && fork_infop != nullptr && (LowCount() > 1 || m_graph.Abundance(step_back[1].m_node) > 1))
                    *fork_infop |= eLeftFork;                

                bool found = false;
                for(auto& back : step_back) {
                    if(back.m_node == m_graph.ReverseComplement(node)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    neighbors.clear();
                    return neighbors;
                }
            }

            return neighbors;
        }        

        vector<Successor> GetReversibleNodeSuccessors(const Node& node, int* numbackp = nullptr, bool check_extension = true) const {
            vector<Successor> neighbors = m_graph.GetNodeSuccessors(node);
            FilterNeighbors(neighbors, check_extension);
            if(numbackp != nullptr)
                *numbackp = 0;
            for(auto& neighbor : neighbors) {
                vector<Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(neighbor.m_node));
                FilterNeighbors(step_back, check_extension);
                bool found = false;
                for(auto& back : step_back) {
                    if(back.m_node == m_graph.ReverseComplement(node)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    neighbors.clear();
                    return neighbors;
                }
                if(numbackp != nullptr)
                    *numbackp = max(*numbackp, (int)step_back.size());
            }
            return neighbors;
        }

        int LowCount() const { return m_low_count; }
        bool GoodNode(const Node& node) const { return m_graph.Abundance(node) >= m_low_count; }
        int HistMin() const { return m_hist_min; }

        // removes noise forks
        void FilterLowAbundanceNeighbors(vector<Successor>& successors, double factor) const {
            // low abundance forks
            if(successors.size() > 1) {
                int abundance = 0;
                for(auto& suc : successors) {
                    abundance += m_graph.Abundance(suc.m_node);
                }
                sort(successors.begin(), successors.end(), [&](const Successor& a, const Successor& b) 
                     {
                         auto abundancea = m_graph.Abundance(a.m_node);
                         auto abundanceb = m_graph.Abundance(b.m_node);
                         if(abundancea == abundanceb)
                             return a.m_nt < b.m_nt;
                         else
                             return abundancea > abundanceb;
                     });
                if(LowCount() == 1 && m_graph.Abundance(successors.front().m_node) > 5) {
                    for(int j = successors.size()-1; j > 0 && m_graph.Abundance(successors.back().m_node) == 1; --j) 
                        successors.pop_back();            
                }
                for(int j = successors.size()-1; j > 0 && m_graph.Abundance(successors.back().m_node) <= m_fraction*abundance; --j) 
                    successors.pop_back();            
            }

            // strand specific noise reduction for Illumina issue of GGT->GG[ACG]
            if(m_graph.GraphIsStranded() && successors.size() > 1) {

                double fraction = factor*m_fraction;
            
                int target = -1;
                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(m_graph.GetNodeSeq(successors[j].m_node).substr(m_graph.KmerLen()-3) == "GGT") 
                        target = j;
                }
                if(target >= 0) {
                    int abundance = m_graph.Abundance(successors[target].m_node);
                    if(abundance > 5) {
                        double am = abundance*(1-m_graph.PlusFraction(successors[target].m_node));
                        for(int j = 0; j < (int)successors.size(); ) {
                            if(m_graph.Abundance(successors[j].m_node)*(1-m_graph.PlusFraction(successors[j].m_node)) < fraction*am)
                                successors.erase(successors.begin()+j);
                            else
                                ++j;
                        }
                    }
                }
            }

        }
        
        bool FilterNeighbors(vector<Successor>& successors, bool check_extension, double factor = 0.1) const {
            bool keep_fork_info = false;

            // low abundance forks
            FilterLowAbundanceNeighbors(successors, factor);

            //not extendable forks
            if(check_extension && successors.size() > 1 && m_graph.Abundance(successors.front().m_node) > 5) {
                for(int i = 0; i < (int)successors.size(); ) {
                    if(ExtendableSuccessor(successors[i], factor)) {
                        ++i;
                    } else {
                        successors.erase(successors.begin()+i);
                        keep_fork_info = true;
                    }
                }
            }

            // strand specific noise reduction for Illumina issue of GGT->GG[ACG] for negative strand and low coverage (the prev loop didn't work)
            if(m_graph.GraphIsStranded() && successors.size() > 1 && (!check_extension || m_graph.Abundance(successors.front().m_node) <= 5)) {

                double fraction = factor*m_fraction;
                int target = -1;

                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(MostLikelySeq(successors[j], 3) == "ACC") 
                        target = j;
                }
                if(target >= 0) {
                    int abundance =  m_graph.Abundance(successors[target].m_node);
                    if(abundance > 5) {
                        double ap = abundance*m_graph.PlusFraction(successors[target].m_node);
                        for(int j = 0; j < (int)successors.size(); ) {
                            if(m_graph.Abundance(successors[j].m_node)*m_graph.PlusFraction(successors[j].m_node) < fraction*ap)
                                successors.erase(successors.begin()+j);
                            else
                                ++j;
                        }
                    }
                }
            }

            // strand balance issue
            if(m_graph.GraphIsStranded() && successors.size() > 1) {

                double fraction = 0.1*m_fraction;
                                      
                bool has_both = false;
                for(int j = 0; !has_both && j < (int)successors.size(); ++j) {
                    double plusf = m_graph.PlusFraction(successors[j].m_node);
                    double minusf = 1.- plusf;
                    has_both = GoodNode(successors[j].m_node) && (min(plusf,minusf) > 0.25);
                }

                if(has_both) {
                    for(int j = 0; j < (int)successors.size(); ) {
                        double plusf = m_graph.PlusFraction(successors[j].m_node);
                        double minusf = 1.- plusf;
                        int abundancej = m_graph.Abundance(successors[j].m_node);
                        if(abundancej > 1 && min(plusf,minusf) < fraction*max(plusf,minusf))
                            successors.erase(successors.begin()+j);
                        else
                            ++j;
                    }
                } 
            }

            return keep_fork_info;
        }

        DBGraph& Graph() { return m_graph; }
        double Fraction() const { return m_fraction; }

        enum EConnectionStatus {eSuccess, eNoConnection, eAmbiguousConnection};

        struct SElement {
            SElement (Successor suc, SElement* link) : m_link(link), m_suc(suc) {}
            struct SElement* m_link;   // previous element  
            Successor m_suc;
        };
        // connects two nodes in a finite number of steps
        pair<TBases<DBGraph>, EConnectionStatus> ConnectTwoNodes(const Node& first_node, const Node& last_node, int steps) const {
        
            pair<TBases<DBGraph>, EConnectionStatus> bases(TBases<DBGraph>(), eNoConnection);

            deque<SElement> storage; // will contain ALL extensions (nothing is deleted)    
            typedef unordered_map<Node, SElement*, typename Node::Hash> TElementMap;  //pointer to its own element OR zero if ambiguous path    
            TElementMap current_elements;

            vector<Successor> successors = m_graph.GetNodeSuccessors(first_node);
            FilterNeighbors(successors, false);
            for(auto& suc : successors) {
                storage.push_back(SElement(suc, 0));
                current_elements[suc.m_node] = &storage.back();
            }

            list<SElement> connections;
            for(int step = 1; step < steps && !current_elements.empty(); ++step) {
                TElementMap new_elements;
                for(auto& el : current_elements) {
                    vector<Successor> successors = m_graph.GetNodeSuccessors(el.first);
                    FilterNeighbors(successors, false);
                    if(el.second == 0) {  // ambiguous path 
                        for(auto& suc : successors) {
                            new_elements[suc.m_node] = 0;
                            if(suc.m_node == last_node) {
                                bases.second = eAmbiguousConnection;
                                return bases;
                            }
                        }
                    } else {
                        for(auto& suc : successors) {
                            storage.push_back(SElement(suc, el.second));
                            if(suc.m_node == last_node) {
                                if(!connections.empty()) {
                                    bases.second = eAmbiguousConnection;
                                    return bases;
                                } else {
                                    connections.push_back(storage.back()); 
                                }                           
                            }
                            pair<typename TElementMap::iterator, bool> rslt = new_elements.insert(make_pair(suc.m_node, &storage.back()));
                            if(!rslt.second || !GoodNode(suc.m_node))
                                rslt.first->second = 0;
                        }
                    }                    
                }
                swap(current_elements, new_elements);
                if(current_elements.size() > m_max_branch)
                    return bases;
            }

            if(connections.empty())
                return bases;

            SElement el = connections.front();
            while(el.m_link != 0) {
                bases.first.push_front(el.m_suc);
                el = *el.m_link;
            }
            bases.first.push_front(el.m_suc);
            bases.second = eSuccess;
            return bases;
        }

        typedef list<TBases<DBGraph>> TBasesList;
        typedef unordered_map<Node, forward_list<typename TBasesList::iterator>, typename Node::Hash> TBranch;  // all 'leaves' will have the same length  
        typedef unordered_map<Node, forward_list<pair<typename TBasesList::iterator, int>>, typename Node::Hash> TLinks;

        void OneStepBranchExtend(TBranch& branch, TBasesList& sequences, TLinks& links) {
            TBranch new_branch;
            for(auto& leaf : branch) {
                vector<Successor> successors = m_graph.GetNodeSuccessors(leaf.first);
                FilterNeighbors(successors, true);
                if(successors.empty()) {
                    for(auto is : leaf.second)
                        is->clear();
                } else {
                    for(int i = successors.size()-1; i >= 0; --i) {
                        auto& lst = new_branch[successors[i].m_node];
                        for(auto is : leaf.second) {
                            if(i > 0) {  // copy sequence if it is a fork                   
                                sequences.push_front(*is);
                                is = sequences.begin();
                                for(int p = 0; p < (int)is->size()-1; ++p)
                                    links[(*is)[p].m_node].emplace_front(is, p);
                            }
                            links[leaf.first].emplace_front(is, is->size()-1);
                            is->push_back(successors[i]);
                            lst.emplace_front(is);
                        }
                    }
                }
            }

            for(auto it_loop = new_branch.begin(); it_loop != new_branch.end(); ) {
                auto it = it_loop++;
                auto rslt = links.find(it->first);
                if(rslt != links.end()) {
                    auto& lst = rslt->second;
                    set<TBases<DBGraph>> seqs; // TODO set of intervals
                    for(auto& link : lst) {
                        if(!link.first->empty()) {
                            if(link.first->back().m_node != it->first) {
                                seqs.emplace(link.first->begin()+link.second+1, link.first->end());
                            } else { // circular extension
                                branch.clear();
                                sequences.clear();
                                return;
                            }
                        }
                    }
                    if(!seqs.empty()) {
                        for(auto is : it->second) {
                            for(auto ex = next(seqs.begin()); ex != seqs.end(); ++ex) {
                                sequences.push_front(*is);
                                sequences.front().insert(sequences.front().end(), ex->begin(), ex->end());
                                for(int p = 0; p < (int)sequences.front().size()-1; ++p)
                                    links[sequences.front()[p].m_node].emplace_front(sequences.begin(), p);
                                new_branch[sequences.front().back().m_node].push_front(sequences.begin());
                            }
                            int l = is->size();
                            is->insert(is->end(), seqs.begin()->begin(), seqs.begin()->end());
                            for(int p = l-1; p < (int)is->size()-1; ++p)
                                links[(*is)[p].m_node].emplace_front(is, p);
                            new_branch[is->back().m_node].push_front(is);
                        }
                        new_branch.erase(it);
                    }
                }
            }

            swap(branch, new_branch);
        }
        

        // sequences, last node, intrusion, node corresponding to intrusion shift, max_le - min_le
        tuple<TLocalVariants, Node, int, Node, int> DiscoverOneSNP(const vector<Successor>& successors, const TVariation& last_chunk, int max_extent) {
            tuple<TLocalVariants, Node, int, Node, int> rslt;
            TBranch extensions;
            TBasesList sequences; // assembled seqs 
            TLinks links;
            int kmer_len = m_graph.KmerLen();

            if(max_extent == 0 || successors.empty())
                return rslt;

            for(auto& suc : successors) {
                sequences.emplace_front(1,suc);
                extensions[suc.m_node].emplace_front(sequences.begin());         
            }

            int max_len = 1;
            int min_len = 1;
            size_t seq_num = sequences.size();
            while(seq_num < m_max_branch && max_len < max_extent) {
                OneStepBranchExtend(extensions, sequences, links);
                max_len = 0;
                min_len = numeric_limits<int>::max();
                seq_num = 0;
                for(auto& seq : sequences) {
                    if(!seq.empty()) {
                        max_len = max(max_len, (int)seq.size());
                        min_len = min(min_len, (int)seq.size());
                        ++seq_num;
                    }
                }

                if(extensions.empty())  // can't extend 
                    return rslt;                
                if(extensions.size() == 1 && min_len >= kmer_len)
                    break;
            }

            if(extensions.size() == 1 && min_len >= kmer_len && max_len <= max_extent) {
                set<char> first_bases;
                for(auto it = sequences.begin(); it != sequences.end(); ) {
                    if(it->empty()) {
                        it = sequences.erase(it);
                    } else {
                        first_bases.insert(it->front().m_nt);
                        ++it;
                    }
                }
                if(first_bases.size() > 1) {  // found snp
                    // clip extra matches from the end
                    int matches = 0;
                    bool all_same = true;
                    while(all_same) {
                        for(auto& seq : sequences) {
                            if(matches == (int)seq.size() || (seq.end()-matches-1)->m_nt != (sequences.front().end()-matches-1)->m_nt) {
                                all_same = false;
                                break;
                            }
                        }
                        if(all_same)
                            ++matches;
                    }
                    if(matches > kmer_len) {
                        int extra = min(matches-kmer_len, max_len);
                        max_len -= extra;
                        min_len -= extra;
                        for(auto& seq : sequences)
                            seq.erase(seq.end()-extra, seq.end());
                    }                

                    // check all nodes  
                    for(auto& seq : sequences) {
                        for(auto& base : seq) {
                            if(!GoodNode(base.m_node))
                                return rslt;
                        }
                    }

                    bool has_empty_variant = false;
                    for(auto& seq : sequences) {
                        if((int)seq.size() == kmer_len) {
                            has_empty_variant = true;
                            break;
                        }
                    }

                    // copy seqs to result  
                    TLocalVariants& seqs = get<0>(rslt);
                    for(auto& seq : sequences) {
                        seqs.push_front(TVariation());
                        for(auto& base : seq)
                            seqs.front().push_back(base.m_nt);
                    }
                    // last node    
                    get<1>(rslt) = sequences.front().back().m_node;
                    // diff
                    get<4>(rslt) = max_len-min_len;
                
                    // check for repeat and report if found         
                    if(has_empty_variant) {
                        for(auto& seq : seqs) // add kmer_len bases     
                            seq.insert(seq.begin(), last_chunk.end()-kmer_len, last_chunk.end());

                        bool all_same = true;
                        int shift = 0;
                        while(all_same) {
                            for(auto& seq : seqs) {
                                if(shift == (int)seq.size()-kmer_len || *(seq.end()-shift-1-kmer_len) != *(seqs.front().end()-shift-1-kmer_len)) {
                                    all_same = false;
                                    break;
                                }
                            } 
                            if(all_same)
                                ++shift;
                        } 
                        if(shift >= kmer_len) {
                            return tuple<TLocalVariants, Node, int, Node, int>();
                        } else {
                            get<2>(rslt) = shift;
                            get<3>(rslt) = (sequences.front().end()-1-shift)->m_node;

                            max_len += shift;
                            min_len += shift;
                            for(auto& seq : seqs) // erase added extra sequence (shit bases remain) 
                                seq.erase(seq.begin(), seq.begin()+kmer_len-shift);
                        }
                    }
                }
            }

            return rslt;            
        }

        //successors by value intentionally
        tuple<TLocalVariants, Node, int> DiscoverSNPCluster(vector<Successor> successors, const TVariation& last_chunk, int max_extent) {
            tuple<TLocalVariants, Node, int> rslt;

            int kmer_len = m_graph.KmerLen();
            const TVariation* last_chunkp = &last_chunk;

            int dist_to_snp = 0;
            while(dist_to_snp < 2*kmer_len && !successors.empty()) {
                tuple<TLocalVariants, Node, int, Node, int> snp_data = DiscoverOneSNP(successors, *last_chunkp, max_extent);
                if(get<0>(snp_data).empty())
                    break;

                TLocalVariants& seqs = get<0>(snp_data);
                Node node = get<1>(snp_data);
                int shift = get<2>(snp_data);
                int diff_len = get<4>(snp_data);

                /*
                cerr << "Last chunk: ";
                for(char c : *last_chunkp)
                    cerr << c;
                cerr << endl;
                cerr << "SNP: " << shift << " " << diff_len << " " << m_graph.GetNodeSeq(node) << endl;
                for(auto& seq: seqs) {
                    cerr << "Chunk: ";
                    for(char c : seq)
                        cerr << c;
                    cerr << endl;
                    auto rseq = seq;
                    ReverseComplementSeq(rseq.begin(), rseq.end());
                    cerr << "RChunk: ";
                    for(char c : rseq)
                        cerr << c;
                    cerr << endl;

                }
                */                                                                                                         

                if(dist_to_snp == 0) {                // first snp in cluster
                    get<0>(rslt) = seqs;
                    get<1>(rslt) = node;
                    get<2>(rslt) = shift;
                } else {
                    int& existing_shift = get<2>(rslt);
                    if(dist_to_snp >= kmer_len+shift &&                           // take into account repeat (if any)
                       dist_to_snp+existing_shift >= kmer_len+diff_len-shift)     // long indels need additional steps before they are connected
                        break;                                                    // no inerference
                                                                                  // combine snps 
                    int len = 0;
                    for(auto& seq: get<0>(rslt)) {
                        seq.erase(seq.end()-shift, seq.end());                    // erase seq for new shift
                        if(existing_shift > 0)
                            seq.erase(seq.begin(), seq.begin()+existing_shift);   // remove existing shift (if any)
                        for(auto it = next(seqs.begin()); it != seqs.end(); ++it) {
                            get<0>(rslt).push_front(seq);
                            get<0>(rslt).front().insert(get<0>(rslt).front().end(), it->begin(), it->end()-shift);
                        }
                        seq.insert(seq.end(), seqs.front().begin(), seqs.front().end()-shift);
                        len = max(len, (int)seq.size());
                    }
                    if(len > max_extent || (size_t)distance(get<0>(rslt).begin(), get<0>(rslt).end()) > m_max_branch)  // longer than threshold ot too many variants
                        return tuple<TLocalVariants, Node, int>();
                    existing_shift = 0;
                    if(shift > 0)
                        node = get<3>(snp_data);
                    get<1>(rslt) = node;
                }
                dist_to_snp = kmer_len;

                //                cerr << "Ext: ";

                bool fork = false;
                while(dist_to_snp < 2*kmer_len) {
                    successors = m_graph.GetNodeSuccessors(node);                    
                    FilterNeighbors(successors, true);
                    fork = successors.size() > 1;
                    if(fork || successors.empty() || !GoodNode(successors.front().m_node))
                        break;                    

                    ++dist_to_snp;
                    node = successors.front().m_node;
                    for(auto& seq: get<0>(rslt))
                        seq.push_back(successors.front().m_nt);

                    //                    cerr << successors.front().m_nt;
                }

                //                cerr << endl;

                if(!fork)
                    break;

                last_chunkp = &get<0>(rslt).front();
            }

            if(dist_to_snp > kmer_len) { // remove extra extension
                for(auto& seq: get<0>(rslt))
                    seq.erase(seq.end()-(dist_to_snp-kmer_len), seq.end());
            }
            
            return rslt;            
        }

        // starting from initial_node assembles the right extension        
        tuple<SContig<DBGraph>, Node, int> ExtendToRight(const Node& initial_node, int allowed_intrusion) { // initial_node may be not owned 
            Node node = initial_node;
            SContig<DBGraph> extension(m_graph);
            int max_extent = m_jump;
            int kmer_len = m_graph.KmerLen();
            int initial_node_intrusion = 0;

            while(true) {
                vector<Successor> successors = m_graph.GetNodeSuccessors(node);                    
                FilterNeighbors(successors, true);
                if(successors.empty()) {                                           // no extensions 
                    break;  
                } else if(successors.size() == 1) {                                // simple extension 
                    Node new_node = successors.front().m_node;
                    if(!GoodNode(new_node))
                        break;
                    vector<Successor> predecessors = m_graph.GetNodeSuccessors(DBGraph::ReverseComplement(new_node));
                    FilterNeighbors(predecessors, true);
                    if(predecessors.size() != 1)                                   // no extensions  or end of unique seq before repeat 
                        break;
                    if(DBGraph::ReverseComplement(predecessors[0].m_node) != node) // no return
                        break;

                    node = new_node;
                    if(m_graph.SetVisited(node)) {                                 // node is available
                        if(extension.m_seq.empty()) {
                            extension.m_seq.InsertNewChunk();
                            extension.m_seq.InsertNewVariant();
                        }
                        extension.m_seq.ExtendTopVariant(successors.front().m_nt);
                    } else {
                        return make_tuple(extension, node, initial_node_intrusion); 
                    }                   
                } else if(!m_allow_snps) {                                        // snps not allowed
                    break;
                } else {                                                          // try snps
                    int last_chunk_len = extension.m_seq.empty() ? 0 : extension.m_seq.back().front().size();
                    TVariation* last_chunkp = nullptr;
                    TVariation last_chunk;
                    if(last_chunk_len >= kmer_len) {
                        last_chunkp = &extension.m_seq.back().front();
                    } else {
                        string initial_node_seq = m_graph.GetNodeSeq(initial_node);
                        last_chunk.insert(last_chunk.end(), initial_node_seq.begin(), initial_node_seq.end());
                        if(last_chunk_len > 0)
                            last_chunk.insert(last_chunk.end(), extension.m_seq.back().front().begin(), extension.m_seq.back().front().end());
                        last_chunkp = &last_chunk;
                    }

                    //                    cerr << "Direct" << endl;

                    auto forward = DiscoverSNPCluster(successors, *last_chunkp, max_extent);
                    TLocalVariants& step = get<0>(forward);
                    int shift = get<2>(forward);

                    if(step.empty())        // no snp
                        break;
                    int step_size = 0;  // not all step seqs have same length
                    for(auto& var : step)
                        step_size = max(step_size, (int)var.size());                

                    // check return
                    vector<Successor> predecessors = m_graph.GetNodeSuccessors(DBGraph::ReverseComplement(get<1>(forward)));

                    //                    cerr << "Lastkmer: " << get<1>(forward).isValid() << " " << m_graph.GetNodeSeq(DBGraph::ReverseComplement(get<1>(forward))) << endl;

                    FilterNeighbors(predecessors, true);
                    if(predecessors.empty())
                        break;
                    TVariation back_chunk;
                    back_chunk.insert(back_chunk.end(), step.front().end()-kmer_len, step.front().end());
                    ReverseComplementSeq(back_chunk.begin(), back_chunk.end());                    

                    //                    cerr << "Backward" << endl;

                    auto backward = DiscoverSNPCluster(predecessors, back_chunk, max_extent);
                    TLocalVariants& step_back = get<0>(backward);

                    if(step_back.empty())        // no snp
                        break;
                    int step_back_size = 0;  // not all step seqs have same length
                    for(auto& var : step_back)
                        step_back_size = max(step_back_size, (int)var.size());  
                    if(step_size != step_back_size)
                        break;
                    for(auto& seq : step_back) 
                        ReverseComplementSeq(seq.begin(), seq.end());

                    if(!equal(last_chunkp->end()-kmer_len, last_chunkp->end()-shift, step_back.front().begin()+shift))
                        break;
                    for(auto& seq : step_back) {
                        seq.erase(seq.begin(), seq.begin()+kmer_len);
                        seq.insert(seq.end(), step.front().end()-kmer_len, step.front().end());
                    }
                    step.sort();
                    step_back.sort();

                    if(step != step_back)
                        break; 

                    // snp is accepted
                    node = get<1>(forward);

                    if(shift > 0) {
                        if(shift >= last_chunk_len) { // extension has no snp and short - pass intrusion (or part of it) to the caller
                            initial_node_intrusion = shift-last_chunk_len;
                            if(initial_node_intrusion > allowed_intrusion) {
                                initial_node_intrusion = 0;
                                break;
                            }
                            if(last_chunk_len > 0)
                                extension.m_seq.pop_back(); 
                        } else {       // shorten previous chunk  
                            extension.m_seq.back().front().erase(extension.m_seq.back().front().end()-shift, extension.m_seq.back().front().end());
                        }
                    }

                    extension.m_seq.InsertNewChunk();                         // empty chunk for variable part
                    for(auto& seq : step) {
                        extension.m_seq.InsertNewVariant();                   // empty seq for new variant    
                        extension.m_seq.ExtendTopVariant(seq.begin(), seq.end()-kmer_len);
                    }
                    extension.m_seq.InsertNewChunk();                         // empty chunk for matching kmer-1 or kmer bases    
                    extension.m_seq.InsertNewVariant(step.front().end()-kmer_len, step.front().end()-1); 

                    CReadHolder rh(false);
                    for(auto& seq : step) {
                        if(shift == 0)                                                                         // last chunk not clipped, nothing added to snp
                            seq.insert(seq.begin(), last_chunkp->end()-(kmer_len-1), last_chunkp->end());
                        else if(shift > last_chunk_len)                                                        // last chunk not clipped, shift bases added to snp
                            seq.insert(seq.begin(), last_chunkp->end()-(kmer_len-1), last_chunkp->end()-shift);
                        else                                                                                   // last chunk clipped, shift bases added to snp
                            seq.insert(seq.begin(), last_chunkp->end()-(kmer_len-1-shift), last_chunkp->end());
                        rh.PushBack(seq);
                    }
                    bool my_snp = true;
                    list<Node> snp_nodes;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend() && my_snp; ++ik) { 
                        Node n = m_graph.GetNode(*ik);
                        if(n.isValid()) {
                            if(m_graph.IsVisited(n))
                                my_snp = false;
                            else
                                snp_nodes.push_back(n);
                        }
                    }
                    if(my_snp && m_graph.SetVisited(node)) {                            // snp belongs to this thread - extend one more base and set all visited
                        extension.m_seq.ExtendTopVariant(step.front().back());
                        for(auto& n : snp_nodes)
                            m_graph.SetVisited(n);                        
                        
                        continue;
                    } else {

                        return make_tuple(extension, node, initial_node_intrusion);
                    }                
                }
            }


            return make_tuple(extension, Node(), initial_node_intrusion);
        }


        // assembles a contig starting from initial_node 
        // min_len - minimal length for accepted contigs
        // changes the state of all used nodes to 'visited' or 'temporary holding'   
        SContig<DBGraph> GetContigForKmer(const Node& initial_node, int min_len) {
            if(m_graph.Abundance(initial_node) < m_hist_min || !GoodNode(initial_node) || !m_graph.SetVisited(initial_node))
                return SContig<DBGraph>(m_graph);

            //node is good and this thread owns it  

            // don't allow intrusion of snps in the initial kmer
            tuple<SContig<DBGraph>, Node, int> to_right = ExtendToRight(initial_node, 0);
            tuple<SContig<DBGraph>, Node, int> to_left = ExtendToRight(DBGraph::ReverseComplement(initial_node), 0);

            SContig<DBGraph> scontig(get<0>(to_left), get<0>(to_right), initial_node, DBGraph::ReverseComplement(get<1>(to_left)), get<1>(to_right), m_graph);
            
            if(!scontig.m_next_left.isValid() && !scontig.m_next_right.isValid() && (int)scontig.LenMin() < min_len) {
                int kmer_len = m_graph.KmerLen();
                for(int i = scontig.m_seq.size()-1; i >= 0; i -= 2) {
                    if(i == (int)scontig.m_seq.size()-1) { // last chunk size >= kmer_len
                        CReadHolder rh(false);
                        rh.PushBack(scontig.m_seq.back().front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                            m_graph.SetTempHolding(m_graph.GetNode(*ik));
                    } else {
                        if((int)scontig.m_seq.ChunkLenMax(i) >= kmer_len) {
                            TVariation seq(scontig.m_seq[i].front().begin(), scontig.m_seq[i].front().end());
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                                m_graph.SetTempHolding(m_graph.GetNode(*ik));
                        }
                        for(auto& variant : scontig.m_seq[i+1]) {
                            TVariation seq(scontig.m_seq[i].front().end()-kmer_len+1, scontig.m_seq[i].front().end()); // all uniq chunks >= kmer_len-1
                            seq.insert(seq.end(), variant.begin(), variant.end());
                            seq.insert(seq.end(), scontig.m_seq[i+2].front().begin(), scontig.m_seq[i+2].front().begin()+kmer_len-1);
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                                m_graph.SetTempHolding(m_graph.GetNode(*ik));
                        }
                    }
                }

                return SContig<DBGraph>(m_graph);
            } else {
                return scontig;
            }
        }

        void CheckRepeats(TContigList<DBGraph>& scontigs) {
            int kmer_len = m_graph.KmerLen();

            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                auto& contig = it->m_seq;

                if(contig.m_left_repeat >= kmer_len && contig.m_left_repeat < (int)contig.LenMin()) {
                    int last_chunk = 0;
                    for(int len = contig.ChunkLenMin(last_chunk); len < contig.m_left_repeat+1; len += contig.ChunkLenMin(++last_chunk));

                    int check_len = contig.m_left_repeat+1;
                    vector<forward_list<Node>> kmers(check_len-kmer_len+1);

                    stack<pair<TVariation*, int>> active_chunks;
                    active_chunks.emplace(&contig[0].front(), 0);
                    deque<TVariation*> current_seqs;
                    while(!active_chunks.empty()) {
                        TVariation* seqp = active_chunks.top().first;
                        int chunk_num = active_chunks.top().second;
                        active_chunks.pop();

                        current_seqs.resize(chunk_num);
                        current_seqs.push_back(seqp);
                        for(int chunk = chunk_num+1; chunk <= last_chunk; ++chunk) {
                            auto it = contig[chunk].begin();
                            current_seqs.push_back(&(*it));
                            for(++it; it != contig[chunk].end(); ++it) 
                                active_chunks.emplace(&(*it), chunk);
                        }
                        list<char> seq;
                        int seq_len = 0;
                        for(auto ics = current_seqs.begin(); seq_len < check_len; ++ics) {
                            int extra = max(0, seq_len+(int)(*ics)->size()-check_len);
                            seq.insert(seq.end(), (*ics)->begin(), (*ics)->end()-extra);
                            seq_len += (*ics)->size()-extra;
                        }                            
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = kmers.size()-1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                            Node node = m_graph.GetNode(*ik);
                            if(find(kmers[pos].begin(), kmers[pos].end(), node) == kmers[pos].end())
                                kmers[pos].push_front(node);
                        }
                    }
                    
                    for( ; contig.m_left_repeat >= kmer_len; --contig.m_left_repeat) {
                        int p = contig.m_left_repeat-kmer_len;

                        bool bad_node = false;
                        for(auto& kmer : kmers[p]) {
                            if(!kmer.isValid() || !GoodNode(kmer)) {
                                bad_node = true;
                                break;
                            }
                        }
                        if(bad_node)
                            break;
                                                 
                        bool no_step = false;
                        for(auto& kmer : kmers[p]) {
                            if(kmer.isValid()) {
                                vector<Successor> successors = m_graph.GetNodeSuccessors(kmer);
                                FilterNeighbors(successors, true);
                                if(successors.empty()) {
                                    no_step = true;
                                    break;
                                }
                                auto& next_lst = kmers[p+1];
                                for(auto& suc : successors) {
                                    if(find_if(next_lst.begin(), next_lst.end(), [suc](const Node& node) {return node == suc.m_node; }) == next_lst.end()) {
                                        no_step = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if(no_step)
                            break;
                    }                    
                }
            

                if(contig.m_right_repeat >= kmer_len && contig.m_right_repeat < (int)contig.LenMin()) {
                    int first_chunk = contig.size()-1;
                    for(int len = contig.ChunkLenMin(first_chunk); len < contig.m_right_repeat+1; len += contig.ChunkLenMin(--first_chunk));

                    int check_len = contig.m_right_repeat+1;
                    vector<forward_list<Node>> kmers(check_len-kmer_len+1);

                    stack<pair<TVariation*, int>> active_chunks;
                    active_chunks.emplace(&contig[contig.size()-1].front(), contig.size()-1);
                    deque<TVariation*> current_seqs;
                    while(!active_chunks.empty()) {
                        TVariation* seqp = active_chunks.top().first;
                        int chunk_num = active_chunks.top().second;
                        active_chunks.pop();

                        if(!current_seqs.empty())
                            current_seqs.erase(current_seqs.begin(), current_seqs.begin()+chunk_num-first_chunk+1);
                        current_seqs.push_front(seqp);
                        for(int chunk = chunk_num-1; chunk >= first_chunk; --chunk) {
                            auto it = contig[chunk].begin();
                            current_seqs.push_front(&(*it));
                            for(++it; it != contig[chunk].end(); ++it) 
                                active_chunks.emplace(&(*it), chunk);
                        }
                        list<char> seq;
                        int seq_len = 0;
                        for(auto ics = current_seqs.rbegin(); seq_len < check_len; ++ics) {
                            int extra = max(0, seq_len+(int)(*ics)->size()-check_len);
                            seq.insert(seq.begin(), (*ics)->begin()+extra, (*ics)->end());
                            seq_len += (*ics)->size()-extra;
                        }                            
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = kmers.size()-1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                            Node node = m_graph.GetNode(*ik);
                            if(find(kmers[pos].begin(), kmers[pos].end(), node) == kmers[pos].end())
                                kmers[pos].push_front(node);
                        }
                    }
                    
                    for( ; contig.m_right_repeat >= kmer_len; --contig.m_right_repeat) {
                        int p = kmers.size()-(contig.m_right_repeat-kmer_len+1);

                        bool bad_node = false;
                        for(auto& kmer : kmers[p]) {
                            if(!kmer.isValid() || !GoodNode(kmer)) {
                                bad_node = true;
                                break;
                            }
                        }
                        if(bad_node)
                            break;

                        bool no_step = false;
                        for(auto& kmer : kmers[p]) {
                            if(kmer.isValid()) {
                                vector<Successor> successors = m_graph.GetNodeSuccessors(DBGraph::ReverseComplement(kmer));
                                FilterNeighbors(successors, true);
                                if(successors.empty()) {
                                    no_step = true;
                                    break;
                                }
                                auto& prev_lst = kmers[p-1];
                                for(auto& suc : successors) {
                                    if(find_if(prev_lst.begin(), prev_lst.end(), [suc](const Node& node) {return node == DBGraph::ReverseComplement(suc.m_node); }) == prev_lst.end()) {
                                        no_step = true;
                                        break;
                                    }
                                }
                            }
                        }                        
                        if(no_step)
                            break;
                    }                    
                }
            }
        }

        void ConnectOverlappingContigs(TContigList<DBGraph>& scontigs) {
            int kmer_len = m_graph.KmerLen();
            unordered_map<Node, forward_list<pair<typename TContigList<DBGraph>::iterator, int>>, typename Node::Hash> kmers;
            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                SContig<DBGraph>& contig = *it;
                if((int)contig.m_seq.ChunkLenMax(0) > kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig.m_seq[0].front());
                    int pos = contig.m_seq.ChunkLenMax(0)-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                        if(pos < (int)contig.m_seq.ChunkLenMax(0)-kmer_len) {
                            Node node = m_graph.GetNode(*ik);
                            if(node.isValid())
                                kmers[node].emplace_front(it, pos);
                        }
                    }
                }
                if(contig.m_seq.size() > 1 && (int)contig.m_seq.ChunkLenMax(contig.m_seq.size()-1) > kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig.m_seq[contig.m_seq.size()-1].front());
                    int pos = contig.m_seq.LenMax()-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                        if(pos > (int)contig.m_seq.LenMax()-(int)contig.m_seq.ChunkLenMax(contig.m_seq.size()-1)) {
                            Node node = m_graph.GetNode(*ik);
                            if(node.isValid())
                                kmers[node].emplace_front(it, pos);
                        }
                    }
                }
            }

            list<tuple<typename TContigList<DBGraph>::iterator, typename TContigList<DBGraph>::iterator, int, int, int>> overlaps; // first contig, second contig, start/end, start/end, len    

            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                SContig<DBGraph>& icontig = *it;

                // right overlap    
                {
                    list<tuple<typename TContigList<DBGraph>::iterator, typename TContigList<DBGraph>::iterator, int, int, int>> contig_overlaps;

                    auto& irchunk = icontig.m_seq.back().front();
                    auto rslt = kmers.find(icontig.BackKmer()); // rightend to left end 
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = hit.second+kmer_len;
                            auto& jlchunk = jt->m_seq.front().front();
                            if(overlap_len > (int)irchunk.size() || overlap_len > (int)jlchunk.size())
                                continue;
                            if(!equal(jlchunk.begin(), jlchunk.begin()+hit.second, irchunk.end()-overlap_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, 1, -1, overlap_len);
                        }
                    }
                    rslt = kmers.find(DBGraph::ReverseComplement(icontig.BackKmer())); // right end to right end   
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = jt->m_seq.LenMax()-hit.second;
                            auto& jrchunk = jt->m_seq.back().front();
                            if(overlap_len > (int)irchunk.size() || overlap_len > (int)jrchunk.size())
                                continue;
                            TVariation seq(irchunk.end()-overlap_len, irchunk.end()-kmer_len);
                            ReverseComplementSeq(seq.begin(), seq.end());
                            if(!equal(seq.begin(), seq.end(), jrchunk.end()-overlap_len+kmer_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, 1, 1, overlap_len);
                        }
                    }
                    if(contig_overlaps.size() == 1)
                        overlaps.splice(overlaps.end(), contig_overlaps);
                }

                //left overlap  
                {
                    list<tuple<typename TContigList<DBGraph>::iterator, typename TContigList<DBGraph>::iterator, int, int, int>> contig_overlaps;

                    auto& ilchunk = it->m_seq.front().front();
                    auto rslt = kmers.find(icontig.FrontKmer()); // left end to right end   
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = jt->m_seq.LenMax()-hit.second;
                            auto& jrchunk = jt->m_seq.back().front();
                            if(overlap_len > (int)ilchunk.size() || overlap_len > (int)jrchunk.size())
                                continue;
                            if(!equal(jrchunk.end()-overlap_len+kmer_len, jrchunk.end(), ilchunk.begin()+kmer_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, -1, 1, overlap_len);
                        }
                    }
                    rslt = kmers.find(DBGraph::ReverseComplement(icontig.FrontKmer())); // left end to left end    
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = hit.second+kmer_len;
                            auto& jlchunk = jt->m_seq.front().front();
                            if(overlap_len > (int)ilchunk.size() || overlap_len > (int)jlchunk.size())
                                continue;
                            TVariation seq(ilchunk.begin()+kmer_len, ilchunk.begin()+overlap_len);
                            ReverseComplementSeq(seq.begin(), seq.end());
                            if(!equal(jlchunk.begin(), jlchunk.begin()+hit.second, seq.begin()))
                                continue;

                            contig_overlaps.emplace_back(it, jt, -1, -1, overlap_len);
                        }
                    }
                    if(contig_overlaps.size() == 1)
                        overlaps.splice(overlaps.end(), contig_overlaps);
                }
            }
                
            for(auto it = overlaps.begin(); it != overlaps.end(); ) {
                auto overlap = *it;
                swap(get<0>(overlap), get<1>(overlap));
                swap(get<2>(overlap), get<3>(overlap));
                auto jt = find(it, overlaps.end(), overlap);
                if(jt == overlaps.end()) {
                    auto tmp = it++;
                    overlaps.erase(tmp);
                } else {
                    overlaps.erase(jt);
                    ++it;
                }
            }                


            for(auto it_loop = overlaps.begin(); it_loop != overlaps.end(); ) {
                auto it = it_loop++;
                auto& overlap = *it;
                int overlap_len = get<4>(overlap);
                auto icontigp = get<0>(overlap);
                auto jcontigp = get<1>(overlap);
                int diri = get<2>(overlap);
                int dirj = get<3>(overlap);
                
                auto NextIBase = [&]() {
                    Node node = diri > 0 ? icontigp->BackKmer(): DBGraph::ReverseComplement(icontigp->FrontKmer());
                    auto forward = m_graph.GetNodeSuccessors(node);
                    FilterNeighbors(forward, true);
                    if(forward.size() == 1) {
                        auto backward = m_graph.GetNodeSuccessors(DBGraph::ReverseComplement(forward.front().m_node));
                        FilterNeighbors(backward, true);
                        if(backward.size() == 1 && DBGraph::ReverseComplement(backward.front().m_node) == node)
                            return forward.front().m_nt;                            
                    }
                    return 'N';
                };
                auto NextJBase = [&]() {
                    return dirj < 0 ? *(jcontigp->m_seq.front().front().begin()+overlap_len) : Complement(*(jcontigp->m_seq.back().front().end()-overlap_len-1));
                };

                
                bool connected;
                if(diri > 0)
                    connected = (icontigp->m_seq.m_right_repeat < kmer_len);
                else
                    connected = (icontigp->m_seq.m_left_repeat < kmer_len);
                if(dirj > 0)
                    connected = connected && (jcontigp->m_seq.m_right_repeat < kmer_len);
                else
                    connected = connected && (jcontigp->m_seq.m_left_repeat < kmer_len);
                connected = connected && (NextIBase() == NextJBase());
                if(connected) {
                    swap(icontigp, jcontigp);
                    swap(diri, dirj);
                    connected = (NextIBase() == NextJBase());
                }
                if(!connected)
                    overlaps.erase(it);
            }
            
            cerr << "Overlap connections: " << overlaps.size() << " " << kmer_len << endl;
                
            while(!overlaps.empty()) {
                auto& overlap = overlaps.front();
                auto icontigp = get<0>(overlap);
                auto jcontigp = get<1>(overlap);
                int diri = get<2>(overlap);
                int dirj = get<3>(overlap);
                int overlap_len = get<4>(overlap);
                
                if(diri > 0) {
                    if(dirj > 0)
                        jcontigp->ReverseComplement();
                    jcontigp->ClipLeft(overlap_len-kmer_len+1);  // AddToRight assumes kmer-1 overlap   
                    icontigp->AddToRight(*jcontigp);
                } else {
                    if(dirj < 0)
                        jcontigp->ReverseComplement();
                    jcontigp->ClipRight(overlap_len-kmer_len+1);  // AddToLeft assumes kmer-1 overlap   
                    icontigp->AddToLeft(*jcontigp);
                }
                overlaps.pop_front();

                for(auto& overlap : overlaps) {
                    if(get<0>(overlap) == jcontigp) {
                        get<0>(overlap) = icontigp;
                        get<2>(overlap) = diri;
                    } else if(get<1>(overlap) == jcontigp) {
                        get<1>(overlap) = icontigp;
                        get<3>(overlap) = diri;
                    }
                }

                scontigs.erase(jcontigp);
            }                
        }


        // Starting from available graph nodes, generates all contigs >= min_len_for_new_seeds. Uses ncores threads.
        TContigList<DBGraph> GenerateNewSeeds(int min_len_for_new_seeds, int ncores, CDBGraphDigger* test_graphdiggerp) {
            //assemble new seeds
            vector<TContigList<DBGraph>> new_seeds_for_threads(ncores);
            list<function<void()>> jobs;
            for(auto& ns : new_seeds_for_threads) {
                jobs.push_back(bind(&CDBGraphDigger::NewSeedsJob, this, ref(ns), min_len_for_new_seeds));
            }
            RunThreads(ncores, jobs);

            //connect fragments 
            Graph().ClearHoldings();
            TContigList<DBGraph> new_seeds = SContig<DBGraph>::ConnectFragments(new_seeds_for_threads, Graph());
            
            int kmer_len = Graph().KmerLen();
            CReadHolder removed_seq(false);
            for(auto iloop = new_seeds.begin(); iloop != new_seeds.end(); ) {
                auto ic = iloop++;
                if((int)ic->LenMin() < min_len_for_new_seeds+2*kmer_len) {
                    removed_seq.PushBack(ic->m_seq[0].front());
                    new_seeds.erase(ic);
                    continue;
                } 

                if(test_graphdiggerp != nullptr) {
                    CReadHolder rh(false);
                    rh.PushBack(ic->m_seq[0].front());
                    double abundance = 0;
                    int knum = 0;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(test_graphdiggerp->Graph().KmerLen()) ; ik != rh.kend(); ++ik, ++knum) {
                        auto node = test_graphdiggerp->Graph().GetNode(*ik);
                        abundance += test_graphdiggerp->Graph().Abundance(node);
                    }
                    if(abundance < knum*test_graphdiggerp->m_hist_min) {
                        removed_seq.PushBack(ic->m_seq[0].front());
                        new_seeds.erase(ic);
                        continue;
                    }
                }

                if(!ic->m_seq.m_circular) {
                    string left(ic->m_seq[0].front().begin(), ic->m_seq[0].front().begin()+2*kmer_len-1);
                    removed_seq.PushBack(left);
                    string right(ic->m_seq[0].front().end()-2*kmer_len+1, ic->m_seq[0].front().end());
                    removed_seq.PushBack(right);
                    ic->ClipLeft(kmer_len);
                    ic->ClipRight(kmer_len); 
                    ic->m_seq.m_left_repeat = kmer_len-1;
                    ic->m_seq.m_right_repeat = kmer_len-1;
                }
            }
            for(CReadHolder::kmer_iterator ik = removed_seq.kbegin(kmer_len) ; ik != removed_seq.kend(); ++ik)            
                Graph().ClearVisited(Graph().GetNode(*ik));            
        
            return new_seeds;
        }
        // Using a longer kmer generates connectors and extenders and improves previously assembled contigs
        // scontigs - contigs (input/output)
        // ncores - number of threads
        void ConnectAndExtendContigs(TContigList<DBGraph>& scontigs, int ncores) {
            vector<TContigList<DBGraph>> extensions_for_jobs(ncores);
            {
                for(auto& contig : scontigs)
                    contig.m_is_taken = 0;

                list<function<void()>> jobs;
                for(auto& ex : extensions_for_jobs) {
                    jobs.push_back(bind(&CDBGraphDigger::ExtendContigsJob, this, ref(scontigs), ref(ex)));
                }
                RunThreads(ncores, jobs);
            }

            TContigList<DBGraph> extensions = SContig<DBGraph>::ConnectFragments(extensions_for_jobs, Graph()); 
            SContig<DBGraph>::ConnectAndExtendContigs(scontigs, extensions, ncores); 

            //stabilize orientation which is random in multithreading 
            {
                for(auto& contig : scontigs)
                    contig.m_is_taken = 0;

                list<function<void()>> jobs;
                for(int thr = 0; thr < ncores; ++thr) {
                    jobs.push_back(bind(&CDBGraphDigger::StabilizeContigJob, this, ref(scontigs)));
                }
                RunThreads(ncores, jobs);
            }
        }
        
        list<array<CReadHolder,2>> ConnectPairs(const list<array<CReadHolder,2>>& mate_pairs, int insert_size, int ncores, bool extend_connected) {
            CStopWatch timer;
            timer.Restart();

            list<array<CReadHolder,2>> paired_reads;
            list<function<void()>> jobs;
            for(auto& reads : mate_pairs) {
                auto& job_input = reads[0];
                paired_reads.push_back(array<CReadHolder,2>({CReadHolder(false), CReadHolder(true)}));
                if(job_input.ReadNum() > 0)  // not empty       
                    jobs.push_back(bind(&CDBGraphDigger::ConnectPairsJob, this, insert_size, ref(job_input), ref(paired_reads.back()), extend_connected));            
            }
            RunThreads(ncores, jobs);
                       
            size_t connected = 0;
            size_t not_connected = 0;
            for(auto& rh : paired_reads) {
                connected += rh[0].ReadNum();
                not_connected += rh[1].ReadNum();
            }
            size_t mates = 0;
            for(auto& rh : mate_pairs)
                mates += rh[0].ReadNum();
            cerr << "Connected: " << connected << " ambiguously connected: " << not_connected/2 << " from " << mates/2 << " mate pairs" << endl;        
            cerr << "Connect pairs in " << timer.Elapsed(); 

            return paired_reads;
        }

        // remove read parts not supported by graph
        uint8_t CheckAndClipReadLite(string& read) {
            int kmer_len = m_graph.KmerLen();
            int rlen = read.size();
            if(rlen < kmer_len) {
                read.clear();
                return 0;
            }

            deque<Node> nodes;
            CReadHolder rh(false);
            rh.PushBack(read);
            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)   // iteration from last kmer to first  
                nodes.push_front(m_graph.GetNode(*ik));            

            vector<int> bases(read.size(), 0);
            for(int ek = 0; ek < (int)nodes.size(); ++ek) {
                Node node = nodes[ek];
                if(node.isValid() && GoodNode(node)) {
                    int left_kmer_end = ek;    // left kmer position on read
                    int right_kmer_end = left_kmer_end+kmer_len-1; // right kmer position on read
                    for(int p = left_kmer_end; p <= right_kmer_end; ++p)
                        bases[p] = 1;
                }
            }

            int left = 0;                                           // first good position    
            int len = 0;                                            // number of consecutive good positions    
            for(int k = 0; k < rlen; ++k) {
                for( ; k < rlen && !bases[k]; ++k);                 // skip bad bases    
                int current_left = k;
                int current_len = 0;
                for( ; k < rlen && bases[k]; ++k, ++current_len);   // count adjacent good bases 
                if(current_len > len) {
                    left = current_left;
                    len = current_len;
                }
            }


            uint8_t color = 0;
            if(len < kmer_len) {
                read.clear();
            } else {
                read = read.substr(left, len); 
                for(int ek = left; ek <= left+len-kmer_len; ++ek) {
                    auto& node = nodes[ek];
                    if(node.isValid())
                        color |= m_graph.GetColor(node);
                }
            }

            return color;
        } 
       
        // Prepares one of the mates of a read pair for connection
        // Finds the longest stretch of the read which could be assembled from both ends and clips the rest
        // read - read input/output
        // nodes - kmers for the remaining part
        // returns left and length for retained segment
        pair<int,int> CheckAndClipRead(string& read, deque<Node>& nodes) {
            int kmer_len = m_graph.KmerLen();

            string lextend = MostLikelyExtension(DBGraph::ReverseComplement(m_graph.GetNode(read.substr(0, kmer_len))), kmer_len);        
            ReverseComplementSeq(lextend.begin(), lextend.end());
            string rextend = MostLikelyExtension(m_graph.GetNode(read.substr(read.size()-kmer_len)), kmer_len);

            deque<Node> extended_nodes;
            CReadHolder rh(false);
            rh.PushBack(lextend+read+rextend);
            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first  
                extended_nodes.push_front(m_graph.GetNode(*ik));

            vector<int> bases(read.size(), 0);
            unsigned read_pos = kmer_len-lextend.size();
            for(int kk = 0; lextend.size()+read_pos+1 < extended_nodes.size() && read_pos < read.size(); ++kk, ++read_pos) {
                Node left = extended_nodes[kk];
                Node node = extended_nodes[kk+1];
                if(!left.isValid() || !GoodNode(left) || !node.isValid() || !GoodNode(node))
                    continue;
                vector<Successor> successors = m_graph.GetNodeSuccessors(left);
                FilterNeighbors(successors, false);
                if(find_if(successors.begin(),successors.end(),[node](const Successor& s){return s.m_node == node;}) == successors.end())                
                    continue;
            
                Node right = m_graph.ReverseComplement(extended_nodes[lextend.size()+read_pos+1]);
                node = m_graph.ReverseComplement(extended_nodes[read_pos+lextend.size()]);
                if(!right.isValid() || !GoodNode(right) || !node.isValid() || !GoodNode(node))
                    continue;
                successors = m_graph.GetNodeSuccessors(right);
                FilterNeighbors(successors, false);
                if(find_if(successors.begin(),successors.end(),[node](const Successor& s){return s.m_node == node;}) == successors.end())
                    continue;                

                bases[read_pos] = 1;
            }        

            int left = 0;             // first kmer position    
            int len = 0;              // number of consecutive good bases    
            for(unsigned k = 0; k < read.size(); ++k) {
                for( ; k < read.size() && !bases[k]; ++k);         // skip bad bases    
                int current_left = k;
                int current_len = 0;
                for( ; k < read.size() && bases[k]; ++k, ++current_len);   // count adjacent good bases 
                if(current_len > len) {
                    left = current_left;
                    len = current_len;
                }
            }

            if(len < kmer_len) {
                left = read.size();
                len = 0;
                read.clear();
                nodes.clear();
            } else {
                read = read.substr(left, len);
                nodes.resize(len-kmer_len+1);
                copy(extended_nodes.begin()+lextend.size()+left, extended_nodes.begin()+lextend.size()+left+len-kmer_len+1, nodes.begin());
            }

            return make_pair(left, len);
        }    
        
    private:

        // one-thread worker for paired reads connection
        // saves reads which were unambiguously connected; extends the ends of ambiguously connected reads and
        // keeps them for future; discards reads which don't have connection
        // insert_size - the maximal limit of the insert length
        // mate_pairs - pairs for connection (one mate after another)
        // paired_reads - [0] connected reads, [1] reads for future connection    
        void ConnectPairsJob(int insert_size, const CReadHolder& mate_pairs, array<CReadHolder,2>& paired_reads, bool extend_connected) {
            if(mate_pairs.ReadNum() < 2)
                return;

            int kmer_len = m_graph.KmerLen();

            for(CReadHolder::string_iterator is = mate_pairs.sbegin(); is != mate_pairs.send(); ++is) {
                string read1 = *is;
                string read2 = *(++is);
                if((int)min(read1.size(),read2.size()) < kmer_len)
                    continue;

                deque<Node> nodes1;
                CheckAndClipRead(read1, nodes1);
                if(read1.empty())
                    continue;
                Node last_node1 = nodes1.back();                        
                
                ReverseComplementSeq(read2.begin(), read2.end());
                deque<Node> nodes2;
                CheckAndClipRead(read2, nodes2);
                if(read2.empty())
                    continue;
                Node first_node2 = nodes2.front();

                int steps = insert_size;
                bool ambiguous = false;
                string read;

                //check for long overlap with extension     
                int hit = find(nodes2.begin(), nodes2.end(), last_node1) - nodes2.begin(); // first kmer position of the hit        
                if(hit < (int)min(nodes1.size(),nodes2.size()) && equal(nodes2.begin(), nodes2.begin()+hit, nodes1.end()-hit-1)) { // overlap
                    // check for circularity
                    pair<TBases<DBGraph>, EConnectionStatus> rslt = ConnectTwoNodes(last_node1, last_node1, steps);
                    if(rslt.second == CDBGraphDigger::eNoConnection)
                        read = read1+read2.substr(hit+kmer_len);
                    else
                        ambiguous = true;                                             
                } else {
                    pair<TBases<DBGraph>, EConnectionStatus> rslt = ConnectTwoNodes(last_node1, first_node2, steps);
                    if(rslt.second == CDBGraphDigger::eAmbiguousConnection) {
                        ambiguous = true; 
                    } else {
                        if(rslt.second == eSuccess) {
                            string r1 = read1;
                            for(auto& suc : rslt.first) {
                                r1.push_back(suc.m_nt);
                            }
                            r1 += read2.substr(kmer_len);
                            rslt = ConnectTwoNodes(DBGraph::ReverseComplement(first_node2), DBGraph::ReverseComplement(last_node1), steps);
                            if(rslt.second == eSuccess) {
                                string seq;
                                for(auto& suc : rslt.first)
                                    seq.push_back(suc.m_nt);
                                ReverseComplementSeq(seq.begin(), seq.end());
                                string r2 = read1.substr(0, read1.size()-kmer_len)+seq+read2;
                                if(r1 == r2)
                                    read = r1;
                            }
                    
                            if(read.empty())
                                ambiguous = true;
                        }
                    }
                }

                if(!read.empty()) { 
                    if(extend_connected) {
                        string lextend =  StringentExtension(DBGraph::ReverseComplement(nodes1.front()), kmer_len).first;
                        ReverseComplementSeq(lextend.begin(), lextend.end());
                        read = lextend+read;
                        read += StringentExtension(nodes2.back(), kmer_len).first;
                    }
                    paired_reads[0].PushBack(read);                      
                } else if(ambiguous) {
                    string lextend =  StringentExtension(DBGraph::ReverseComplement(nodes1.front()), kmer_len).first;
                    ReverseComplementSeq(lextend.begin(), lextend.end());
                    paired_reads[1].PushBack(lextend+read1);
                    read2 += StringentExtension(nodes2.back(), kmer_len).first;
                    ReverseComplementSeq(read2.begin(), read2.end());
                    paired_reads[1].PushBack(read2);
                }
            }                                
        }

        // one-thread worker for generating new seeds
        // returns contigs sequences which are either >= min_len or are known fragments
        // contigs - generated contigs
        // min_len - minimal length for acceptable contigs
        void NewSeedsJob(TContigList<DBGraph>& contigs, int min_len) {
            for(auto it = Graph().Begin(); it != Graph().End(); ++it) {
                SContig<DBGraph> contig = GetContigForKmer(it, min_len);
                if(!contig.m_seq.empty())
                    contigs.push_back(contig);
            }
        }

        void StabilizeContigJob(TContigList<DBGraph>& scontigs) {
            for(auto& contig : scontigs) {
                if(contig.m_is_taken.Set(1))  // grab contig
                    contig.SelectMinDirection();
            }
        }

        // one-thread worker for generating connectors and extenders for previously assembled contigs
        // scontigs - contigs (input/output)
        // extensions - generated sequences
        void ExtendContigsJob(TContigList<DBGraph>& scontigs, TContigList<DBGraph>& extensions) {
            for(auto& contig : scontigs) {
                if(contig.m_seq.m_circular || !contig.m_is_taken.Set(1))  // grab contig
                    continue;

                int kmer_len = contig.m_kmer_len;
                int chunks = contig.m_seq.size();
                
                if(contig.m_seq.m_right_repeat < kmer_len) {
                    Node takeoff_node = contig.BackKmer();
                    if(takeoff_node.isValid() && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer 
                        int allowed_intrusion = max(0, (int)contig.m_seq.ChunkLenMax(chunks-1)-kmer_len);
                        tuple<SContig<DBGraph>, Node, int> extension = ExtendToRight(takeoff_node, allowed_intrusion);
                        if(!get<0>(extension).m_seq.empty() || get<1>(extension).isValid()) { // extension could be empty - starting kmer + landing kmer  
                            bool skip = false;
                            int intrusion = get<2>(extension);
                            if(intrusion > 0 && !get<1>(extension).isValid()) { // there is intrusion and no further extension
                                CContigSequence& ext_seq = get<0>(extension).m_seq;
                                int ext_chunks = ext_seq.size();
                                int last_chunk = ext_seq.ChunkLenMax(ext_chunks-1);
                                if(last_chunk < kmer_len && (int)ext_seq.LenMax()-last_chunk-(int)ext_seq.ChunkLenMax(ext_chunks-2) < intrusion) // last chunk and snp will be clipped resulting in shorter sequence
                                    skip = true;
                                string back_kmer_seq(contig.m_seq.back().front().end()-intrusion-kmer_len, contig.m_seq.back().front().end()-intrusion);
                                TKmer back_kmer(back_kmer_seq);
                                if(!Graph().GetNode(back_kmer).isValid()) // new back kmer is not in the graph
                                    skip = true;
                            }

                            if(!skip) {
                                contig.ClipRight(intrusion);
                                SContig<DBGraph> sc(&contig, 1, contig.BackKmer(), get<0>(extension), get<1>(extension), Graph());
                                extensions.push_back(sc); 
                            }
                        }
                    }                
                }
                                                 
                if(contig.m_seq.m_left_repeat < kmer_len) {
                    Node takeoff_node = DBGraph::ReverseComplement(contig.FrontKmer());
                    if(takeoff_node.isValid() && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer     
                        int allowed_intrusion = max(0, (int)contig.m_seq.ChunkLenMax(0)-kmer_len);
                        tuple<SContig<DBGraph>, Node, int> extension = ExtendToRight(takeoff_node, allowed_intrusion);
                        if(!get<0>(extension).m_seq.empty() || get<1>(extension).isValid()) { // extension could be empty - starting kmer + landing kmer   
                            bool skip = false;
                            int intrusion = get<2>(extension);
                            if(intrusion > 0 && !get<1>(extension).isValid()) { // there is intrusion and no further extension
                                CContigSequence& ext_seq = get<0>(extension).m_seq;
                                int ext_chunks = ext_seq.size();
                                int last_chunk = ext_seq.ChunkLenMax(ext_chunks-1);
                                if(last_chunk < kmer_len && (int)ext_seq.LenMax()-last_chunk-(int)ext_seq.ChunkLenMax(ext_chunks-2) < intrusion) // last chunk and snp will be clipped resulting in shorter sequence
                                    skip = true;
                                string front_kmer_seq(contig.m_seq.front().front().begin()+intrusion, contig.m_seq.front().front().begin()+intrusion+kmer_len);
                                TKmer front_kmer(front_kmer_seq);
                                if(!Graph().GetNode(front_kmer).isValid()) // new back kmer is not in the graph
                                    skip = true;
                            }

                            if(!skip) {
                                contig.ClipLeft(intrusion);
                                SContig<DBGraph> sc(&contig, -1, DBGraph::ReverseComplement(contig.FrontKmer()), get<0>(extension), get<1>(extension), Graph());
                                sc.ReverseComplement();
                                extensions.push_back(sc); 
                            }
                        }
                    }
                }
            }
        }        


        DBGraph& m_graph;
        double m_fraction;
        int m_jump;
        int m_hist_min;
        int m_low_count;
        size_t m_max_branch;
        bool m_allow_snps;
    };

} // namespace
#endif /* _GraphDigger_ */

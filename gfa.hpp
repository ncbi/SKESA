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

#ifndef _GFA_
#define _GFA_

#include "DBGraph.hpp"
#include "graphdigger.hpp"
#include "genetic_code.hpp"
#include <random>

using namespace std;
namespace DeBruijn {
    typedef CDBHashGraph DBGraph;
    typedef DBGraph::Node Node;
    typedef DBGraph::Successor Successor;
    typedef CDBGraphDigger<DBGraph> GraphDigger;

    // unlimited counter
    struct SVarNum {
        SVarNum(uint32_t n = 0) : m_data(1, n) {}
        SVarNum& operator+=(const SVarNum& other) {
            uint64_t overflow = 0;
            for(unsigned i = 0; i < other.m_data.size() || overflow > 0; ++i) {
                if(i == m_data.size())
                    m_data.push_back(0);
                overflow += m_data[i];
                if(i < other.m_data.size())
                    overflow += other.m_data[i];
                m_data[i] = overflow;
                overflow >>= 32;
            }
            return *this;
        }
        bool operator<(const SVarNum& other) const {
            int n = m_data.size();
            while(m_data[n-1] == 0 && n > 1)
                --n;
            int othern = other.m_data.size();
            while(other.m_data[othern-1] == 0 && othern > 1)
                --othern;
            if(n != othern) {
                return n < othern;
            } else {
                for(int i = n-1; i >= 0; --i) {
                    if(m_data[i] != other.m_data[i])
                        return m_data[i] < other.m_data[i];
                }
                return false;
            }
        }
        string ToString() const {
            // double dabble
            int nbits = 32*m_data.size();
            int nbcd = nbits/3;
            int smin = nbcd-2;
            vector<uint8_t> bcd(nbcd);

            for(int i = m_data.size()-1; i >= 0; --i) {
                for(int j = 0; j < 32; ++j) {
                    for(int k = smin; k < nbcd; ++k)
                        bcd[k] += (bcd[k] >= 5) ? 3 : 0;
                    if(bcd[smin] >= 8)
                        smin -= 1;

                    for(int k = smin; k < nbcd-1; ++k) {
                        bcd[k] <<= 1;
                        bcd[k] &= 0xF;
                        bcd[k] |= (bcd[k+1] >= 8);
                    }

                    int shifted_in = (m_data[i] & (1 << (31-j))) ?  1 : 0;
                    bcd[nbcd-1] <<= 1;
                    bcd[nbcd-1] &= 0xF;
                    bcd[nbcd-1] |= shifted_in;                    
                }
            }

            auto first = find_if(bcd.begin(), bcd.end(), [](uint8_t c) { return c != 0; });
            if(first == bcd.end())
                return "0";

            int prec = 6;
            int num = bcd.end()-first;
            string txt;
            if(num <= prec) {
                for(auto p = first; p != bcd.end(); ++p)
                    txt.push_back(*p+'0');
            } else {
                txt.push_back(*first+'0');
                txt.push_back('.');
                for(auto p = first+1; p < first+prec; ++p)
                    txt.push_back(*p+'0');
                txt += "e+"+to_string(num-1);
            }

            return txt;
        }

        vector<uint32_t> m_data;
    };

    struct SegBase {
        forward_list<Node> m_left_kmers;  // kmers which left ends are at this position 
        forward_list<Node> m_right_kmers; // kmer(s) which right ends are at this position
        int m_fork = eNoFork;
        char m_nt;
        bool operator==(const SegBase& other) const { return m_nt == other.m_nt; }
        bool operator!=(const SegBase& other) const { return m_nt != other.m_nt; }
        bool operator<(const SegBase& other) const { return m_nt < other.m_nt; }
    };

    class SegSeq : public deque<SegBase> {
    public:
        SegSeq() {}
        template <typename Base>
        SegSeq(const vector<Base>& seq) {
            for(auto& base : seq) {
                emplace_back();
                back().m_nt = base.m_nt;
            }
        }
        void RightExtend(char c, const Node& node) {
            SegBase base;
            base.m_nt = c;
            base.m_right_kmers.push_front(node);
            push_back(base);
        }
        void LeftExtend(char c, const Node& node) {
            SegBase base;
            base.m_nt = c;
            base.m_left_kmers.push_front(node);
            push_front(base);
        }
        void ClipLeft(int len) { erase(begin(), begin()+len); }
        void ClipRight(int len) { erase(end()-len, end()); }
        void SplitForksAndKmersFrom(SegSeq& other) {
            if(size() != other.size())
                throw runtime_error("Can't SplitForksAndKmersFrom from different size SegSeq");
            for(size_t p = 0; p < size(); ++p) {
                auto& base = (*this)[p];
                auto& other_base = other[p];
                base.m_fork |= other_base.m_fork;
                base.m_left_kmers.splice_after(base.m_left_kmers.before_begin(), other_base.m_left_kmers);
                base.m_left_kmers.sort();
                base.m_left_kmers.unique();
                base.m_right_kmers.splice_after(base.m_right_kmers.before_begin(), other_base.m_right_kmers);
                base.m_right_kmers.sort();
                base.m_right_kmers.unique();
            }

        }
        void ReverseComplement() {
            reverse(this->begin(), this->end());
            for(auto& base : *this) {
                base.m_nt = Complement(base.m_nt);
                base.m_left_kmers.swap(base.m_right_kmers);
                for(Node& node : base.m_left_kmers)  
                    node = DBGraph::ReverseComplement(node);
                for(Node& node : base.m_right_kmers)                        
                    node = DBGraph::ReverseComplement(node);
            }
        }
        SegSeq substr(int pos, int len = numeric_limits<int>::max()) const { 
            SegSeq x;
            x.insert(x.end(), this->begin()+pos, this->begin()+pos+min(len,(int)this->size()-pos));
            return x;
        }
        SegSeq& operator+=(const SegSeq& other) {
            insert(this->end(), other.begin(), other.end());
            return *this;
        }
        SegSeq operator+(const SegSeq& other) const {
            SegSeq x = *this;
            x += other;
            return x;
        }
    };


    struct GFASegment;
    typedef typename list<GFASegment>::iterator GFAIterator;
    struct SGFAIteratorHash { size_t operator()(const GFAIterator& i) const { return std::hash<void*>()(&*i); }};
    typedef unordered_set<GFAIterator, SGFAIteratorHash> GFAIteratorUSet;

    struct GFASegment {
        GFASegment(const SegSeq& seq = SegSeq()) : m_seq(seq) {}
        bool operator==(const GFASegment& other) const { return m_seq == other.m_seq; }// lexigraphical order, nothing else 
        bool operator!=(const GFASegment& other) const { return m_seq != other.m_seq; }// lexigraphical order, nothing else 
        bool operator<(const GFASegment& other) const { return m_seq < other.m_seq; }  // lexigraphical order, nothing else 
        //during graph cleaning paths are aften copied; this hash is going to be used for identifying paths consisting of the 'same' segments
        size_t Hash() const { return std::hash<const void*>()(m_copy_of == nullptr ? this : m_copy_of); }
        bool LeftSingle() const { return !m_left_connections.empty() && next(m_left_connections.begin()) == m_left_connections.end(); }
        bool LeftFork() const { return !m_left_connections.empty() && next(m_left_connections.begin()) != m_left_connections.end(); }
        int LeftConnectionsNum() const { return distance(m_left_connections.begin(), m_left_connections.end()); }
        bool RightSingle() const { return !m_right_connections.empty() && next(m_right_connections.begin()) == m_right_connections.end(); }
        bool RightFork() const { return !m_right_connections.empty() && next(m_right_connections.begin()) != m_right_connections.end(); }
        int RightConnectionsNum() const { return distance(m_right_connections.begin(), m_right_connections.end()); }
        void ReverseComplement() {
            m_seq.ReverseComplement();
            swap(m_left_connections, m_right_connections);
            swap(m_left_kmer_count, m_right_kmer_count);
            swap(m_left_len, m_right_len);
            swap(m_left_check, m_right_check);
        }

        SegSeq m_seq;
        forward_list<GFAIterator> m_left_connections;
        forward_list<GFAIterator> m_right_connections;
        GFASegment* m_copy_of = nullptr;
        GFAIterator m_copy_ofi;
        size_t m_kmer_count = 0;
        size_t m_left_kmer_count = 0;
        size_t m_right_kmer_count = 0;
        int m_num = 0;
        int m_group = 0;
        int m_left_len = 0; // reused for not aligned and extension length      
        int m_right_len = 0; 
        bool m_left_check = false;
        bool m_right_check = false;
        bool m_marked_for_erase = false;
        bool m_cyclical = false;
        int m_frame = -1; // frame of the left end (0-based position%3); used only for protein targets in which case it is >= 0
    };

    struct Position {
        Position() {}
        Position(GFAIterator segmp, int pos) : m_segmp(segmp), m_pos(pos) {}
        GFAIterator m_segmp;
        int m_pos;
    };
    struct Path {
        //paths including differnt copies of the same sements are equal
        bool operator==(const Path& other) const {
            if(m_len != other.m_len || m_left != other.m_left || m_right != other.m_right || m_segments.size() != other.m_segments.size())
                return false;

            for(unsigned i = 0; i < m_segments.size(); ++i) {
                auto ptr = m_segments[i]->m_copy_of == nullptr ? &(*m_segments[i]) : m_segments[i]->m_copy_of;
                auto other_ptr = other.m_segments[i]->m_copy_of == nullptr ? &(*other.m_segments[i]) : other.m_segments[i]->m_copy_of;
                if(ptr != other_ptr)
                    return false;
            }

            return true;
        }
        struct Hash {
            size_t operator()(const Path& path) const {
                size_t h = 0;
                for(auto segi : path.m_segments)
                    h ^= segi->Hash();
                return h;
            }
        };
        const Position* CurrentPosition() const {
            if(m_segments.empty() || (m_current_seg == 0 && m_current_pos.m_pos < m_left) ||
               (m_current_seg == (int)m_segments.size()-1 && m_current_pos.m_pos > m_right)) {
                return nullptr;
            } else {
                return &m_current_pos;
            }
        }
        const Position* StepRight() {
            if(m_segments.empty() || CurrentPosition() == nullptr)
                return nullptr;
            if(++m_current_pos.m_pos == (int)m_current_pos.m_segmp->m_seq.size() && m_current_seg < (int)m_segments.size()-1) {
                m_current_pos.m_segmp = m_segments[++m_current_seg];
                m_current_pos.m_pos = 0;
            }
            return CurrentPosition();
        }
        const Position* StepLeft() {
            if(m_segments.empty() || CurrentPosition() == nullptr)
                return nullptr;
            if(--m_current_pos.m_pos < 0 && m_current_seg > 0) {
                m_current_pos.m_segmp = m_segments[--m_current_seg];
                m_current_pos.m_pos = m_current_pos.m_segmp->m_seq.size()-1;
            }
            return CurrentPosition();
        }
        const Position* JumpToRightEnd() {
            if(m_segments.empty())
                return nullptr;
            m_current_seg = m_segments.size()-1;
            m_current_pos.m_segmp = m_segments.back();
            m_current_pos.m_pos = m_right;
            return CurrentPosition();
        }
        const Position* JumpToLeftEnd() {
            if(m_segments.empty())
                return nullptr;
            m_current_seg = 0;
            m_current_pos.m_segmp = m_segments.front();
            m_current_pos.m_pos = m_left;
            return CurrentPosition();
        }
        string Sequence() const {
            string seq;
            for(unsigned i = 0; i < m_segments.size(); ++i) {
                int left = (i == 0) ? m_left : 0;
                int right = (i == m_segments.size()-1) ? m_right : m_segments[i]->m_seq.size()-1;
                for(int j = left; j <= right; ++j)
                    seq.push_back(m_segments[i]->m_seq[j].m_nt);
            }
            return seq;
        }
        int Length() const { return m_len; }
        Position LeftEnd() const {
            Position l;
            l.m_segmp = m_segments.front();
            l.m_pos = m_left;
            return l;
        }
        Position RightEnd() const {
            Position r;
            r.m_segmp = m_segments.back();
            r.m_pos = m_right;
            return r;
        }
        void ClipRight(int clip) {
            m_len -= clip;
            while(clip > 0) {
                if(m_right+1 <= clip) {
                    clip -= m_right+1;
                    m_segments.pop_back();
                    m_right = m_segments.back()->m_seq.size()-1;
                    if(m_current_seg == (int)m_segments.size()) {
                        --m_current_seg;
                        m_current_pos.m_segmp = m_segments.back();
                        m_current_pos.m_pos = m_right;
                    }
                } else {
                    m_right -= clip;
                    clip = 0;
                    if(m_current_seg == (int)m_segments.size()-1)
                        m_current_pos.m_pos = min(m_right, m_current_pos.m_pos);
                }
            }
        }
        void ClipLeft(int clip) {
            m_len -= clip;
            while(clip > 0) {
                int slen = m_segments.front()->m_seq.size()-m_left;
                if(slen <= clip) {
                    clip -= slen;
                    m_segments.pop_front();
                    m_left = 0;
                    if(m_current_seg == 0) {
                        m_current_pos.m_segmp = m_segments.front();
                        m_current_pos.m_pos = 0;
                    } else {
                        --m_current_seg;
                    }
                } else {
                    m_left += clip;
                    clip = 0;
                    if(m_current_seg == 0)
                        m_current_pos.m_pos = max(m_left, m_current_pos.m_pos);
                }
            }            
        }
        bool IntactPath() const {
            for(auto it = m_segments.begin(); it != prev(m_segments.end()); ++it) {
                if(find((*it)->m_right_connections.begin(), (*it)->m_right_connections.end(), *next(it)) == (*it)->m_right_connections.end())
                    return false;
            }
            return true;
        }

        deque<GFAIterator> m_segments;                      // all segments included in path (deque because we want a simple copy constructor)
        int m_left;                                      // starting position in leftmost segment
        int m_right;                                     // ending position in rigthmost segment (included)

        Position m_current_pos;                          // current point
        int m_current_seg;                               // segment for current point
        int m_len;                                     
    };


    typedef deque<uint64_t> TReadPosInfo;
    typedef CKmerHashMap<tuple<TReadPosInfo, TReadPosInfo, SAtomic<uint8_t>>, 8> TReadPos;
    typedef unordered_map<const GFAIterator, tuple<list<GFAIterator>,list<GFAIterator>>, SGFAIteratorHash> TCopyInfo;

    void RecalculateCopyInfo(TCopyInfo& copies, const GFAIteratorUSet& erased) {
        for(auto icopy = copies.begin(); icopy != copies.end(); ) {
            auto& checked = get<0>(icopy->second);
            auto& not_checked = get<1>(icopy->second);
            checked.erase(remove_if(checked.begin(), checked.end(), [&](GFAIterator p) { return erased.count(p); }), checked.end());                  // remove erased from list
            not_checked.erase(remove_if(not_checked.begin(), not_checked.end(), [&](GFAIterator p) { return erased.count(p); }), not_checked.end());  // remove erased from list
            if(checked.empty() && not_checked.empty()) {      // remove empty entry
                icopy = copies.erase(icopy); 
            } else if(erased.count(icopy->first)) {           // remove erased from keys 
                GFAIterator firstp;
                if(checked.empty()) {
                    firstp = not_checked.front();
                    not_checked.pop_front();
                } else {
                    firstp = checked.front();
                    checked.pop_front();
                }
                firstp->m_copy_of = nullptr;
                for(GFAIterator p : checked) {
                    p->m_copy_of = &(*firstp);
                    p->m_copy_ofi = firstp;
                }
                for(GFAIterator p : not_checked) {
                    p->m_copy_of = &(*firstp);
                    p->m_copy_ofi = firstp;
                }
                if(!checked.empty() || !not_checked.empty()) {
                    if(copies.size() == copies.max_load_factor()*copies.bucket_count())
                        throw runtime_error("Unexpected rehash");
                    auto& rslt = copies[firstp];
                    get<0>(rslt).splice(get<0>(rslt).end(), checked);
                    get<1>(rslt).splice(get<1>(rslt).end(), not_checked);
                }
                icopy = copies.erase(icopy);
            } else {
                ++icopy;
            }
        }
    }

    class GFAGraph;
    typedef list<GFAGraph> TGFACollection;
    typedef map<string, tuple<string, int, unordered_set<uint32_t>, SAtomic<uint8_t>>> TTargets;  // [accession], seq, min_len, words, centinel

    class GFAGraph : public list<GFASegment> {            
    protected:
        string m_acc;
        deque<CKmerHashCount::Index> m_ksignature;
        size_t m_size = 0;
        int m_score = 0;
        int m_kmer_len;
        int m_max_num = 0;
        int m_id = 0;
        SAtomic<uint8_t> m_sentinel = 0;
        bool m_is_aa = false;

    public:
        GFAGraph(const string& acc, int kmer_len) : m_acc(acc), m_kmer_len(kmer_len) {}
        void ReverseComplement() {
            for(auto& seg : *this)
                seg.ReverseComplement();
        }
        deque<CKmerHashCount::Index>& KSignature() { return m_ksignature; }
        const deque<CKmerHashCount::Index>& KSignature() const { return m_ksignature; }
        int& Score() { return m_score; }
        int Score() const { return m_score; }
        SAtomic<uint8_t>& Sentinel() { return m_sentinel; }
        int& ID() { return m_id; }
        size_t& Size() { return m_size; }
        int MaxNum() const { return m_max_num; }

        GFAGraph& operator=(const GFAGraph& other) {
            if(this != &other) {
                clear();
                map<const GFASegment*, GFAIterator> other_to_copy;
                for(auto& seg : other) {
                    push_back(seg);
                    other_to_copy[&seg] = prev(end());
                }
                m_acc = other.m_acc;
                m_ksignature = other.m_ksignature;
                m_size = other.m_size;
                m_score = other.m_score;
                m_kmer_len = other.m_kmer_len;
                m_max_num = other.m_max_num;
                m_id = other.m_id;
                m_sentinel = other.m_sentinel;
        
                for(auto& seg : *this) {
                    seg.m_num = ++m_max_num;
                    for(auto& lc : seg.m_left_connections)
                        lc = other_to_copy[&(*lc)];
                    for(auto& rc : seg.m_right_connections)
                        rc = other_to_copy[&(*rc)];
                    if(seg.m_copy_of != nullptr) {
                        seg.m_copy_ofi = other_to_copy[seg.m_copy_of];
                        seg.m_copy_of = &(*seg.m_copy_ofi);
                    }
                }        
            }

            return *this;
        }

        GFAGraph(const GFAGraph& other) {
            *this = other;
        }
    
        //checks only ksignature
        bool IsSubGraphOf(const GFAGraph& other) const {
            if(KSignature().size() > other.KSignature().size())
                return false;

            deque<CKmerHashCount::Index> intersect(KSignature().size());
            auto iend = set_intersection (KSignature().begin(), KSignature().end(),
                                          other.KSignature().begin(), other.KSignature().end(),
                                          intersect.begin());
            return (iend == intersect.end());
        }

        void CutToChunks() {
            for(auto& segm : *this) {
                if(segm.m_seq.back().m_fork & eRightFork) {
                    for(auto rc : segm.m_right_connections) 
                        rc->m_seq.front().m_fork |= eRightBranch;                
                }
                if(segm.m_seq.front().m_fork & eLeftFork) {
                    for(auto lc : segm.m_left_connections)
                        lc->m_seq.back().m_fork |= eLeftBranch;
                }
            }
            for(auto it = begin(); it != end(); ++it) {
                auto& segm = *it;
                auto& seq = segm.m_seq;
                for(auto ib_loop = seq.begin(); ib_loop != seq.end(); ) {
                    auto ib = ib_loop++;
                    if(ib != seq.begin() && (ib->m_fork & eLeftFork)) { // not first and left fork  
                        int len = ib-seq.begin();
                        push_front(seq.substr(0, len));
                        ++m_size;
                        seq.erase(seq.begin(), seq.begin()+len);
                        front().m_seq.back().m_fork |= eLeftBranch;
                        front().m_group = segm.m_group;
                        front().m_cyclical = segm.m_cyclical;
                        front().m_num = ++m_max_num;
                        if(segm.m_frame >= 0) {
                            front().m_frame = segm.m_frame;
                            segm.m_frame = (segm.m_frame+len)%3;
                        }
                        TransferLeftLinks(it, begin());
                        LinkSegments(begin(), it);
                    }
                    if(ib_loop != seq.end() && (ib->m_fork & eRightFork)) { // not last and right fork  
                        int len = ib-seq.begin()+1;
                        push_front(seq.substr(0, len));
                        ++m_size;
                        seq.erase(seq.begin(), seq.begin()+len);
                        seq.front().m_fork |= eRightBranch;
                        front().m_group = segm.m_group;
                        front().m_cyclical = segm.m_cyclical;
                        front().m_num = ++m_max_num;
                        if(segm.m_frame >= 0) {
                            front().m_frame = segm.m_frame;
                            segm.m_frame = (segm.m_frame+len)%3;
                        }
                        TransferLeftLinks(it, begin());
                        LinkSegments(begin(), it);
                    }
                }
            }                 
        }

        void RemovePath(Path& path, bool toright, TCopyInfo& copies) {
            int segn = path.m_segments.size();
            GFAIteratorUSet erased;
            SplitInletsBeforeClip(path, toright, copies);
            if(toright) {
                int clip = path.m_right+1;
                for(int j = segn-2; j > 0; --j) {
                    auto& segi = path.m_segments[j];
                    auto& connections = segi->m_right_connections;
                    if(distance(connections.begin(), connections.end()) < 2) {
                        erased.insert(segi);
                        clip += segi->m_seq.size();
                    } else {
                        break;
                    }
                }  

                UnLinkSegments(path.m_segments[segn-2], path.m_segments[segn-1]);
                path.ClipRight(clip);
            } else {
                int clip = path.m_segments.front()->m_seq.size()-path.m_left;
                for(int j = 1; j < segn-1; ++j) {
                    auto& segi = path.m_segments[j];
                    auto& connections = segi->m_left_connections;
                    if(distance(connections.begin(), connections.end()) < 2) {
                        erased.insert(segi);
                        clip += segi->m_seq.size();
                    } else {
                        break;
                    }
                }

                UnLinkSegments(path.m_segments[0], path.m_segments[1]);
                path.ClipLeft(clip);
            }
        
            for(auto segi : erased) {
                if(!segi->m_marked_for_erase) {
                    RemoveLinksToSegment(segi);
                    segi->m_left_connections.clear();
                    segi->m_right_connections.clear();
                    segi->m_marked_for_erase = true;
                }
            }        
        }

        //path includes first base after clip
        void SplitInletsBeforeClip(Path& path, bool toright, TCopyInfo& copies) {
            int segn = path.m_segments.size();
            unordered_map<GFAIterator, GFAIterator, SGFAIteratorHash> copy_info; // iterator to orig -> iterator to copy  
            if(toright) {
                int left_inlet = find_if(path.m_segments.begin()+1, path.m_segments.end()-1, [](GFAIterator i) {return i->LeftFork();}) - path.m_segments.begin(); //ends not included
                if(left_inlet >= segn-1)
                    return;

                // copy segments; remember iterator; change num of copy (main graph)        
                for(int i = left_inlet; i < segn-1; ++i) {
                    auto iorig = path.m_segments[i];
                    copy_info[iorig] = PushSegmentBack(*iorig);
                    if(back().m_copy_of == nullptr) {
                        back().m_copy_of = &(*iorig);
                        back().m_copy_ofi = iorig;
                    }
                    if(back().m_left_check && back().m_right_check)
                        get<1>(copies[back().m_copy_ofi]).push_back(prev(end()));
                    else
                        get<0>(copies[back().m_copy_ofi]).push_back(prev(end()));
                }

                //remove left connection to copy subpath
                copy_info[path.m_segments[left_inlet]]->m_left_connections.remove(path.m_segments[left_inlet-1]);

                //fix connections
                for(int i = left_inlet; i < segn-1; ++i) {
                    auto iorig = path.m_segments[i];
                    auto iprev_orig = path.m_segments[i-1];
                    auto icpy = copy_info[iorig];

                    //strip orig from all left connections not directly in path
                    for(auto& lc : icpy->m_left_connections) {
                        if(lc != iprev_orig)
                            UnLinkSegments(lc, iorig);                             // unlink orig if not directly in path                            
                        auto rslt = copy_info.find(lc);
                        if(rslt == copy_info.end()) {                              //from outside
                            lc->m_right_connections.push_front(icpy);              // link copy
                        } else {                                                   //from copied
                            if(lc != iprev_orig) {
                                lc->m_right_connections.push_front(icpy);          // relink orig to copy  
                                icpy->m_left_connections.push_front(rslt->second); // create same link in copied path (this copied segment is now connected to both paths)
                            } else {
                                lc = rslt->second;                                 // link copy to previous copied segment
                            }
                        }
                    }

                    //duplicate all right connections in copy
                    for(auto& rc : icpy->m_right_connections) {
                        auto rslt = copy_info.find(rc);
                        if(rslt != copy_info.end())                    // link to copied segment   
                            rc = rslt->second;
                        else                                        
                            rc->m_left_connections.push_front(icpy);   // connect to copy                        
                    }
                }
            } else {
                int right_inlet = find_if(path.m_segments.rbegin()+1, path.m_segments.rend()-1, [](GFAIterator i) {return i->RightFork();}) - path.m_segments.rbegin(); //ends not included
                right_inlet = segn-1-right_inlet; //in left to right direction
                if(right_inlet < 1)
                    return;

                // copy segments; remember iterator; change num of copy (main graph)        
                for(int i = 1; i <= right_inlet; ++i) {
                    auto iorig = path.m_segments[i];
                    copy_info[iorig] = PushSegmentBack(*iorig);
                    if(back().m_copy_of == nullptr) {
                        back().m_copy_of = &(*iorig);
                        back().m_copy_ofi = iorig;
                    }
                    if(back().m_left_check && back().m_right_check)
                        get<1>(copies[back().m_copy_ofi]).push_back(prev(end()));
                    else
                        get<0>(copies[back().m_copy_ofi]).push_back(prev(end()));
                }
                

                //remove right connection to copy subpath
                copy_info[path.m_segments[right_inlet]]->m_right_connections.remove(path.m_segments[right_inlet+1]);

                //fix connections
                for(int i = 1; i <= right_inlet; ++i) {
                    auto iorig = path.m_segments[i];
                    auto inext_orig = path.m_segments[i+1];
                    auto icpy = copy_info[iorig];

                    //strip orig from all right connections not in path
                    for(auto& rc : icpy->m_right_connections) {
                        if(rc != inext_orig)
                            UnLinkSegments(iorig, rc);                // unlink orig if not in path 
                        auto rslt = copy_info.find(rc);
                        if(rslt == copy_info.end()) {
                            rc->m_left_connections.push_front(icpy);
                        } else {
                            if(rc != inext_orig){
                                rc->m_left_connections.push_front(icpy);
                                icpy->m_right_connections.push_front(rslt->second);
                            } else {
                                rc = rslt->second;
                            }
                        }

                        /*
                        if(rc != inext_orig)
                            UnLinkSegments(iorig, rc);                // unlink orig if not in path 
                        auto rslt = copy_info.find(rc);
                        if(rslt != copy_info.end())                   // link to copied segment    
                            rc = rslt->second;
                        else                                          // transfer right link to copy       
                            rc->m_left_connections.push_front(icpy);
                        */                        
                    }

                    //duplicate all left connections in copy
                    for(auto& lc : icpy->m_left_connections) {
                        auto rslt = copy_info.find(lc);
                        if(rslt != copy_info.end())                    // link to copied segment   
                            lc = rslt->second;
                        else                                        
                            lc->m_right_connections.push_front(icpy);   // connect to copy                        
                    }
                }     
            }
        }

        void GenerateKmers(DBGraph& dbg_graph) {
            m_ksignature.clear();
            int kmer_len = dbg_graph.KmerLen();

            for(auto it = begin(); it != end(); ++it) {
                auto& segm = *it;
                auto& seq = segm.m_seq;
                int seq_len = seq.size();
                string s;
                for(auto& base : seq) {
                    base.m_left_kmers.clear();
                    base.m_right_kmers.clear();
                    s.push_back(base.m_nt);
                }

                if(seq_len >= kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(s);
                    int lpos = seq_len-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len); ik != rh.kend(); ++ik, --lpos) { // iteration from last kmer to first    
                        auto kmer = dbg_graph.GetNode(*ik);
                        if(kmer.isValid()) {
                            auto& lk = seq[lpos].m_left_kmers;
                            auto& rk = seq[lpos+kmer_len-1].m_right_kmers;
                            lk.push_front(kmer);
                            rk.push_front(kmer);
                            m_ksignature.push_back(kmer);                            
                        }
                    }
                }

                if(!segm.m_left_connections.empty()) {
                    int l = min(kmer_len-1, seq_len);
                    list<Path> lpaths = Expand(it, 0, kmer_len-1, l-1);
                    for(auto& path : lpaths) {
                        if(path.Length() >= kmer_len) {
                            CReadHolder rh(false);
                            rh.PushBack(path.Sequence());
                            int rpos = l-1;
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len); ik != rh.kend(); ++ik, --rpos) { // iteration from last kmer to first    
                                auto kmer = dbg_graph.GetNode(*ik);
                                if(kmer.isValid()) {
                                    auto& rk = seq[rpos].m_right_kmers;
                                    if(find(rk.begin(), rk.end(), kmer) == rk.end()) {
                                        rk.push_front(kmer);
                                        m_ksignature.push_back(kmer); 
                                    }                       
                                }
                            }
                        }
                    }
                }

                if(!segm.m_right_connections.empty()) {
                    int l = min(kmer_len-1, seq_len);
                    list<Path> rpaths = Expand(it, seq_len-1, l-1, kmer_len-1);
                    for(auto& path : rpaths) {
                        if(path.Length() >= kmer_len) {
                            CReadHolder rh(false);
                            rh.PushBack(path.Sequence());
                            int lpos = seq_len-l+path.Length()-kmer_len;
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len); ik != rh.kend(); ++ik, --lpos) { // iteration from last kmer to first    
                                auto kmer = dbg_graph.GetNode(*ik);
                                if(kmer.isValid()) {
                                    auto& lk = seq[lpos].m_left_kmers;
                                    if(find(lk.begin(), lk.end(), kmer) == lk.end()) {
                                        lk.push_front(kmer);
                                        m_ksignature.push_back(kmer); 
                                    }
                                }                       
                            }
                        }
                    }                
                }
                
            }

            std::sort(m_ksignature.begin(), m_ksignature.end());
            m_ksignature.erase(std::unique(m_ksignature.begin(),m_ksignature.end()), m_ksignature.end());

            //calculate coverage;   
            for(auto& seg : *this) {
                seg.m_kmer_count = 0;
                for(auto& base : seg.m_seq) {
                    size_t lcount = 0;
                    size_t rcount = 0;
                    for(Node& node : base.m_left_kmers)
                        lcount += dbg_graph.Abundance(node);
                    for(Node& node : base.m_right_kmers)
                        rcount += dbg_graph.Abundance(node);
                    seg.m_kmer_count += max(lcount, rcount);
                }
            }            
                                              
            EnumerateSegments();
         }
        
        void GenerateKmersAndScores(DBGraph& dbg_graph) {
            GenerateKmers(dbg_graph);
            CalculateChainLength();
            CalculateCoverageAndEnumerateSegments();
        }

        Path ExpandEdgeToMax(GFAIterator left, GFAIterator right) {
            Path path;
            path.m_segments.push_back(left);
            path.m_segments.push_back(right);
            int len = left->m_seq.size()+right->m_seq.size();

            //expand to right   
            while(!path.m_segments.back()->m_right_connections.empty()) {
                auto& last = path.m_segments.back();
                GFAIterator imax = last->m_right_connections.front();
                int maxlen = imax->m_right_len+imax->m_seq.size();
                size_t maxkmer = imax->m_right_kmer_count+imax->m_kmer_count;
                for(auto it = next(last->m_right_connections.begin()); it != last->m_right_connections.end(); ++it) {
                    GFAIterator i = *it;
                    int ilen = i->m_right_len+i->m_seq.size();
                    size_t ikmer = i->m_right_kmer_count+i->m_kmer_count;
                    if(ikmer < maxkmer) {
                        continue;
                    } else if(ikmer > maxkmer) {        //maxcount  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(ilen < maxlen) {
                        continue;
                    } else if(ilen > maxlen) {          //longest   
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(i->m_seq < imax->m_seq) { //alphabetical  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    }
                }
                path.m_segments.push_back(imax);
                len += imax->m_seq.size();
            }

            //expand to left    
            while(!path.m_segments.front()->m_left_connections.empty()) {
                auto& first = path.m_segments.front();
                GFAIterator imax = first->m_left_connections.front();
                int maxlen = imax->m_left_len+imax->m_seq.size();
                size_t maxkmer = imax->m_left_kmer_count+imax->m_kmer_count;
                for(auto it = next(first->m_left_connections.begin()); it != first->m_left_connections.end(); ++it) {
                    GFAIterator i = *it;
                    int ilen = i->m_left_len+i->m_seq.size();
                    size_t ikmer = i->m_left_kmer_count+i->m_kmer_count;
                    if(ikmer < maxkmer) {
                        continue;
                    } else if(ikmer > maxkmer) {        //maxcount  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(ilen < maxlen) {
                        continue;
                    } else if(ilen > maxlen) {          //longest   
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(i->m_seq < imax->m_seq) { //alphabetical  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    }
                }
                path.m_segments.push_front(imax);
                len += imax->m_seq.size();
            } 
	        
            path.m_left = 0;
            path.m_right = path.m_segments.back()->m_seq.size()-1;
            path.m_current_pos.m_segmp =  path.m_segments.front();
            path.m_current_pos.m_pos = 0;
            path.m_current_seg = 0;
            path.m_len = len;
	
            return path;
        }
       
        Path ExpandToMax(GFAIterator start) {
            Path path;
            path.m_segments.push_back(start);
            int len = start->m_seq.size();
        
            //expand to right   
            while(!path.m_segments.back()->m_right_connections.empty()) {
                auto& last = path.m_segments.back();
                GFAIterator imax = last->m_right_connections.front();
                int maxlen = imax->m_right_len+imax->m_seq.size();
                size_t maxkmer = imax->m_right_kmer_count+imax->m_kmer_count;
                for(auto it = next(last->m_right_connections.begin()); it != last->m_right_connections.end(); ++it) {
                    GFAIterator i = *it;
                    int ilen = i->m_right_len+i->m_seq.size();
                    size_t ikmer = i->m_right_kmer_count+i->m_kmer_count;
                    if(ikmer < maxkmer) {
                        continue;
                    } else if(ikmer > maxkmer) {        //maxcount  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(ilen < maxlen) {
                        continue;
                    } else if(ilen > maxlen) {          //longest   
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(i->m_seq < imax->m_seq) { //alphabetical  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    }
            }
                path.m_segments.push_back(imax);
                len += imax->m_seq.size();
            }

            //expand to left    
            while(!path.m_segments.front()->m_left_connections.empty()) {
                auto& first = path.m_segments.front();
                GFAIterator imax = first->m_left_connections.front();
                int maxlen = imax->m_left_len+imax->m_seq.size();
                size_t maxkmer = imax->m_left_kmer_count+imax->m_kmer_count;
                for(auto it = next(first->m_left_connections.begin()); it != first->m_left_connections.end(); ++it) {
                    GFAIterator i = *it;
                    int ilen = i->m_left_len+i->m_seq.size();
                    size_t ikmer = i->m_left_kmer_count+i->m_kmer_count;
                    if(ikmer < maxkmer) {
                        continue;
                    } else if(ikmer > maxkmer) {        //maxcount  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(ilen < maxlen) {
                        continue;
                    } else if(ilen > maxlen) {          //longest   
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    } else if(i->m_seq < imax->m_seq) { //alphabetical  
                        imax = i;
                        maxlen = ilen;
                        maxkmer = ikmer;
                    }
                }
                path.m_segments.push_front(imax);
                len += imax->m_seq.size();
            } 

            path.m_left = 0;
            path.m_right = path.m_segments.back()->m_seq.size()-1;
            path.m_current_pos.m_segmp =  path.m_segments.front();
            path.m_current_pos.m_pos = 0;
            path.m_current_seg = 0;
            path.m_len = len;
            
            return path;
        }

        list<Path> ExpandToFork(GFAIterator start, bool toright, size_t maxp) {
            list<Path> expansion(1); 
            auto pathp = &expansion.back(); 
            pathp->m_segments.push_back(start);

            size_t total = 1;
            if(toright) {
                pathp->m_current_pos.m_segmp = start;
                pathp->m_current_pos.m_pos = 0;
                pathp->m_current_seg = 0;
                pathp->m_right = start->m_seq.size()-1;
                pathp->m_left = 0;

                bool keep_doing = true;
                while(keep_doing && total <= maxp) {
                    keep_doing = false;
                    for(auto& path : expansion) {
                        auto& segs = path.m_segments; 
                        if(segs.back()->LeftFork())
                            continue;
                        auto rcp = &segs.back()->m_right_connections;
                        for(auto it = rcp->begin(); it != rcp->end(); ++it) {
                            keep_doing = true;
                            auto pathp = &path;
                            if(next(it) != rcp->end()) {
                                expansion.push_front(path);
                                pathp = &expansion.front();
                                ++total;
                            }
                            pathp->m_segments.push_back(*it);
                        }
                    }
                }
            } else {
                bool keep_doing = true;
                while(keep_doing && total <= maxp) {
                    keep_doing = false;
                    for(auto& path : expansion) {
                        auto& segs = path.m_segments;                    
                        if(segs.front()->RightFork())
                            continue;
                        auto lcp = &path.m_segments.front()->m_left_connections;
                        for(auto it = lcp->begin(); it != lcp->end(); ++it) {
                            keep_doing = true;
                            auto pathp = &path;
                            if(next(it) != lcp->end()) {
                                expansion.push_front(path);
                                pathp = &expansion.front();
                                ++total;
                            }
                            pathp->m_segments.push_front(*it);
                        }
                    }
                }
            }

            for(auto it_loop = expansion.begin(); it_loop != expansion.end(); ) {
                auto it = it_loop++;
                auto& path = *it;
                if(toright) {
                    if(!path.m_segments.back()->LeftFork()) {
                        expansion.erase(it);
                        continue;
                    }
                    path.m_current_pos.m_pos = 0;
                    path.m_current_seg = 0;
                } else {
                    if(!path.m_segments.front()->RightFork()) {
                        expansion.erase(it);
                        continue;
                    }
                    path.m_current_pos.m_pos = path.m_segments.back()->m_seq.size()-1;
                    path.m_current_seg = path.m_segments.size()-1;
                }
                path.m_current_pos.m_segmp = start;
                path.m_left = 0;
                path.m_right = path.m_segments.back()->m_seq.size()-1;
                path.m_len = 0;
                for(auto iseg :path.m_segments )
                    path.m_len += iseg->m_seq.size();  
            }

            return expansion;            
        }
            
        list<Path> Expand(GFAIterator start, int start_pos, int to_left, int to_right, size_t maxp = numeric_limits<size_t>::max(), bool repeat_check = false) const {
            list<Path> expansion(1); 
            auto pathp = &expansion.back(); 
            pathp->m_current_pos.m_segmp = start;
            pathp->m_current_pos.m_pos = start_pos;
            pathp->m_current_seg = 0;

            pathp->m_segments.push_back(start);
            pathp->m_right = min(start_pos+to_right, (int)start->m_seq.size()-1);
            pathp->m_left = max(start_pos-to_left, 0);
            size_t total = 1;
                
            //expand to right   
            bool keep_doing = true;
            pathp->m_len = pathp->m_right-start_pos; // rlen 
            
            while(keep_doing && total <= maxp) {
                keep_doing = false;
                for(auto& path : expansion) {
                    if(path.m_len < to_right) {
                        auto& segs = path.m_segments;                    
                        if(repeat_check && segs.back()->m_cyclical && segs.back()->RightFork())                       //can go out of cycle
                            continue;
                        auto rcp = &segs.back()->m_right_connections;
                        for(auto it = rcp->begin(); it != rcp->end(); ++it) {
                            if(repeat_check && (*it)->m_cyclical && find(segs.begin(), segs.end(), *it) != segs.end()) // full cycle
                                continue;

                            keep_doing = true;
                            auto pathp = &path;
                            if(next(it) != rcp->end()) {
                                expansion.push_front(path);
                                pathp = &expansion.front();
                                ++total;
                            }
                            pathp->m_segments.push_back(*it);
                            pathp->m_right = min(to_right-pathp->m_len-1, (int)(*it)->m_seq.size()-1);
                            pathp->m_len += pathp->m_right+1;
                        }
                    }
                }
            }

            //expand to left    
            keep_doing = true;
            for(auto& path : expansion)
                path.m_len = start_pos-pathp->m_left; //llen

            while(keep_doing && total <= maxp) {
                keep_doing = false;
                for(auto& path : expansion) {
                    if(path.m_len < to_left) {
                        auto& segs = path.m_segments;                    
                        if(repeat_check && segs.front()->m_cyclical && segs.front()->LeftFork())                       //can go out of cycle
                            continue;
                        auto lcp = &path.m_segments.front()->m_left_connections;
                        for(auto it = lcp->begin(); it != lcp->end(); ++it) {
                            if(repeat_check && (*it)->m_cyclical && find(segs.begin(), segs.end(), *it) != segs.end()) // full cycle
                                continue;

                            keep_doing = true;
                            auto pathp = &path;
                            if(next(it) != lcp->end()) {
                                expansion.push_front(path);
                                pathp = &expansion.front();
                                ++total;
                            }
                            pathp->m_segments.push_front(*it);
                            ++(pathp->m_current_seg);
                            pathp->m_left = max((int)(*it)->m_seq.size()-(to_left-pathp->m_len), 0);
                            pathp->m_len += (*it)->m_seq.size()-pathp->m_left;
                        }
                    }
                }
            }

            for(auto& path : expansion) {
                if(path.m_segments.size() == 1)
                    path.m_len = path.m_right-path.m_left+1;
                else
                    path.m_len = path.m_segments.front()->m_seq.size()-path.m_left+path.m_right+1;
                for(int i = 1; i < (int)path.m_segments.size()-1; ++i)
                    path.m_len += path.m_segments[i]->m_seq.size();            
            }

            return expansion;
        } 
    
        SVarNum NumberOfVariants() {
            SVarNum total = 0;

            unordered_map<GFAIterator, SVarNum, SGFAIteratorHash> counts;
            for(auto it = begin(); it != end(); ++it) {
                if(it->m_left_connections.empty())
                    counts[it] = 1;
            }

            while(!counts.empty()) {
                unordered_map<GFAIterator, SVarNum, SGFAIteratorHash> new_counts;
                for(auto& count : counts) {
                    for(auto rc : count.first->m_right_connections) {
                        new_counts[rc] += count.second;
                    }
                    if(count.first->m_right_connections.empty()) {
                        total += count.second;
                    }
                }
                counts.swap(new_counts);
            }
        
            return total;
        }
    
        bool CheckConnections() const {
            for(auto it = begin(); it != end(); ++it) {
                for(auto rc : it->m_right_connections) {
                    if(count(rc->m_left_connections.begin(), rc->m_left_connections.end(), it) != 1)
                        return false;
                }
                for(auto rc : it->m_left_connections) {
                    if(count(rc->m_right_connections.begin(), rc->m_right_connections.end(), it) != 1)
                        return false;
                }
            }
            
            return true;
        }

        static void RemoveLinksToSegment(GFAIterator it) {
            for(auto lc : it->m_left_connections)
                lc->m_right_connections.remove(it);
            for(auto rc : it->m_right_connections)
                rc->m_left_connections.remove(it);                
        } 
        static void LinkSegments(GFAIterator left, GFAIterator right) {
            if(find(left->m_right_connections.begin(), left->m_right_connections.end(), right) == left->m_right_connections.end()) { // link if not already linked          
                left->m_right_connections.push_front(right);
                right->m_left_connections.push_front(left);
            }
        }
        static void UnLinkSegments(GFAIterator left, GFAIterator right) {
            left->m_right_connections.remove(right);
            right->m_left_connections.remove(left);
        }
        static void TransferRightLinks(GFAIterator from, GFAIterator to) {
            for(auto rc : from->m_right_connections) {
                rc->m_left_connections.remove(from);
                LinkSegments(to, rc);
            }
            from->m_right_connections.clear();
        }
        static void TransferLeftLinks(GFAIterator from, GFAIterator to) {
            for(auto lc : from->m_left_connections) {
                lc->m_right_connections.remove(from);
                LinkSegments(lc, to);
            }
            from->m_left_connections.clear();
        }
        void RemoveSegment(GFAIterator it) {
            --m_size;
            RemoveLinksToSegment(it);
            erase(it);
        }
        GFAIterator PushSegmentBack(const GFASegment& segm) {
            ++m_size;
            push_back(segm);
            back().m_num = ++m_max_num;
            return prev(end());
        }
        GFAIterator PushSegmentFront(const GFASegment& segm) {
            ++m_size;
            push_front(segm);
            front().m_num = ++m_max_num;
            return begin();
        }

        void CalculateChainLength() {
            stack<GFAIterator> current_lsegs;
            stack<GFAIterator> current_rsegs;
            for(auto iseg = begin(); iseg != end(); ++iseg) {
                if(iseg->m_left_connections.empty())
                    current_lsegs.push(iseg);
                if(iseg->m_right_connections.empty())
                    current_rsegs.push(iseg);
                iseg->m_left_len = 0;                
                iseg->m_right_len = 0;
            }

            while(!current_lsegs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_lsegs.empty()) {
                    GFAIterator iseg = current_lsegs.top();
                    current_lsegs.pop();
                    for(auto jseg : iseg->m_right_connections) {
                        int next_len = iseg->m_left_len+iseg->m_seq.size();
                        if(next_len > jseg->m_left_len) {
                            jseg->m_left_len = next_len;
                            next_step_segs.push(jseg);
                        }
                    }
                }
                std::swap(current_lsegs, next_step_segs);
            }

            while(!current_rsegs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_rsegs.empty()) {
                    GFAIterator iseg = current_rsegs.top();
                    current_rsegs.pop();
                    for(auto jseg : iseg->m_left_connections) {
                        int next_len = iseg->m_right_len+iseg->m_seq.size();
                        if(next_len > jseg->m_right_len) {
                            jseg->m_right_len = next_len;
                            next_step_segs.push(jseg);
                        }
                    }
                }
                std::swap(current_rsegs, next_step_segs);
            }
        }

        void ClipToCodons() {
            bool keep_doing = true;
            while(keep_doing) {
                keep_doing = false;
                for(auto it_loop = begin(); it_loop != end(); ) {
                    auto it = it_loop++;
                    if(it->m_frame < 0)
                        throw runtime_error("Unknown frame");                                        
                    if(it->m_left_connections.empty() && it->m_frame > 0) {
                        int lclip = 3-it->m_frame;
                        if((int)it->m_seq.size() > lclip) {
                            it->m_seq.ClipLeft(lclip);
                            it->m_frame = 0;
                        } else {
                            RemoveSegment(it);
                            keep_doing = true;
                            continue;
                        }
                    }
                    if(it->m_right_connections.empty() && (it->m_frame+it->m_seq.size())%3 > 0) {
                        int rclip = (it->m_frame+it->m_seq.size())%3;
                        if((int)it->m_seq.size() > rclip) {
                            it->m_seq.ClipRight(rclip);
                        } else {
                            RemoveSegment(it);
                            keep_doing = true;
                            continue;
                        }
                    }
                }
            }
        }

        void TranslateToAA(const GeneticCode& gcode, DBGraph& dbg_graph) {

            //move forks to codon boundaries
            for(auto iloop = begin(); iloop != end(); ) {
                auto iseg = iloop++;
                size_t seg_len = iseg->m_seq.size();
                int frame = iseg->m_frame;
                iseg->m_frame = -1;
                if(iseg->LeftFork() && iseg->m_frame > 0) {
                    size_t shift = 3-iseg->m_frame;
                    if(seg_len > shift || iseg->m_right_connections.empty()) { // long or right end
                        auto delta = iseg->m_seq.substr(0, min(shift, seg_len));
                        for(auto lcloop = iseg->m_left_connections.begin(); lcloop != iseg->m_left_connections.end(); ) {
                            auto lc = *lcloop++;
                            auto newseg = PushSegmentFront(GFASegment(delta));
                            newseg->m_group = iseg->m_group;
                            UnLinkSegments(lc, iseg);
                            LinkSegments(lc, newseg);
                            LinkSegments(newseg, iseg);
                        }
                    } else {  // short and in the middle
                        for(auto lc : iseg->m_left_connections) {
                            for(auto rc : iseg->m_right_connections) {
                                auto newseg = PushSegmentFront(*iseg);
                                LinkSegments(lc, newseg);
                                LinkSegments(newseg, rc);
                            }
                        }
                    }
                    if(seg_len > shift) {
                        iseg->m_seq.ClipLeft(shift);
                        frame = 0;           // needed for next step
                        seg_len -= shift;
                    } else {
                        RemoveSegment(iseg);
                        continue;
                    }
                }

                if(iseg->RightFork() && (frame+seg_len)%3 > 0) {
                    size_t shift = (frame+seg_len)%3;
                    if(seg_len > shift || iseg->m_left_connections.empty()) { // long or left end
                        auto delta = iseg->m_seq.substr(seg_len-shift);
                        for(auto rcloop = iseg->m_right_connections.begin(); rcloop != iseg->m_right_connections.end(); ) {
                            auto rc = *rcloop++;
                            auto newseg = PushSegmentFront(GFASegment(delta));
                            newseg->m_group = iseg->m_group;
                            UnLinkSegments(iseg, rc);
                            LinkSegments(newseg, rc);
                            LinkSegments(iseg, newseg);
                        } 
                    } else {  // short and in the middle
                        for(auto lc : iseg->m_left_connections) {
                            for(auto rc : iseg->m_right_connections) {
                                auto newseg = PushSegmentFront(*iseg);
                                LinkSegments(lc, newseg);
                                LinkSegments(newseg, rc);
                            }
                        }
                    }
                    if(seg_len > shift) {
                        iseg->m_seq.ClipRight(shift);
                        seg_len -= shift;
                    } else {
                        RemoveSegment(iseg);
                        continue;
                    }
                }
            }
            MergeSimpleLinks();  //after this all segments contain whole codons

            GenerateKmers(dbg_graph);
            for(auto& segm : *this) {
                string prot_seq;
                for(auto& base : segm.m_seq)
                    prot_seq.push_back(base.m_nt);
                prot_seq = gcode.Translate(prot_seq, false);

                SegSeq segm_seq;
                for(size_t p = 0; p < prot_seq.size(); ++p) {
                    segm_seq.emplace_back();
                    SegBase& base = segm_seq.back();
                    base.m_nt = prot_seq[p];
                    for(int frame = 0; frame < 3; ++frame) {
                        size_t nucp = 3*p+frame;
                        base.m_left_kmers.splice_after(base.m_left_kmers.before_begin(), segm.m_seq[nucp].m_left_kmers);
                        base.m_right_kmers.splice_after(base.m_right_kmers.before_begin(), segm.m_seq[nucp].m_right_kmers);
                    }
                }
                std::swap(segm_seq, segm.m_seq);
            }
            m_is_aa = true;
            MergeRedundantLinks();
            
            //calculate coverage;   
            for(auto& seg : *this) {
                seg.m_kmer_count = 0;
                for(auto& base : seg.m_seq) {
                    size_t lcount = 0;
                    size_t rcount = 0;
                    for(Node& node : base.m_left_kmers)
                        lcount += dbg_graph.Abundance(node);                    
                    for(Node& node : base.m_right_kmers)
                        rcount += dbg_graph.Abundance(node);                    
                    seg.m_kmer_count += max(lcount, rcount);
                }
            }            
            CalculateCoverageAndEnumerateSegments();
        }

        void SnpsToAmbig() {
            for(auto iseg = this->begin(); iseg != this->end(); ++iseg) {
                auto& seg = *iseg;
                if(seg.LeftConnectionsNum() == 2) {
                    for(auto it_loop = seg.m_left_connections.begin(); it_loop != seg.m_left_connections.end(); ) {
                        auto it = it_loop++;
                        auto is = *it;
                        if(!is->m_left_connections.empty() || !is->RightSingle())
                            break;
                        string ambig(1,is->m_seq.back().m_nt);
                        for(auto jt = seg.m_left_connections.begin(); jt != seg.m_left_connections.end(); ++jt) {
                            auto js = *jt;
                            if(jt == it || js->m_seq.size() < is->m_seq.size())
                                continue;
                            if(js->m_seq.size() == 1 || is->m_seq.size() == 1 || equal(is->m_seq.rbegin()+1, is->m_seq.rend(), js->m_seq.rbegin()+1, [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; })) {
                                ambig.push_back(js->m_seq.back().m_nt);
                                std::sort(ambig.begin(), ambig.end());
                                ambig.erase(std::unique(ambig.begin(), ambig.end()), ambig.end());
                                if(ambig.size() > 1)
                                    js->m_seq.back().m_nt = ToAmbiguousIUPAC[ambig];
                                js->m_seq.back().m_left_kmers.splice_after(js->m_seq.back().m_left_kmers.before_begin(), is->m_seq.back().m_left_kmers);
                                js->m_seq.back().m_right_kmers.splice_after(js->m_seq.back().m_right_kmers.before_begin(), is->m_seq.back().m_right_kmers);
                                js->m_kmer_count = max(js->m_kmer_count, is->m_kmer_count);
                                RemoveSegment(is);
                                break;
                            }
                        }
                    } 
                }

                if(seg.RightConnectionsNum() == 2) {
                    for(auto it_loop = seg.m_right_connections.begin(); it_loop != seg.m_right_connections.end(); ) {
                        auto it = it_loop++;
                        auto is = *it;
                        if(!is->m_right_connections.empty() || !is->LeftSingle())
                            break;
                        string ambig(1,is->m_seq.front().m_nt);
                        for(auto jt = seg.m_right_connections.begin(); jt != seg.m_right_connections.end(); ++jt) {
                            auto js = *jt;
                            if(jt == it || js->m_seq.size() < is->m_seq.size())
                                continue;                            
                            if(js->m_seq.size() == 1 || is->m_seq.size() == 1 || equal(is->m_seq.begin()+1, is->m_seq.end(), js->m_seq.begin()+1, [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; })) {
                                ambig.push_back(js->m_seq.front().m_nt);
                                std::sort(ambig.begin(), ambig.end());
                                ambig.erase(std::unique(ambig.begin(), ambig.end()), ambig.end());
                                if(ambig.size() > 1)
                                    js->m_seq.front().m_nt = ToAmbiguousIUPAC[ambig];
                                js->m_seq.front().m_left_kmers.splice_after(js->m_seq.front().m_left_kmers.before_begin(), is->m_seq.front().m_left_kmers);
                                js->m_seq.front().m_right_kmers.splice_after(js->m_seq.front().m_right_kmers.before_begin(), is->m_seq.front().m_right_kmers);
                                js->m_kmer_count = max(js->m_kmer_count, is->m_kmer_count);
                                RemoveSegment(is);
                                break;
                            }
                        }
                    }
                }                
            }
            MergeForks();

            for(GFAIterator iseg = begin(); iseg != end(); ) {
                bool rsnp = false;
                auto inext = (iseg->m_right_connections.empty() || iseg->m_right_connections.front()->m_right_connections.empty()) ? 
                               end() : iseg->m_right_connections.front()->m_right_connections.front();
                for(auto rc : iseg->m_right_connections) {
                    rsnp = rc->LeftSingle() && rc->RightSingle() && rc->m_seq.size() == 1 && rc->m_right_connections.front() == inext;
                    if(!rsnp)
                        break;
                }
                if(rsnp && iseg->RightConnectionsNum() == inext->LeftConnectionsNum()) {
                    auto inext = iseg->m_right_connections.front()->m_right_connections.front();
                    string ambig;
                    auto& seq = iseg->m_seq;
                    seq.emplace_back();
                    size_t kcount = 0;
                    for(auto irc = iseg->m_right_connections.begin(); irc != iseg->m_right_connections.end(); ) {
                        auto rc = *irc++;
                        auto& base = rc->m_seq[0];
                        ambig.push_back(base.m_nt);
                        seq.back().m_left_kmers.splice_after(seq.back().m_left_kmers.before_begin(), base.m_left_kmers);
                        seq.back().m_right_kmers.splice_after(seq.back().m_right_kmers.before_begin(), base.m_right_kmers);
                        kcount = max(kcount, rc->m_kmer_count);
                        RemoveSegment(rc);
                    }
                    std::sort(ambig.begin(), ambig.end());
                    seq.back().m_nt = ToAmbiguousIUPAC[ambig];
                    seq += inext->m_seq;
                    iseg->m_kmer_count += kcount+inext->m_kmer_count;
                    iseg->m_right_kmer_count = inext->m_right_kmer_count;
                    iseg->m_right_len = inext->m_right_len;
                    TransferRightLinks(inext, iseg);
                    RemoveSegment(inext);
                } else {
                    ++iseg;
                }
            }            
        }    

        int RemoveShortChains(int minlen, bool mark_only = false) {
            int maxlen = 0;
            for(auto it_loop = this->begin(); it_loop != this->end(); ) {
                auto it = it_loop++;
                auto& seg = *it;
                int len = seg.m_left_len+seg.m_seq.size()+seg.m_right_len;
                maxlen = max(maxlen, len);
                if(len < minlen) {
                    if(mark_only) {
                        if(!it->m_marked_for_erase) {
                            it->m_marked_for_erase = true;
                            RemoveLinksToSegment(it);
                            it->m_left_connections.clear();
                            it->m_right_connections.clear();
                        }
                    } else {
                        RemoveSegment(it);
                    }                    
                }
            }
            return maxlen;
        }

        void MergeRedundantDuplicates() {
            bool keep_doing = true;
            while(keep_doing) {
                keep_doing = false;
                for(auto iseg = begin(); iseg != end(); ++iseg) {
                    if(iseg->RightFork()) {
                        for(auto it = iseg->m_right_connections.begin(); it != iseg->m_right_connections.end(); ++it) {
                            auto is = *it;
                            if(is == iseg)     // self      
                                continue;
                            auto ilc = is->m_left_connections;
                            ilc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                            for(auto jt = next(it); jt != iseg->m_right_connections.end(); ) {
                                auto js = *jt++;
                                if(js == iseg) // self      
                                    continue;
                                if(is->m_right_check != js->m_right_check) // mixed status
                                    continue;
                                if((is->m_copy_of != nullptr && is->m_copy_of == js->m_copy_of) || is->m_copy_of == &(*js) || js->m_copy_of == &(*is)) {  // copies of the same
                                    auto jlc = js->m_left_connections;
                                    jlc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                                    // only     if i and j connect back to the same segments    
                                    if(ilc != jlc)
                                        continue;

                                    keep_doing = true; 
                                    js->m_marked_for_erase = true;
                                    TransferRightLinks(js, is);
                                    RemoveLinksToSegment(js);
                                    js->m_left_connections.clear();
                                }
                            }
                        }
                    }
                    if(iseg->LeftFork()) {
                        for(auto it = iseg->m_left_connections.begin(); it != iseg->m_left_connections.end(); ++it) {
                            auto is = *it;
                            if(is == iseg)     // self      
                                continue;
                            auto irc = is->m_right_connections;
                            irc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                            for(auto jt = next(it); jt != iseg->m_left_connections.end(); ) {
                                auto js = *jt++;
                                if(js == iseg) // self      
                                    continue;
                                if(is->m_left_check != js->m_left_check) // mixed status
                                    continue;
                                if((is->m_copy_of != nullptr && is->m_copy_of == js->m_copy_of) || is->m_copy_of == &(*js) || js->m_copy_of == &(*is)) {  // copies of the same
                                    auto jrc = js->m_right_connections;
                                    jrc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                                    // only     if i and j connect back to the same segments    
                                    if(irc != jrc)
                                        continue;

                                    keep_doing = true;
                                    js->m_marked_for_erase = true;
                                    TransferLeftLinks(js, is);
                                    RemoveLinksToSegment(js);
                                    js->m_right_connections.clear();
                                }
                            }
                        }
                    }
                }
            }
        }

        void RemoveRedundantPaths(size_t maxp) {
            for (auto it = begin(); it != end(); ++it) {
                if(it->m_left_connections.empty()) {
                    list<Path> direct_paths = ExpandToFork(it, true, maxp);
                    for(auto& dpath : direct_paths) {
                        int len = dpath.Length();
                        auto iback = dpath.m_segments.back();
                        list<Path> reverse_paths = Expand(iback, iback->m_seq.size()-1, len-1, 0, maxp);
                        for(auto& rpath : reverse_paths) {
                            int dnum = dpath.m_segments.size();
                            int rnum = rpath.m_segments.size();
                            if(len != rpath.Length() || rpath.m_segments[rnum-2] == dpath.m_segments[dnum-2])
                                continue;
                            if(dpath.Sequence() == rpath.Sequence()) {
                                for(int i = dnum-2; i >= 0; --i) {
                                    UnLinkSegments(dpath.m_segments[i], dpath.m_segments[i+1]);                                    
                                    if(!dpath.m_segments[i]->m_right_connections.empty())
                                        break;
                                }
                            }
                        }
                    }
                }
                if(it->m_right_connections.empty()) {
                    list<Path> direct_paths = ExpandToFork(it, false, maxp);
                    for(auto& dpath : direct_paths) {
                        int len = dpath.Length();
                        list<Path> reverse_paths = Expand(dpath.m_segments.front(), 0, 0, len-1, maxp);
                        for(auto& rpath : reverse_paths) {
                            if(len != rpath.Length() || rpath.m_segments[1] == dpath.m_segments[1])
                                continue;
                            if(dpath.Sequence() == rpath.Sequence()) {
                                for(int i = 1; i < (int)dpath.m_segments.size(); ++i) {
                                    UnLinkSegments(dpath.m_segments[i-1], dpath.m_segments[i]);                                    
                                    if(!dpath.m_segments[i]->m_left_connections.empty())
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }

        GFAIterator ReduceGraph(int min_len, TCopyInfo& copies, GFAIterator inext) {
            MergeRedundantDuplicates();
            if(min_len > 0) {
                CalculateChainLength();
                RemoveShortChains(min_len, true);   
            }

            GFAIteratorUSet erased;
            for(auto iseg = begin(); iseg != end(); ++iseg) {
                if(iseg->m_marked_for_erase)
                    erased.insert(iseg);
            }
            RecalculateCopyInfo(copies, erased);
            for(auto iseg : erased) {
                if(iseg == inext)
                    ++inext;
                RemoveSegment(iseg);
            }         

            return inext;
        }

        int AssignGroupNumber() {
            for(auto& seg : *this)
                seg.m_group = 0;
               
            int group = 0;
            for(auto& seg : *this) {
                if(seg.m_group > 0)
                    continue;
                seg.m_group = ++group;
                list<GFAIterator> links;
                links.insert(links.end(), seg.m_left_connections.begin(), seg.m_left_connections.end());
                links.insert(links.end(), seg.m_right_connections.begin(), seg.m_right_connections.end());
                while(!links.empty()) {
                    if(links.front()->m_group == 0) {
                        links.front()->m_group = group;
                        links.insert(links.end(), links.front()->m_left_connections.begin(), links.front()->m_left_connections.end());
                        links.insert(links.end(), links.front()->m_right_connections.begin(), links.front()->m_right_connections.end());
                    }
                    links.pop_front();
                }
            }
            
            return group;
        }

        void TrimGroups(double fraction, int minlen, unordered_set<const GFASegment*>* erasedp = nullptr, GFAIterator* ip = nullptr) {
            map<int, int> group_len;
            for(auto& seg : *this)
                group_len[seg.m_group] = max(group_len[seg.m_group], seg.m_left_len+(int)seg.m_seq.size()+seg.m_right_len);
            
            for(auto it_loop = this->begin(); it_loop != this->end(); ) {
                auto it = it_loop++;
                auto& seg = *it;
                int glen = group_len[seg.m_group];
                int seg_len = seg.m_left_len+seg.m_seq.size()+seg.m_right_len;
                if(seg_len < minlen || seg_len < fraction*glen) {
                    if(erasedp != nullptr)
                        erasedp->insert(&(*it));
                    if(ip != nullptr && *ip == it)
                        ++(*ip);
                    RemoveSegment(it);
                }
            }
        }

        TGFACollection SplitGroups() { // destroys source   
            TGFACollection splitted;
            map<int, TGFACollection::iterator> group_to_graph;

            while(!this->empty()) {
                int group = this->front().m_group;
                if(!group_to_graph.count(group)) {
                    splitted.emplace_front(Target(), KmerLen());
                    group_to_graph[group] = splitted.begin();
                }
                auto& dest = *group_to_graph[group];
                dest.splice(dest.end(), *this, this->begin());
            }

            for(auto& graph : splitted)
                graph.Size() = graph.size();

            return splitted;
        }

        void MergeSimpleLinks() {
            for(auto it_loop = this->begin(); it_loop != this->end(); ) {
                auto it = it_loop++;
                auto& seg = *it;
                if(seg.LeftSingle()) {  // exactly one link  
                    auto left_seg = seg.m_left_connections.front();
                    if(left_seg != it && left_seg->RightSingle()) { // exactly one link and not self 
                        left_seg->m_seq += seg.m_seq;
                        TransferRightLinks(it, left_seg);
                        RemoveSegment(it);
                        continue;
                    }
                }
                if(seg.RightSingle()) { // exactly one link        
                    auto right_seg = seg.m_right_connections.front();
                    if(right_seg != it && right_seg->LeftSingle()) {  // exactly one link and not self       
                        right_seg->m_seq = seg.m_seq+right_seg->m_seq;
                        right_seg->m_frame = seg.m_frame;
                        TransferLeftLinks(it, right_seg);
                        RemoveSegment(it);
                        continue;
                    }
                }                
            }
        }

        void MergeForks() {
            // merge simple links   
            MergeSimpleLinks();

            //merge redundand forks     
            bool keep_doing = true;
            while(keep_doing) {
                keep_doing = false;
                for(auto iseg = this->begin(); iseg != this->end(); ++iseg) {
                    auto& seg = *iseg;

                    if(seg.LeftFork()) {  // fork  
                        for(auto it_loop = seg.m_left_connections.begin(); it_loop != seg.m_left_connections.end(); ) {
                            auto it = it_loop++;
                            auto is = *it;
                            if(!is->m_left_connections.empty() || is->RightFork()) 
                                continue;
                            for(auto jt = seg.m_left_connections.begin(); jt != seg.m_left_connections.end(); ++jt) {
                                auto js = *jt;
                                if(jt == it || js->m_seq.size() < is->m_seq.size())
                                    continue;
                                if(equal(is->m_seq.rbegin(), is->m_seq.rend(), js->m_seq.rbegin(), [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; })) {
                                    RemoveSegment(is);
                                    break;
                                }
                            }
                        }
                    }

                    while(seg.LeftFork()) { // fork        
                        for(auto it = seg.m_left_connections.begin(); it != seg.m_left_connections.end(); ++it) {
                            auto is = *it;
                            if(is == iseg)     // self      
                                continue;
                            auto irc = is->m_right_connections;
                            irc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                            for(auto jt = next(it); jt != seg.m_left_connections.end(); ) {
                                auto js = *jt++;
                                if(js == iseg) // self  
                                    continue;

                                auto jrc = js->m_right_connections;
                                jrc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                                // only if i and j connect back to the same segments    
                                if(irc != jrc)
                                    continue;
                                
                                int ilen = is->m_seq.size();
                                int jlen = js->m_seq.size();
                                int len = min(ilen, jlen);
                                auto rslt = mismatch(is->m_seq.rbegin(), is->m_seq.rbegin()+len, js->m_seq.rbegin(), [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; });
                                int matches = rslt.first-is->m_seq.rbegin();
                                if(matches > 0) {
                                    if(ilen > matches) { // cut i and link parts 
                                        auto newi = PushSegmentBack(GFASegment(is->m_seq.substr(0, ilen-matches)));
                                        newi->m_group = seg.m_group;
                                        newi->m_cyclical = is->m_cyclical;
                                        TransferLeftLinks(is, newi);
                                        is->m_seq.ClipLeft(ilen-matches);
                                        if(is->m_frame >= 0) {
                                            newi->m_frame = is->m_frame;
                                            is->m_frame = (is->m_frame+ilen-matches)%3;
                                        }
                                        LinkSegments(newi, is);
                                    }
                                    if(jlen > matches) { // cut j and keep first part and link it to i          
                                        auto newj = PushSegmentBack(GFASegment(js->m_seq.substr(0, jlen-matches)));
                                        newj->m_group = seg.m_group;
                                        newj->m_cyclical = js->m_cyclical;
                                        if(js->m_cyclical)
                                            is->m_cyclical = true;
                                        newj->m_frame = js->m_frame;
                                        TransferLeftLinks(js, newj);
                                        LinkSegments(newj, is);
                                        js->m_seq.ClipLeft(jlen-matches);
                                    } else {             // move j links to i           
                                        if(js->m_cyclical)
                                            is->m_cyclical = true;
                                        TransferLeftLinks(js, is);
                                    }
                                    // transfer fork info and kmers from j to i  
                                    is->m_seq.SplitForksAndKmersFrom(js->m_seq);
                                    // delete j
                                    RemoveSegment(js);
                                    keep_doing = true;
                                }
                            }
                        }
                        if(seg.LeftSingle()) {      
                            auto left_seg = seg.m_left_connections.front();
                            if(left_seg->RightSingle()) { // only one connection left - merge    
                                seg.m_seq = left_seg->m_seq+seg.m_seq;
                                seg.m_frame = left_seg->m_frame;
                                TransferLeftLinks(left_seg, iseg);
                                RemoveSegment(left_seg);
                                continue;
                            }                                   
                        }
                        break;                            
                    }

                    if(seg.RightFork()) {  // fork   
                        for(auto it_loop = seg.m_right_connections.begin(); it_loop != seg.m_right_connections.end(); ) {
                            auto it = it_loop++;
                            auto is = *it;
                            if(!is->m_right_connections.empty() || is->LeftFork())   
                                continue;
                            for(auto jt = seg.m_right_connections.begin(); jt != seg.m_right_connections.end(); ++jt) {
                                auto js = *jt;
                                if(jt == it || js->m_seq.size() < is->m_seq.size())
                                    continue;                            
                                if(equal(is->m_seq.begin(), is->m_seq.end(), js->m_seq.begin(), [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; })) {
                                    RemoveSegment(is);
                                    break;
                                }
                            }
                        }
                    }

                    while(seg.RightFork()) { // fork     
                        for(auto it = seg.m_right_connections.begin(); it != seg.m_right_connections.end(); ++it) {
                            auto is = *it;
                            if(is == iseg)     // self      
                                continue;
                            auto ilc = is->m_left_connections;
                            ilc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                            for(auto jt = next(it); jt != seg.m_right_connections.end(); ) {
                                auto js = *jt++;
                                if(js == iseg) // self      
                                    continue;

                                auto jlc = js->m_left_connections;
                                jlc.sort([](const GFAIterator& a, const GFAIterator& b){ return &(*a) < &(*b); }); // compare physical pointers 
                                // only if i and j connect back to the same segments    
                                if(ilc != jlc)
                                    continue;

                                int ilen = is->m_seq.size();
                                int jlen = js->m_seq.size();
                                int len = min(ilen, jlen);
                                auto rslt = mismatch(is->m_seq.begin(), is->m_seq.begin()+len, js->m_seq.begin(), [](const SegBase& a, const SegBase& b) { return a.m_nt == b.m_nt; });
                                int matches = rslt.first-is->m_seq.begin();
                                if(matches > 0) {
                                    if(ilen > matches) { // cut i and link parts            
                                        auto newi = PushSegmentBack(GFASegment(is->m_seq.substr(matches)));
                                        newi->m_group = seg.m_group;
                                        newi->m_cyclical = is->m_cyclical;
                                        if(is->m_frame >= 0)
                                            newi->m_frame = (is->m_frame+matches)%3;
                                        TransferRightLinks(is, newi);
                                        is->m_seq.ClipRight(ilen-matches);
                                        LinkSegments(is, newi);
                                    }
                                    if(jlen > matches) { // cut j and keep second part and link it to i         
                                        auto newj = PushSegmentBack(GFASegment(js->m_seq.substr(matches)));
                                        newj->m_group = seg.m_group;
                                        newj->m_cyclical = js->m_cyclical;
                                        if(js->m_frame >= 0)
                                            newj->m_frame = (js->m_frame+matches)%3;
                                        if(js->m_cyclical)
                                            is->m_cyclical = true;
                                        TransferRightLinks(js, newj);
                                        LinkSegments(is, newj);
                                        js->m_seq.ClipRight(jlen-matches);
                                    } else {             // move j links to i           
                                        if(js->m_cyclical)
                                            is->m_cyclical = true;
                                        TransferRightLinks(js, is);
                                    }
                                    // transfer fork info and kmers from j to i 
                                    is->m_seq.SplitForksAndKmersFrom(js->m_seq);
                                    // delete j 
                                    RemoveSegment(js);
                                    keep_doing = true;
                                }
                            }
                        }
                        if(seg.RightSingle()) {        
                            auto right_seg = seg.m_right_connections.front();
                            if(right_seg->LeftSingle()) { // only one connection left - merge    
                                seg.m_seq += right_seg->m_seq;
                                TransferRightLinks(right_seg, iseg);
                                RemoveSegment(right_seg);
                                continue;
                            }
                        }
                        break;                            
                    }   
                }
            } 

            MergeSimpleLinks();
        }

        void MergeRedundantLinks() {
            MergeForks();
            CalculateChainLength();            
        }

        GFAIteratorUSet LeftBranchSegments(GFAIterator it) {
            GFAIteratorUSet branch;
            GFAIteratorUSet inlets;

            stack<GFAIterator> segments;
            segments.push(it);
            branch.insert(it);

            while(!segments.empty()) {
                auto iseg = segments.top();
                segments.pop();
                for(auto cn : iseg->m_left_connections) {
                    if(branch.insert(cn).second) {
                        segments.push(cn);
                        for(auto inl : cn->m_right_connections)
                            inlets.insert(inl);                        
                    }
                }
            }

            for(auto inl : inlets) {
                if(!branch.count(inl)) {
                    branch.clear();
                    return branch;
                }
            }

            return branch;
        }

        GFAIteratorUSet RightBranchSegments(GFAIterator it) {
            GFAIteratorUSet branch;
            GFAIteratorUSet inlets;

            stack<GFAIterator> segments;
            segments.push(it);
            branch.insert(it);

            while(!segments.empty()) {
                auto iseg = segments.top();
                segments.pop();
                for(auto cn : iseg->m_right_connections) {
                    if(branch.insert(cn).second) {
                        segments.push(cn);
                        for(auto inl : cn->m_left_connections)
                            inlets.insert(inl);                        
                    }
                }
            }

            for(auto inl : inlets) {
                if(!branch.count(inl)) {
                    branch.clear();
                    return branch;
                }
            }

            return branch;
        }
        
        bool RemoveHair(DBGraph& dbg_graph, double eps) {
            GenerateKmersAndScores(dbg_graph);

            bool deleted = false;

            for(auto iseg = begin(); iseg != end(); ++iseg) {
                for(auto iloop = iseg->m_left_connections.begin(); iloop != iseg->m_left_connections.end(); ) {
                    auto jseg = *iloop++;
                    int jcount = jseg->m_kmer_count+jseg->m_left_kmer_count;
                    if(eps*iseg->m_left_kmer_count > jcount) {
                        auto branch = LeftBranchSegments(jseg);
                        if(!branch.empty()) {
                            deleted = true;
                            if(jseg->RightFork()) {
                                UnLinkSegments(jseg, iseg);
                            } else {
                                for(auto i : branch)
                                    RemoveSegment(i);
                            }
                        }
                    }
                }
                for(auto iloop = iseg->m_right_connections.begin(); iloop != iseg->m_right_connections.end(); ) {
                    auto jseg = *iloop++;
                    int jcount = jseg->m_kmer_count+jseg->m_right_kmer_count;
                    if(eps*iseg->m_right_kmer_count > jcount) {
                        auto branch = RightBranchSegments(jseg);
                        if(!branch.empty()) {
                            deleted = true;
                            if(jseg->LeftFork()) {
                                UnLinkSegments(iseg, jseg);
                            } else {
                                for(auto i : branch)
                                    RemoveSegment(i);
                            }
                        }
                    }
                }
            }

            return deleted;
        }

        void EnumerateSegments() {
            m_max_num = 0;
            for(auto& seg : *this)
                seg.m_num = ++m_max_num;
            m_size = m_max_num;
        }

        void CalculateCoverageAndEnumerateSegments() {
            stack<GFAIterator> current_lsegs;
            stack<GFAIterator> current_rsegs;
            for(auto iseg = begin(); iseg != end(); ++iseg) {
                if(iseg->m_left_connections.empty())
                    current_lsegs.push(iseg);
                if(iseg->m_right_connections.empty())
                    current_rsegs.push(iseg);
                iseg->m_left_kmer_count = 0;                
                iseg->m_right_kmer_count = 0;
            }

            while(!current_lsegs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_lsegs.empty()) {
                    GFAIterator iseg = current_lsegs.top();
                    current_lsegs.pop();
                    for(auto jseg : iseg->m_right_connections) {
                        size_t next_count = iseg->m_left_kmer_count+iseg->m_kmer_count;
                        if(next_count > jseg->m_left_kmer_count) {
                            jseg->m_left_kmer_count = next_count;
                            next_step_segs.push(jseg);
                        }
                    }
                }
                std::swap(current_lsegs, next_step_segs);
            }

            while(!current_rsegs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_rsegs.empty()) {
                    GFAIterator iseg = current_rsegs.top();
                    current_rsegs.pop();
                    for(auto jseg : iseg->m_left_connections) {
                        size_t next_count = iseg->m_right_kmer_count+iseg->m_kmer_count;
                        if(next_count > jseg->m_right_kmer_count) {
                            jseg->m_right_kmer_count = next_count;
                            next_step_segs.push(jseg);
                        }
                    }
                }
                std::swap(current_rsegs, next_step_segs);
            }

            EnumerateSegments();
        }

        SegSeq ExtendToFirstFork(Node node, unordered_set<Node, typename Node::Hash>& used_node_indexes, GraphDigger& graphdigger) const {
            SegSeq s;
            while(true) {
                vector<Successor> successors = graphdigger.Graph().GetNodeSuccessors(node);
                graphdigger.FilterNeighbors(successors, true);
                if(successors.empty())
                    return s;                            
                if(successors.size() != 1)
                    return s;
                Node rev_node = DBGraph::ReverseComplement(successors[0].m_node);
                vector<Successor> predecessors = graphdigger.Graph().GetNodeSuccessors(rev_node);
                graphdigger.FilterNeighbors(predecessors, true);
                
                if(predecessors.size() == 1 && 
                   predecessors[0].m_node == DBGraph::ReverseComplement(node) &&
                   used_node_indexes.insert(successors[0].m_node.DropStrand()).second) { // good extension      
                    
                    node = successors[0].m_node;
                    SegBase base;
                    base.m_nt = successors[0].m_nt;
                    base.m_right_kmers.push_front(node);
                    s.push_back(base);
                } else {
                    return s;
                }
            }
        }

        void ExtendToFirstFork(GraphDigger& graphdigger) {
            map<int, unordered_set<Node, typename Node::Hash>> used_node_indexes;
            
            for(auto& seg : *this) {
                auto& used_ind = used_node_indexes[seg.m_group] ;
                for(auto& base : seg.m_seq) {
                    for(Node& node : base.m_left_kmers)                        
                        used_ind.insert(node.DropStrand());
                    for(Node& node : base.m_right_kmers)                        
                        used_ind.insert(node.DropStrand());
                }
            }        
        
            for(auto it = begin(); it != end(); ++it) {
                auto& seg = *it;
                int seq_len = seg.m_seq.size();
                if(seg.m_right_connections.empty()) {
                    list<Path> paths = Expand(it, seq_len-1, m_kmer_len-1, 0);
                    if(paths.size() == 1 && paths.front().Length() == m_kmer_len) {
                        TKmer kmer(paths.front().Sequence());
                        auto node = graphdigger.Graph().GetNode(kmer);
                        SegSeq ext = ExtendToFirstFork(node, used_node_indexes[seg.m_group], graphdigger);
                        seg.m_seq += ext;
                    }
                }
                if(seg.m_left_connections.empty()) {
                    list<Path> paths = Expand(it, 0, 0, m_kmer_len-1);
                    if(paths.size() == 1 && paths.front().Length() == m_kmer_len) {
                        TKmer kmer(paths.front().Sequence());
                        auto node = graphdigger.Graph().GetNode(kmer);
                        node = DBGraph::ReverseComplement(node);
                        SegSeq ext = ExtendToFirstFork(node, used_node_indexes[seg.m_group], graphdigger);
                        ext.ReverseComplement();
                        ext += seg.m_seq;
                        ext.swap(seg.m_seq);
                    }
                }
            }
        }
    
        string SegId(const GFASegment& seg) { 
            return m_acc+":"+to_string(seg.m_group)+":"+to_string(seg.m_num); 
        };

        void PrintGFA(ostream& out) {
            for(auto& seg : *this) {
                out << "S\t" << SegId(seg) << "\t";
                for(auto& base : seg.m_seq)
                    out << base.m_nt;
                double coverage = seg.m_kmer_count;
                if(m_is_aa)
                    coverage /= 3;
                out << "\tKC:i:" << coverage << "\n";
            }
            for(auto& seg : *this) {
                /*
                for(auto& lc : seg.m_left_connections)
                    out << "L\t" << SegId(seg) << "\t-\t" << SegId(*lc) << "\t-\t0M" << "\n"; 
                */                                   
                for(auto& rc : seg.m_right_connections)
                    out << "L\t" << SegId(seg) << "\t+\t" << SegId(*rc) << "\t+\t0M" << "\n";                                    
            }
        }
        void PrintAllVariants(ostream& out, unsigned max_variants) {
            typedef tuple<size_t, int> TPathScore;
            unordered_map<GFAIterator, set<TPathScore>, SGFAIteratorHash> seg_scores;
            size_t max_stack = 2*max_variants; // give some extra in case of identical sequences

            //propagate score from right to left
            stack<GFAIterator> current_segs;
            for(GFAIterator it = begin(); it != end(); ++it) {
                if(it->m_right_connections.empty()) { // start from right extrems
                    current_segs.push(it);
                    TPathScore iscore(it->m_kmer_count, it->m_seq.size());
                    seg_scores[it].insert(iscore);
                }
            }
            while(!current_segs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_segs.empty()) {
                    GFAIterator rseg = current_segs.top();
                    current_segs.pop();
                    auto& rscores = seg_scores[rseg];
                    for(auto lseg : rseg->m_left_connections) {
                        auto& lscores = seg_scores[lseg];
                        auto lsize = lscores.size();
                        TPathScore min_score;
                        if(lsize > 0)
                            min_score = *lscores.begin();
                        for(auto s : rscores) 
                            lscores.emplace(get<0>(s)+lseg->m_kmer_count, get<1>(s)+lseg->m_seq.size());
                        while(lscores.size() > max_stack)
                            lscores.erase(lscores.begin());
                        if(lscores.size() > lsize || *lscores.begin() > min_score)
                            next_step_segs.push(lseg);
                    }
                }
                std::swap(current_segs, next_step_segs);
            }

            //find top left ends
            map<TPathScore, list<GFAIterator>> ordered_lefts;
            for(auto it = begin(); it != end(); ++it) {
                if(it->m_left_connections.empty()) {
                    for(auto s : seg_scores[it]) 
                        ordered_lefts[s].push_back(it);
                }
            }
            while(ordered_lefts.size() > max_stack)
                ordered_lefts.erase(ordered_lefts.begin());

            //output variants
            stack<tuple<GFAIterator, TPathScore, GFAIterator>> remaining_segments;  // segment, currnet score, left connection
            //max score on top
            for(auto& sl : ordered_lefts) {
                for(auto iseg : sl.second)
                    remaining_segments.emplace(iseg, sl.first, end());
            }

            unsigned total = 0;
            list<GFAIterator> path;
            TPathScore path_score;
            set<string> seqs;
            while(!remaining_segments.empty() && total < max_variants) {
                GFAIterator iseg = get<0>(remaining_segments.top());
                TPathScore score = get<1>(remaining_segments.top());
                GFAIterator lseg = get<2>(remaining_segments.top());
                remaining_segments.pop();

                while(!path.empty() && path.back() != lseg)
                    path.pop_back();
                path.push_back(iseg);
                if(lseg == end())
                    path_score = score;

                if(iseg->m_right_connections.empty()) { // end
                    string seq;
                    for(auto i : path) {
                        for(auto& base : i->m_seq)
                            seq.push_back(base.m_nt);
                    }
                    if(seqs.insert(seq).second) {
                        out << ">" << Target() << ":" << front().m_group << ":" << ++total << ":" << get<0>(path_score);
                        for(auto i : path)
                            out << " " << i->m_num;
                        out << "\n" << seq << "\n";
                    }
                } else {                               // keep extending
                    TPathScore remaining_score(get<0>(score)-iseg->m_kmer_count, get<1>(score)-iseg->m_seq.size());
                    for(GFAIterator it : iseg->m_right_connections) {
                        if(seg_scores[it].count(remaining_score))
                            remaining_segments.emplace(it, remaining_score, iseg);
                    }
                }
            }
        }
        void PrintAllTranslatedVariants(ostream& out, unsigned max_variants, GeneticCode& genetic_code) {
            typedef tuple<size_t, int> TPathScore;
            unordered_map<GFAIterator, set<TPathScore>, SGFAIteratorHash> seg_scores;
            size_t max_stack = 2*max_variants; // give some extra in case of identical sequences

            //propagate score from right to left
            stack<GFAIterator> current_segs;
            for(GFAIterator it = begin(); it != end(); ++it) {
                if(it->m_right_connections.empty()) { // start from right extrems
                    current_segs.push(it);
                    TPathScore iscore(it->m_kmer_count, it->m_seq.size());
                    seg_scores[it].insert(iscore);
                }
            }
            while(!current_segs.empty()) {
                stack<GFAIterator> next_step_segs;
                while(!current_segs.empty()) {
                    GFAIterator rseg = current_segs.top();
                    current_segs.pop();
                    auto& rscores = seg_scores[rseg];
                    for(auto lseg : rseg->m_left_connections) {
                        auto& lscores = seg_scores[lseg];
                        auto lsize = lscores.size();
                        TPathScore min_score;
                        if(lsize > 0)
                            min_score = *lscores.begin();
                        for(auto s : rscores) 
                            lscores.emplace(get<0>(s)+lseg->m_kmer_count, get<1>(s)+lseg->m_seq.size());
                        while(lscores.size() > max_stack)
                            lscores.erase(lscores.begin());
                        if(lscores.size() > lsize || *lscores.begin() > min_score)
                            next_step_segs.push(lseg);
                    }
                }
                std::swap(current_segs, next_step_segs);
            }

            //find top left ends
            map<TPathScore, list<GFAIterator>> ordered_lefts;
            for(auto it = begin(); it != end(); ++it) {
                if(it->m_left_connections.empty()) {
                    for(auto s : seg_scores[it]) 
                        ordered_lefts[s].push_back(it);
                }
            }
            while(ordered_lefts.size() > max_variants)
                ordered_lefts.erase(ordered_lefts.begin());

            //output variants
            stack<tuple<GFAIterator, TPathScore, GFAIterator>> remaining_segments;  // segment, currnet score, left connection
            //max score on top
            for(auto& sl : ordered_lefts) {
                for(auto iseg : sl.second)
                    remaining_segments.emplace(iseg, sl.first, end());
            }

            unsigned total = 0;
            list<GFAIterator> path;
            TPathScore path_score;
            set<string> seqs;
            while(!remaining_segments.empty() && total < max_variants) {
                GFAIterator iseg = get<0>(remaining_segments.top());
                TPathScore score = get<1>(remaining_segments.top());
                GFAIterator lseg = get<2>(remaining_segments.top());
                remaining_segments.pop();

                while(!path.empty() && path.back() != lseg)
                    path.pop_back();
                path.push_back(iseg);
                if(lseg == end())
                    path_score = score;

                if(iseg->m_right_connections.empty()) { // end
                    string seq;
                    for(auto i : path) {
                        for(auto& base : i->m_seq)
                            seq.push_back(base.m_nt);
                    }
                    seq = genetic_code.Translate(seq, false);
                    if(seqs.insert(seq).second) {
                        out << ">" << Target() << ":" << front().m_group << ":" << ++total << ":" << get<0>(path_score);
                        for(auto i : path)
                            out << " " << i->m_num;
                        out << "\n" << seq << "\n";
                    }
                } else {                               // keep extending
                    TPathScore remaining_score(get<0>(score)-iseg->m_kmer_count, get<1>(score)-iseg->m_seq.size());
                    for(GFAIterator it : iseg->m_right_connections) {
                        if(seg_scores[it].count(remaining_score))
                            remaining_segments.emplace(it, remaining_score, iseg);
                    }
                }
            }
        }
        void PrintSelectedVariants(ostream& out) {
            typedef pair<GFAIterator, GFAIterator> TItP;
            struct STwoGFAIteratorsHash { size_t operator()(const TItP& p) const { return SGFAIteratorHash()(p.first)^SGFAIteratorHash()(p.second); } };

            map<string, tuple<int,int,size_t>> variants;
            GFAIteratorUSet included_segments;

            typedef tuple<size_t, int> TPathScore;
            unordered_map<TItP, TPathScore, STwoGFAIteratorsHash> included_edges;

            for(auto it = begin(); it != end(); ++it) {
                if(included_segments.insert(it).second) {
                    Path path = ExpandToMax(it);
                    TPathScore score(it->m_left_kmer_count+it->m_kmer_count+it->m_right_kmer_count, it->m_left_len+it->m_seq.size()+it->m_right_len);
                    variants[path.Sequence()] = make_tuple(it->m_num, it->m_num, get<0>(score));
                    included_segments.insert(path.m_segments.begin(), path.m_segments.end());
                    for(int i = 0; i < (int)path.m_segments.size()-1; ++i) {
                        TItP edge(path.m_segments[i], path.m_segments[i+1]);
                        included_edges[edge] = max(included_edges[edge], score);
                    }
                }
            }

            for(auto left = begin(); left != end(); ++left) {
                for(auto right : left->m_right_connections) {
                    TItP edge(left, right);
                    TPathScore score(left->m_left_kmer_count+left->m_kmer_count+right->m_kmer_count+right->m_right_kmer_count,
                                     left->m_left_len+left->m_seq.size()+right->m_seq.size()+right->m_right_len);
                    if(included_edges[edge] < score) {
                        included_edges[edge] = score;
                        Path path = ExpandEdgeToMax(left, right);
                        variants[path.Sequence()] = make_tuple(left->m_num, right->m_num, get<0>(score));
                        for(int i = 0; i < (int)path.m_segments.size()-1; ++i) {
                            TItP newedge(path.m_segments[i], path.m_segments[i+1]);
                            included_edges[newedge] = max(included_edges[newedge], score);
                        }
                    }
                }
            }

            for(auto& var : variants)
                out <<  ">Contig_" << Target() << ":" << front().m_group << ":" << get<0>(var.second) << ":" << get<1>(var.second) << ":" << get<2>(var.second) << "\n" << var.first << "\n";
        }

        const string& Target() const { return m_acc; }
        int KmerLen() const { return m_kmer_len; }

        void ScoreGraph(unordered_set<uint32_t>& target_words, int word_size) {
            GFAIterator best = end();
            for(auto it = begin(); it != end(); ++it) {
                if(it->m_left_connections.empty()) {
                    if(best == end() || it->m_kmer_count+it->m_right_kmer_count > best->m_kmer_count+best->m_right_kmer_count)
                        best = it;
                }
            }
            auto path = ExpandToMax(best);
            string consensus = path.Sequence();
	        
            for(int p = 0; p <= (int)consensus.size()-word_size; ++p) {
                string seed = consensus.substr(p, word_size);
                uint32_t word = 0;
                for(char c : seed) {
                    word = word << 2;
                    word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                }
                m_score += target_words.count(word);
            }                
        }    
    };


    class Spider;
    typedef list<Spider> TSpiderCollection;

    class Spider : public GFAGraph {
    private:

        set<Node> m_end_kmers;
        DBGraph* m_graphp;
        double m_fraction;

    public:
        Spider(const map<string, string>& contigs, DBGraph& graph, double fraction, const string& acc) : GFAGraph(acc, graph.KmerLen()), m_graphp(&graph), m_fraction(fraction) {
            int kmer_len = m_graphp->KmerLen();
            for(auto& contig : contigs) {
                auto& seq = contig.second;
                int len = seq.size();                
                Node rnode = m_graphp->GetNode(seq.substr(len-kmer_len));
                if(rnode.isValid())
                    m_end_kmers.insert(rnode);
                Node lnode = m_graphp->GetNode(seq.substr(0, kmer_len));
                if(lnode.isValid())
                    m_end_kmers.insert(lnode.ReverseComplement());                    
            }
        }
        Spider(const list<Node>& lkmers, const list<Node>& rkmers, DBGraph& graph, double fraction, const string& acc) : GFAGraph(acc, graph.KmerLen()), m_graphp(&graph), m_fraction(fraction) {
            for(auto& kmer : rkmers) // right ends of contigs
                m_end_kmers.insert(kmer);
            for(auto& kmer : lkmers) // left ends of contigs
                m_end_kmers.insert(kmer.ReverseComplement());
        }
        void DetectCycles() {
            for(auto& seg : *this) {
                if(seg.m_right_connections.empty() || seg.m_left_connections.empty())
                    continue;
                unordered_set<GFASegment*> visited;
                stack<GFASegment*> vaiting;
                vaiting.push(&seg);
                while(!vaiting.empty() && !seg.m_cyclical) {
                    auto& current = *vaiting.top();
                    visited.insert(vaiting.top());
                    vaiting.pop();
                    for(auto iseg : current.m_right_connections) {
                        if(&(*iseg) == &seg) {
                            seg.m_cyclical = true;
                            break;
                        }
                        if(visited.count(&(*iseg)))
                           continue;
                        vaiting.push(&(*iseg));
                    }
                }
            }            
        }
        void UpdateEndKmers() {
            int kmer_len = m_graphp->KmerLen();
            for(auto segmi = begin(); segmi != end(); ++segmi) {
                if(segmi->m_right_connections.empty()) {
                    int len = segmi->m_seq.size();
                    auto& rkmers = segmi->m_seq[len-1].m_right_kmers;
                    rkmers.clear();
                    list<Path> paths = Expand(segmi, len-1, kmer_len-1, 0);
                    for(auto& path : paths) {
                        if(path.Length() == kmer_len) {
                            auto node = m_graphp->GetNode(path.Sequence());
                            if(node.isValid() && find(rkmers.begin(), rkmers.end(), node) == rkmers.end())
                                rkmers.push_front(node);
                        }
                    }
                }
                if(segmi->m_left_connections.empty()) {
                    auto& lkmers = segmi->m_seq[0].m_left_kmers;
                    lkmers.clear();
                    list<Path> paths = Expand(segmi, 0, 0, kmer_len-1);
                    for(auto& path : paths) {
                        if(path.Length() == kmer_len) {
                            auto node = m_graphp->GetNode(path.Sequence());
                            if(node.isValid() && find(lkmers.begin(), lkmers.end(), node) == lkmers.end())
                                lkmers.push_front(node);
                        }
                    }
                }
            }

            set<Node> end_kmers;
            for(auto& seg : *this) {
                if(seg.m_left_connections.empty()) {
                    for(Node& node : seg.m_seq.front().m_left_kmers) {
                        if(m_end_kmers.count(node))
                            end_kmers.insert(node);
                    }
                }
                if(seg.m_right_connections.empty()) {
                    for(Node& node : seg.m_seq.back().m_right_kmers) {
                        if(m_end_kmers.count(node.ReverseComplement()))
                            end_kmers.insert(node.ReverseComplement());
                    }
                }
            }
            std::swap(end_kmers, m_end_kmers);
        }
        bool ResetEndKmers(const multimap<Node,string>& lacc, const multimap<Node,string>& racc) {
            bool reset = false;
            for(auto it_loop = m_end_kmers.begin(); it_loop != m_end_kmers.end(); ) {
                auto it = it_loop++;
                if(!racc.count(*it) && !lacc.count(it->ReverseComplement())) {
                    reset = true;
                    m_end_kmers.erase(it);
                }               
            }

            return reset;
        }
        static int EndsIntersect(const set<Node>& ends1, const set<Node>& ends2) {
            int minsize = min(ends1.size(), ends2.size());
            vector<Node> intersect(minsize);
            return (set_intersection(ends1.begin(), ends1.end(), ends2.begin(), ends2.end(), intersect.begin()) - intersect.begin());
        }
        bool EndsIncludedIn(const Spider& other) {
            return (EndsIntersect(EndKmers(), other.EndKmers()) == (int)EndKmers().size());
        }
        void ConnectOneEnd(const Node& node, int ext_len) {
            unordered_map<Node, Position, typename Node::Hash> links; // p >= 0 position of RIGHT kmer end; p < 0 len-|p| is position of LEFT kmer end
            list<GFAIterator> active_segms;

            string kmer = m_graphp->GetNodeSeq(node);
            GraphDigger graph_digger(*m_graphp, m_fraction, 0, 0, false);

            emplace_front();
            begin()->m_num = ++m_max_num;
            ++m_size;
            for(int k = 0; k < KmerLen()-1; ++k) {
                front().m_seq.emplace_back();
                front().m_seq.back().m_nt = kmer[k];
            } 
            front().m_seq.RightExtend(kmer.back(), node);
            front().m_seq.front().m_left_kmers.push_front(node);

            active_segms.push_front(begin());
            links[node] = Position(begin(), KmerLen()-1);
            for(int i = 0; i < ext_len && OneRightStep(links, active_segms, node, graph_digger); ++i) {}
            
            //remove forks from start/end kmers
            for(auto iseg = begin(); iseg != end(); ++iseg) {
                auto& seq = iseg->m_seq;                

                auto& lk = seq.front().m_left_kmers;
                if(iseg->m_left_connections.empty() && find(lk.begin(), lk.end(), node) != lk.end()) {
                    list<Path> paths = Expand(iseg, 0, 0, 2*(KmerLen()-1));  // inlet kmers are connected by the last base only
                        for(auto& path : paths) {
                            for(int i = 1; i < (int)path.m_segments.size(); ++i) {
                                auto isegl = path.m_segments[i-1];
                                auto isegr = path.m_segments[i];
                                for(auto ilc = isegr->m_left_connections.begin(); ilc != isegr->m_left_connections.end(); ) {
                                    auto lc = *(ilc++);
                                    if(lc != isegl)
                                        UnLinkSegments(lc, isegr);
                                }
                            }
                        }
                }
                               
                auto& rk = seq.back().m_right_kmers;
                if(iseg->m_right_connections.empty() && !rk.empty()) {
                    auto revnode = rk.front().ReverseComplement();
                    if(m_end_kmers.find(revnode) != m_end_kmers.end()) {
                        list<Path> paths = Expand(iseg, iseg->m_seq.size()-1, KmerLen()-1, 0);
                        for(auto& path : paths) {
                            for(int i = (int)path.m_segments.size()-2; i >= 0; --i) {
                                auto isegl = path.m_segments[i];
                                auto isegr = path.m_segments[i+1];
                                for(auto irc = isegl->m_right_connections.begin(); irc != isegl->m_right_connections.end(); ) {
                                    auto rc = *(irc++);
                                    if(rc != isegr)
                                        UnLinkSegments(isegl, rc); 
                                }
                            }
                        }
                    }
                }
                
            }                        
            
            /*
            MergeForks();
            //prevent inclusion of large contig chunks (>=kmer) in loops
            for(auto iseg = begin(); iseg != end(); ++iseg) {
                auto& seq = iseg->m_seq;
                auto& lk = seq[0].m_left_kmers;
                if(iseg->m_left_connections.empty() && find(lk.begin(), lk.end(), node) != lk.end()) {
                    if(seq.size() == 1 && iseg->RightSingle()) {
                        auto inext = iseg->m_right_connections.front();
                        if(inext == iseg)
                            break;
                        seq += inext->m_seq;
                        TransferRightLinks(inext, iseg);
                        RemoveLinksToSegment(inext);
                        RemoveSegment(inext);
                    }
                    break;
                }
            }
            */                       

            RemoveLooseEnds();
            // keep forks for cleaning           MergeForks(); 
            //there is only one group which will be 0            AssignGroupNumber();
            GenerateKmers(*m_graphp);
        }
        void CalculateDistanceToStartingEnds() {
            for(auto& seg : *this) {
                seg.m_left_len = numeric_limits<int>::max();
                seg.m_right_len = numeric_limits<int>::max();
                if(seg.m_left_connections.empty()) {
                    for(Node& node : seg.m_seq.front().m_left_kmers) {
                        if(m_end_kmers.count(node))
                            seg.m_left_len = 0;
                    }
                }
                if(seg.m_right_connections.empty()) {
                    for(Node& node : seg.m_seq.back().m_right_kmers) {
                        if(m_end_kmers.count(node.ReverseComplement()))
                            seg.m_right_len = 0;
                    }
                }
            }

            bool keep_doing = true;
            while(keep_doing) {
                keep_doing = false;
                for(auto& seg : *this) {
                    int min_left = numeric_limits<int>::max();
                    for(auto& lc : seg.m_left_connections) {
                        if(lc->m_left_len != numeric_limits<int>::max())
                            min_left = min(min_left, lc->m_left_len+(int)lc->m_seq.size()); 
                    }                                          
                    if(min_left < seg.m_left_len) {
                        seg.m_left_len = min_left;
                        keep_doing = true;
                    }
                        
                    int min_right = numeric_limits<int>::max();
                    for(auto& rc : seg.m_right_connections) {
                        if(rc->m_right_len != numeric_limits<int>::max())
                            min_right = min(min_right, rc->m_right_len+(int)rc->m_seq.size());
                    }                    
                    if(min_right < seg.m_right_len) {
                        seg.m_right_len = min_right;
                        keep_doing = true;
                    }
                }
            }
        }
        bool RemoveLooseEnds() {
            bool rslt = false;
            CalculateDistanceToStartingEnds();

            for(auto it_loop = begin(); it_loop != end(); ) {
                auto it = it_loop++;
                if(it->m_left_len == numeric_limits<int>::max() || it->m_right_len == numeric_limits<int>::max()) {
                    RemoveSegment(it);
                    rslt = true;
                }
            }

            return rslt;
        }
        void MergeRedundantDuplicates() {            
            GFAGraph::MergeRedundantDuplicates();
            CalculateDistanceToStartingEnds();

            for(auto segi = begin(); segi != end(); ++segi) {
                if(!segi->m_marked_for_erase && (segi->m_left_len == numeric_limits<int>::max() || segi->m_right_len == numeric_limits<int>::max())) {
                    RemoveLinksToSegment(segi);
                    segi->m_left_connections.clear();
                    segi->m_right_connections.clear();
                    segi->m_marked_for_erase = true;
                }
            }
        }
        void MergeForks() {
            //mask ends form merging
            for(auto& segm : *this) {
                if(segm.m_left_connections.empty()) {
                    for(Node& node : segm.m_seq.front().m_left_kmers) {
                        if(m_end_kmers.count(node))
                            segm.m_seq.front().m_nt = 0; 
                    }
                }
                if(segm.m_right_connections.empty()) {
                    for(Node& node : segm.m_seq.back().m_right_kmers) {
                        if(m_end_kmers.count(node.ReverseComplement()))
                            segm.m_seq.back().m_nt = 0;
                    }
                }
            }
            GFAGraph::MergeForks();
            //restore ends
            for(auto& segm : *this) {
                if(segm.m_seq.front().m_nt == 0 && !segm.m_seq.front().m_left_kmers.empty()) { //empty check for 1bp segments
                    Node& node = segm.m_seq.front().m_left_kmers.front();
                    string kmer = m_graphp->GetNodeSeq(node);
                    segm.m_seq.front().m_nt = kmer.front();
                }
                if(segm.m_seq.back().m_nt == 0 && !segm.m_seq.back().m_right_kmers.empty()) { //empty check for 1bp segments
                    Node& node = segm.m_seq.back().m_right_kmers.front();
                    string kmer = m_graphp->GetNodeSeq(node);
                    segm.m_seq.back().m_nt = kmer.back();
                }
            }            
        }
        bool RemoveHair(DBGraph& dbg_graph, double eps) { // Spider doesn't have loose ends
            return false;
        }
        void RemoveRedundantPaths(size_t maxp) {}         // Spider doesn't have loose ends
        bool Absorb(Spider& other) {
            map<Node, GFAIterator> other_lends;
            map<Node, GFAIterator> other_rends;
            for(auto iseg = other.begin(); iseg != other.end(); ++iseg) {
                if(iseg->m_left_connections.empty()) {
                    for(Node& node : iseg->m_seq.front().m_left_kmers) {
                        if(other.m_end_kmers.count(node))
                            other_lends[node] = iseg;
                    }
                }
                if(iseg->m_right_connections.empty()) {
                    for(Node& node : iseg->m_seq.back().m_right_kmers) {
                        if(other.m_end_kmers.count(node.ReverseComplement()))
                            other_rends[node] = iseg;
                    }
                }
            }

            int same = 0;
            int opposite = 0;
            for(auto iseg = begin(); iseg != end(); ++iseg) {                
                if(iseg->m_left_connections.empty()) {
                    for(Node& node : iseg->m_seq.front().m_left_kmers) {
                        if(other_lends.count(node))
                            ++same;
                        else if(other_rends.count(node.ReverseComplement()))
                            ++opposite;
                    }
                }
                if(iseg->m_right_connections.empty()) {
                    for(Node& node : iseg->m_seq.back().m_right_kmers) {
                        if(other_rends.count(node))
                            ++same;
                        else if(other_lends.count(node.ReverseComplement()))
                            ++opposite;
                    }
                }
            }
            if(same == 0 && opposite == 0)
                return false;

            if(opposite > same)
                ReverseComplement();

            int group = front().m_group;
            for(auto i = begin(); i != end(); ++i) {
                if(i->m_left_connections.empty()) {
                    auto iseg = i;
                    auto base = iseg->m_seq.front();
                    for(Node& node : base.m_left_kmers) {
                        auto rslt = other_lends.find(node);
                        if(rslt != other_lends.end()) {
                            if(iseg->m_seq.size() > 1) {
                                emplace_front();
                                front().m_group = group;
                                front().m_seq.LeftExtend(base.m_nt, node);
                                iseg->m_seq.pop_front();
                                LinkSegments(begin(), iseg);
                                iseg = begin();
                            }
                            auto iother = rslt->second;
                            if(iother->m_seq.size() > 1) {
                                iother->m_seq.pop_front();
                                LinkSegments(iseg, iother);                        
                            } else {
                                TransferRightLinks(iother, iseg);
                                other.RemoveSegment(iother);
                            }
                            break;
                        }
                    }
                }

                if(i->m_right_connections.empty()) {
                    auto iseg = i;
                    auto base = iseg->m_seq.back();
                    for(Node& node : base.m_right_kmers) {
                        auto rslt = other_rends.find(node);
                        if(rslt != other_rends.end()) {
                            if(iseg->m_seq.size() > 1) {
                                emplace_front();
                                front().m_group = group;
                                front().m_seq.RightExtend(base.m_nt, node);
                                iseg->m_seq.pop_back();
                                LinkSegments(iseg, begin());
                                iseg = begin();
                            }
                            auto iother = rslt->second;
                            if(iother->m_seq.size() > 1) {
                                iother->m_seq.pop_back();
                                LinkSegments(iother, iseg);                        
                            } else {
                                TransferLeftLinks(iother, iseg);
                                other.RemoveSegment(iother);
                            }
                        }
                    }
                }
            }

            for(auto iseg = other.begin(); iseg != other.end(); ++iseg)
                iseg->m_group = group;
            splice(end(), other);
            m_end_kmers.insert(other.m_end_kmers.begin(), other.m_end_kmers.end());

            MergeForks();
            EnumerateSegments();

            return true;
        }
        TSpiderCollection SplitGroups() { // destroys source   
            TSpiderCollection splitted;
            map<int, TSpiderCollection::iterator> group_to_graph;

            map<int, list<Node>> lkmers; // left ends of contigs
            map<int, list<Node>> rkmers; // right ends of contigs
            for(auto& segm : *this) {
                if(segm.m_left_connections.empty()) {
                    for(Node& node : segm.m_seq.front().m_left_kmers) {
                        if(m_end_kmers.count(node))
                            rkmers[segm.m_group].push_back(node);
                    }
                }
                if(segm.m_right_connections.empty()) {
                    for(Node& node : segm.m_seq.back().m_right_kmers) {
                        if(m_end_kmers.count(node.ReverseComplement()))
                            lkmers[segm.m_group].push_back(node);
                    }
                }
            }

            while(!this->empty()) {
                int group = this->front().m_group;
                if(!group_to_graph.count(group)) {
                    splitted.emplace_front(lkmers[group], rkmers[group], *m_graphp, m_fraction, Target());
                    group_to_graph[group] = splitted.begin();
                }
                auto& dest = *group_to_graph[group];
                dest.splice(dest.end(), *this, this->begin());
            }

            for(auto& graph : splitted)
                graph.Size() = graph.size();

            return splitted;
        }
        void DeleteLEnd(const Node& lkmer) {
            m_end_kmers.erase(lkmer.ReverseComplement());
        }
        void DeleteREnd(const Node& rkmer) {
            m_end_kmers.erase(rkmer);
        }
        const set<Node>& EndKmers() const { return m_end_kmers; }
        int Connections() const { return m_end_kmers.size(); }
        
    private:
        bool OneRightStep(unordered_map<Node, Position, typename Node::Hash>& links, list<GFAIterator>& active_segms, const Node& initial_node, GraphDigger& graph_digger) {
            if(active_segms.empty()) 
                return false;

            //extend active segments by 1bp
            for(auto i_loop = active_segms.begin(); i_loop != active_segms.end(); ) {
                auto iseg = i_loop++; // iterator to iterator
                auto& segm = **iseg;
                const Node& node = segm.m_seq.back().m_right_kmers.front();
                auto successors = graph_digger.GetReversibleNodeSuccessors(node);

                if(successors.empty()) {
                    active_segms.erase(iseg);
                } else if(successors.size() == 1) {
                    segm.m_seq.RightExtend(successors[0].m_nt, successors[0].m_node);
                } else {
                    (*iseg)->m_seq.back().m_fork |= eRightFork;
                    for(auto& suc : successors)
                        BranchToRight(suc, *iseg, active_segms);
                    active_segms.erase(iseg);                          //remove from active
                }
            }

            //stop path if a start is reached
            for(auto i_loop = active_segms.begin(); i_loop != active_segms.end(); ) {
                auto iseg = i_loop++; // iterator to iterator
                auto& segm = **iseg;
                const Node& node = segm.m_seq.back().m_right_kmers.front();
                if(node == initial_node)
                    active_segms.erase(iseg);
            }

            //mark left forks
            for(auto iseg = active_segms.begin(); iseg != active_segms.end(); ++iseg) {
                auto& segm = **iseg;
                const Node& node = segm.m_seq.back().m_right_kmers.front();
                auto revnode = node.ReverseComplement();
                //left forks
                if(graph_digger.GetReversibleNodeSuccessors(revnode).size() > 1) {
                    int right_kmer_end = segm.m_seq.size()-1;
                    int left_kmer_end = right_kmer_end-KmerLen()+1;
                    if(left_kmer_end >= 0) {                          // still same segment
                        segm.m_seq[left_kmer_end].m_fork |= eLeftFork;
                    } else {                                          //short segment need to expand (multiple matches possible)
                        list<Path> paths = Expand(*iseg, right_kmer_end, KmerLen()-1, 0);
                        string kmer = m_graphp->GetNodeSeq(segm.m_seq.back().m_right_kmers.front());
                        for(auto& path : paths) {
                            if(path.Sequence() == kmer) {
                                GFAIterator destp = path.m_segments.front();
                                left_kmer_end = path.m_left;
                                destp->m_seq[left_kmer_end].m_fork |= eLeftFork;
                            }
                        }
                    }
                }
            }

            //remember kmer positions and collapse identical kmers
            for(auto i_loop = active_segms.begin(); i_loop != active_segms.end(); ) {
                auto iseg = i_loop++; // iterator to iterator
                GFAIterator segp = *iseg;
                Node& node = segp->m_seq.back().m_right_kmers.front();
                Position right_kmer_end(segp, segp->m_seq.size()-1);
                auto rslt = links.emplace(node, right_kmer_end);                          //remember kmer for collapsing
                if(!rslt.second) {                                                        //kmer exists - collaps
                    Position dest = rslt.first->second;
                    if(dest.m_pos < 0)                                                    //negative destination segment - find right kmer end
                        dest = FindRightKmerEnd(node, dest);
                    GFAIterator destp = dest.m_segmp;
                    int pos = dest.m_pos;
                    //kmers are connected by the last right base; we hope that MergeForls will collaps the rest
                    CutOffLeft(destp, pos, links);
                    segp->m_seq.pop_back();                                               //remove last base; may result in empty segment       
                    if(segp->m_seq.empty()) {
                        TransferLeftLinks(segp, destp);
                        RemoveSegment(segp);
                    } else {
                        LinkSegments(segp, destp);
                    }
                    active_segms.erase(iseg);
                }
            }            

            //check if ends are reached (no self loop)
            for(auto i_loop = active_segms.begin(); i_loop != active_segms.end(); ) {
                auto iseg = i_loop++;
                auto& segm = **iseg;
                const Node& node = segm.m_seq.back().m_right_kmers.front();
                auto revnode = node.ReverseComplement();
                if(revnode != initial_node && m_end_kmers.find(revnode) != m_end_kmers.end())
                    active_segms.erase(iseg);
            }            

            return !active_segms.empty();
        }
        void BranchToRight(const Successor& suc, GFAIterator segp, list<GFAIterator>& active_segms) {
            emplace_front();                                 //insert and extend new segment
            begin()->m_num = ++m_max_num;
            ++m_size;
            front().m_seq.RightExtend(suc.m_nt, suc.m_node);
            active_segms.push_front(begin());              //put new segment in active list
            LinkSegments(segp, begin());                    //graph link
        }
        void BranchToLeft(const Successor& suc, GFAIterator segp, list<GFAIterator>& active_segms) {
            emplace_front();                                 //insert and extend new segment
            begin()->m_num = ++m_max_num;
            ++m_size;
            front().m_seq.LeftExtend(Complement(suc.m_nt), suc.m_node.ReverseComplement());
            active_segms.push_front(begin());              //put new segment in active list
            LinkSegments(begin(),segp);                     //graph link
        }
        Position FindRightKmerEnd(Node& node, const Position& left_kmer_pos) { //left_kmer_pos negative; returns normal right
            GFAIterator destp = left_kmer_pos.m_segmp;
            int left_kmer_end = destp->m_seq.size()+left_kmer_pos.m_pos;
            Position pos(destp, left_kmer_end+KmerLen()-1);
            if(pos.m_pos >= (int)destp->m_seq.size()) { //short segment need to expand          
                list<Path> paths = Expand(destp, left_kmer_end, 0, KmerLen()-1);
                string kmer = m_graphp->GetNodeSeq(node);
                for(auto& path : paths) {
                    if(path.Sequence() == kmer) {
                        pos.m_segmp = path.m_segments.back();
                        pos.m_pos = path.m_right;
                        break;
                    }
                }
            }
            return pos;
        }
        Position FindLeftKmerEnd(Node& node, const Position& right_kmer_pos) { //right_kmer_pos positive; returns negative left
            GFAIterator destp = right_kmer_pos.m_segmp;
            int right_kmer_end = right_kmer_pos.m_pos;
            int left_kmer_end = right_kmer_end-KmerLen()+1;
            Position pos(destp, left_kmer_end-(int)destp->m_seq.size());
            if(left_kmer_end < 0) {  //short segment need to expand
                list<Path> paths = Expand(destp, right_kmer_end, KmerLen()-1, 0);
                string kmer = m_graphp->GetNodeSeq(node);
                for(auto& path : paths) {
                    if(path.Sequence() == kmer) {
                        pos.m_segmp = path.m_segments.front();
                        pos.m_pos = path.m_left-(int)pos.m_segmp->m_seq.size();
                        break;
                    }
                } 
            }
            return pos;
        }
        void CutOffLeft(GFAIterator destp, int pos, unordered_map<Node, Position, typename Node::Hash>& links) { //keeps positions >= pos in the original segment; moves previous positions to a new segment
            if(pos == 0) //nothing to cut
                return;

            //create and link new segment for left part 
            emplace_front();
            auto leftp = begin();
            leftp->m_num = ++m_max_num;
            ++m_size;
            TransferLeftLinks(destp, leftp);
            LinkSegments(leftp, destp);

            //move sequence before pos      
            copy(destp->m_seq.begin(), destp->m_seq.begin()+pos, back_inserter(leftp->m_seq));
            destp->m_seq.erase(destp->m_seq.begin(), destp->m_seq.begin()+pos);

            UpdateLinks(leftp, false, destp, true, links);
        }
        void CutOffRight(GFAIterator destp, int pos, unordered_map<Node, Position, typename Node::Hash>& links) { //keeps positions <= pos in the original segment; moves next positions to a new segment
            if(pos == (int)destp->m_seq.size()-1) //nothing to cut
                return;

            //create and link new segment for right part 
            emplace_front();
            auto rightp = begin();
            rightp->m_num = ++m_max_num;
            ++m_size;
            TransferRightLinks(destp, rightp);
            LinkSegments(destp, rightp);

            //move sequence after pos      
            copy(destp->m_seq.begin()+pos+1, destp->m_seq.end(), back_inserter(rightp->m_seq));
            destp->m_seq.erase(destp->m_seq.begin()+pos+1, destp->m_seq.end());

            UpdateLinks(destp, true, rightp, false, links);
        }
        void UpdateLinks(GFAIterator leftp, bool check_left, GFAIterator rightp, bool check_right, unordered_map<Node, Position, typename Node::Hash>& links) {
            int left_len = leftp->m_seq.size();
            int right_len = rightp->m_seq.size();
            int orig_len = left_len+right_len;

            if(check_left) {
                auto& lkmers = leftp->m_seq.front().m_left_kmers;
                if(!lkmers.empty()) {
                    auto firstb = links.find(lkmers.front());
                    if(firstb != links.end() && firstb->second.m_segmp == leftp && firstb->second.m_pos == -orig_len) //change first one only if belongs to leftp     
                        firstb->second.m_pos = -left_len;
                }
            }
            for(int p = 0; p < left_len; ++p) {
                auto& rkmers = leftp->m_seq[p].m_right_kmers;
                if(!rkmers.empty())
                    links[rkmers.front()] = Position(leftp, p); 
                auto& lkmers = leftp->m_seq[p].m_left_kmers;
                if((p > 0 || !check_left) && !lkmers.empty())
                    links[lkmers.front()] = Position(leftp, p-left_len);
            }
            
            for(int p = 0; p < right_len; ++p) {
                auto& rkmers = rightp->m_seq[p].m_right_kmers;
                if((p < right_len-1 || !check_right) && !rkmers.empty())
                    links[rkmers.front()] = Position(rightp, p); 
                auto& lkmers = rightp->m_seq[p].m_left_kmers;
                if(!lkmers.empty())
                    links[lkmers.front()] = Position(rightp, p-right_len);
            }
            if(check_right) {
                auto& rkmers = rightp->m_seq.back().m_right_kmers;
                if(!rkmers.empty()) {
                    auto lastb = links.find(rkmers.front());
                    if(lastb != links.end() && lastb->second.m_segmp == rightp && lastb->second.m_pos == orig_len-1) //change last one only if belongs to rightp          
                        lastb->second.m_pos = rightp->m_seq.size()-1;
                }
            }
        }
    };


    template<class Collection>
    void EnumerateCollection(Collection& collection) {
        map<string, int> target_groups;
        for(auto& graph : collection) {
            graph.EnumerateSegments();
            int group = ++target_groups[graph.Target()];
            for(auto& seg : graph)
                seg.m_group = group;
        }
    } 

    template<class Collection>
    void SortCollection(Collection& collection) {
        typedef typename Collection::value_type TGraphType;
        collection.sort([](const TGraphType& a, const TGraphType& b)
                        {
                            if(a.Target() != b.Target()) {
                                return a.Target() < b.Target();
                            } else if(a.Score() != b.Score()) {
                                return a.Score() > b.Score();
                            } else if(a != b) {
                                return a < b;             //lexigraphical order of segments
                            } else {                      // all segmets are same at this poit - check connection order
                                auto ia = a.begin();
                                auto ib = b.begin();
                                for( ; ia != a.end(); ++ia, ++ib) {
                                    int ar = distance(ia->m_right_connections.begin(), ia->m_right_connections.end());
                                    int br = distance(ib->m_right_connections.begin(), ib->m_right_connections.end());
                                    if(ar != br)
                                        return ar < br;
                                    auto iar = ia->m_right_connections.begin();
                                    auto ibr = ib->m_right_connections.begin();
                                    for( ; iar != ia->m_right_connections.end(); ++iar, ++ibr) {
                                        if((*iar)->m_seq != (*ibr)->m_seq)
                                            return (*iar)->m_seq < (*ibr)->m_seq;
                                    }
                                    
                                    int al = distance(ia->m_left_connections.begin(), ia->m_left_connections.end());
                                    int bl = distance(ib->m_left_connections.begin(), ib->m_left_connections.end());
                                    if(al != bl)
                                        return al < bl;
                                    auto ial = ia->m_left_connections.begin();
                                    auto ibl = ib->m_left_connections.begin();
                                    for( ; ial != ia->m_left_connections.end(); ++ial, ++ibl) {
                                        if((*ial)->m_seq != (*ibl)->m_seq)
                                            return (*ial)->m_seq < (*ibl)->m_seq;
                                    }                                         
                                }

                                return false;
                            }
                        });
    }

    template<class Collection>
    void RemoveRedundantGraphs(Collection& collection) {
        SortCollection(collection);

        unordered_map<CKmerHashCount::Index, list<typename Collection::iterator>, CKmerHashCount::Index::Hash> kmer_to_graph;
        for(auto it = collection.begin(); it != collection.end(); ++it) {
            for(auto& kmer : it->KSignature())
                kmer_to_graph[kmer].push_back(it); 
        }

        list<typename Collection::iterator> erased;
        set<GFAGraph*> erasedp;
        for(auto it = collection.begin(); it != collection.end(); ++it) {
            for(auto jt : kmer_to_graph[it->KSignature().front()]) {
                if(it == jt || erasedp.count(&(*jt))) 
                    continue;
                if(jt->KSignature().size() == it->KSignature().size() && it->Score() > jt->Score())
                    continue;
                if(it->IsSubGraphOf(*jt)) {
                    cerr << jt->Target() << ":" << jt->front().m_group << " " << jt->KSignature().size() << " kills " << it->Target() << ":" << it->front().m_group << " " << it->KSignature().size() << endl;
                    erased.push_back(it);
                    erasedp.insert(&(*it));
                    break;
                }
            }
        }
        for(auto it : erased)
            collection.erase(it);

        EnumerateCollection(collection);
    }

    void RemoveSpiderSubGraphs(TSpiderCollection& collection) {
        EnumerateCollection(collection);

        SortCollection(collection);
        unordered_map<CKmerHashCount::Index, list<TSpiderCollection::iterator>, CKmerHashCount::Index::Hash> kmer_to_graph;
        for(auto it = collection.begin(); it != collection.end(); ++it) {
            for(auto& kmer : it->KSignature())
                kmer_to_graph[kmer].push_back(it); 
        }
        struct IHash { size_t operator()(TSpiderCollection::iterator it) const { return std::hash<void*>()(&(*it)); } };
        unordered_set<TSpiderCollection::iterator, IHash> erased;
        for(auto it = collection.begin(); it != collection.end(); ++it) {
            for(auto jt : kmer_to_graph[it->KSignature().front()]) {
                if(it == jt || erased.count(jt)) 
                    continue;
                if(it->IsSubGraphOf(*jt) && it->EndsIncludedIn(*jt)) {
                    cerr << jt->Target() << ":" << jt->front().m_group << " " << jt->KSignature().size() << " kills " << it->Target() << ":" << it->front().m_group << " " << it->KSignature().size() << endl;
                    erased.insert(it);
                    break;
                }
            }            
        }
        for(auto it : erased)
            collection.erase(it);

    }
 
        
    template <typename E>
    class CBlockedIndex {
    public:
        CBlockedIndex(size_t block_size = 0) : m_block_size(block_size), m_centinel(0) {}
        void Reset(size_t block_size) {
            m_storage.clear();
            m_block_size = block_size;
            m_centinel = 0;
        }
        void GetBlock(E*& blockp, size_t& indexb) {
            while(!m_centinel.Set(1, 0));    //grab deque       
            m_storage.emplace_back();
            m_storage.back().resize(m_block_size);
            blockp = m_storage.back().data();
            indexb = (m_storage.size()-1)*m_block_size;
            m_centinel = 0;                 //release deque     
        }
        size_t BlockSize() const { return m_block_size; }
        E& operator[](size_t index) { return m_storage[index/m_block_size][index%m_block_size]; }
        const E& operator[](size_t index) const { return m_storage[index/m_block_size][index%m_block_size]; }
        
    private:
        deque<vector<E>> m_storage;
        size_t m_block_size;
        SAtomic<uint8_t> m_centinel;
    };

    template<class Collection>
    class GraphCleaner {
    private:
        typedef typename Collection::value_type Graph;
        Collection& m_gfa_collection;
        TTargets* m_targetsp;
        DBGraph& m_graph;
        int m_kmer_len;
        double m_fraction;
        GraphDigger m_graph_digger;
        list<array<CReadHolder,2>>& m_raw_reads;
        double m_entropy_level;
        bool m_no_reads;
        bool m_no_pairs;
        int m_not_aligned_len;
        int m_not_aligned_count;
        int m_aligned_count;
        int m_maxp;
        int m_ncores;

        int m_read_length;
        int m_insert_length = 0;
        int m_insert_max = 0;
        int m_insert_min = 0;
        TReadPos m_read_pos;
        int m_graph_chunks = 8;
        list<array<deque<uint8_t>,2>> m_read_colors;        
        CBlockedIndex<CReadHolder::string_iterator> m_selected_reads_index;
        size_t m_selected_reads_shift = 11;           // first bit direction (0 is plus); next 10 bits for position in read; others number of read
        size_t m_selected_reads_mask = (1ULL << m_selected_reads_shift) - 1;
        size_t m_max_allowed_read_length = 1ULL << (m_selected_reads_shift-1);
        void IndexToRead(uint64_t index, bool& reversed, int& pos, CReadHolder::string_iterator& is) const {
            reversed = index&1;
            pos = (index&m_selected_reads_mask) >> 1;
            is = m_selected_reads_index[index >> m_selected_reads_shift];
        }
        CReadHolder::string_iterator IndexToRead(uint64_t index) const { return m_selected_reads_index[index >> m_selected_reads_shift]; }
        SAtomic<size_t> m_kmers_in_assembly = 0;
        mutex m_out_mutex;
        
    public:
        GraphCleaner(Collection& gfa_collection, TTargets* targetsp, DBGraph& graph, double fraction, double entropy_level, int not_aligned_len, int not_aligned_count, int aligned_count, int maxp, 
                     bool no_reads, bool no_pairs, list<array<CReadHolder,2>>& raw_reads, int ncores) : 
            m_gfa_collection(gfa_collection), m_targetsp(targetsp), m_graph(graph), m_kmer_len(graph.KmerLen()), m_fraction(fraction), m_graph_digger(graph, fraction, 0, 0, false), m_raw_reads(raw_reads), 
            m_entropy_level(entropy_level), m_no_reads(no_reads), m_no_pairs(no_pairs), 
            m_not_aligned_len(not_aligned_len), m_not_aligned_count(not_aligned_count), m_aligned_count(aligned_count), m_maxp(maxp), m_ncores(ncores) {

            if(m_no_reads && m_no_pairs)
                return;

            EstimateReads();
            
            m_read_pos = TReadPos(m_kmer_len);

            CStopWatch timer;
            timer.Restart();
            ColorKmers();
            cerr << "Kmers colored in: " << timer.Elapsed();

            timer.Restart();
            ClipReads();
            cerr << "Reads clipped in: " << timer.Elapsed(); 

            {
                Collection filtered_graphs;
                Collection raw_graphs;
                raw_graphs.splice(raw_graphs.end(), m_gfa_collection);

                while(!raw_graphs.empty()) {
                    uint8_t color_id = raw_graphs.front().ID();
                    auto first = raw_graphs.begin();
                    auto last = first;
                    size_t i = 0;
                    for( ; last != raw_graphs.end() && last->ID() == color_id; ++last, ++i);
                    m_gfa_collection.splice(m_gfa_collection.end(), raw_graphs, first, last);
                    cerr << "Filtering " << i << " graphs" << endl;

                    timer.Restart();
                    InitKmerHash();
                    cerr << "Init kmer hash in " << timer.Elapsed(); 

                    timer.Restart();
                    size_t selected_reads = IndexReads();
                    cerr << "Selected reads: " << selected_reads << " indexed in: " << timer.Elapsed(); 
                    cerr << "Read hash elements: " << m_kmers_in_assembly.Load() << " Read hash size: " << m_read_pos.TableSize() << endl; 
                    size_t total = 0;
                    for(auto it = m_read_pos.Begin(); it != m_read_pos.End(); ++it)
                        total += get<0>(*it.GetMapped()).size();
                    cerr << "Total entries: " << total << endl;                
            
                    timer.Restart();
                    AlignReads();
                    cerr << "Align reads in " << timer.Elapsed(); 

                    filtered_graphs.splice(filtered_graphs.end(), m_gfa_collection);
                }
                m_gfa_collection.splice(m_gfa_collection.end(), filtered_graphs);
            } 
        }
            

    private:
        void EstimateReads() {
            double total_seq = 0;
            size_t total_reads = 0;
            size_t mates = 0;
            for(auto& reads : m_raw_reads) {
                total_seq += reads[0].TotalSeq()+reads[1].TotalSeq();
                total_reads += reads[0].ReadNum()+reads[1].ReadNum();
                mates += reads[0].ReadNum();
            }
            m_read_length = 0;
            if(total_reads > 0)
                m_read_length = total_seq/total_reads+0.5;
            if(mates > 0) {
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
                    list<array<CReadHolder,2>> connected_mate_pairs = m_graph_digger.ConnectPairs(mate_pairs, long_insert_size, m_ncores, false);
                    CReadHolder connected_mates(false);
                    for(auto& mp : connected_mate_pairs) {
                        for(CReadHolder::string_iterator is = mp[0].sbegin(); is != mp[0].send(); ++is)
                            connected_mates.PushBack(is);
                    }

                    m_insert_length = connected_mates.N50();
                    m_insert_max = connected_mates.NXX(0.15);
                    m_insert_min = connected_mates.NXX(0.85);
                }
            }

            cerr << "Read length: " << m_read_length << " Inserts: " << m_insert_max << "/" << m_insert_length << "/" << m_insert_min << endl;
        }

        void ColorKmers() {
            int total_graphs = m_gfa_collection.size();
            int graph_chunk_size = max(2*m_ncores, total_graphs/m_graph_chunks+1);
            int i = 0;
            for(auto& gfa_graph : m_gfa_collection) {
                gfa_graph.ID() = i++/graph_chunk_size;
                gfa_graph.Sentinel() = 0;
            }

            list<function<void()>> jobs;            
            for(int thr = 0; thr < m_ncores; ++thr) {
                jobs.push_back(bind(&GraphCleaner::ColorKmersJob, this));
            }            
            RunThreads(m_ncores, jobs);

            for(auto& gfa_graph : m_gfa_collection)
                gfa_graph.Sentinel() = 0;
        }
        void ColorKmersJob() {
            for(auto& gfa_graph : m_gfa_collection) {
                if(!gfa_graph.Sentinel().Set(1, 0))
                    continue;

                uint8_t color = (1 << gfa_graph.ID());
                for(auto& index : gfa_graph.KSignature()) {
                    Node node(index, Node::ePlus);
                    m_graph.SetColor(node, color);
                }
            }
        }

        void ClipReads() {
            int64_t total_reads = 0;
            int64_t total_seq = 0;
            for(const auto& reads : m_raw_reads) {
                total_reads += reads[0].ReadNum()+reads[1].ReadNum(); 
                total_seq += reads[0].TotalSeq()+reads[1].TotalSeq();
            }

            list<function<void()>> jobs;
            for(auto& job_input : m_raw_reads) {
                m_read_colors.emplace_back();
                jobs.push_back(bind(&GraphCleaner::ClipReadsJob, this, ref(job_input), ref(m_read_colors.back())));  
            }      
            
            RunThreads(m_ncores, jobs);

            int64_t total_reads_after = 0;
            int64_t total_seq_after = 0;
            for(const auto& reads : m_raw_reads) {
                total_reads_after += reads[0].ReadNum()+reads[1].ReadNum(); 
                total_seq_after+= reads[0].TotalSeq()+reads[1].TotalSeq();
            }
            cerr << "Raw reads: " << total_reads << " Raw sequence: " << total_seq << " Reads after clip: " << total_reads_after << " Sequence after clip: " << total_seq_after << endl;
        }
        void ClipReadsJob(array<CReadHolder,2>& reads, array<deque<uint8_t>,2>& colors) {
            array<CReadHolder,2> clipped_reads{CReadHolder(true), CReadHolder(false)};
            for(auto is = reads[0].sbegin(); is != reads[0].send(); ++is) {
                auto is1 = is;
                auto is2 = ++is;
                string read1 = *is1;
                uint8_t color1 = m_graph_digger.CheckAndClipReadLite(read1);
                string read2 = *is2;
                uint8_t color2 = m_graph_digger.CheckAndClipReadLite(read2);
                if(color1 || color2) {
                    clipped_reads[0].PushBack(read1);
                    colors[0].push_back(color1);
                    clipped_reads[0].PushBack(read2);
                    colors[0].push_back(color2);
                }
            }
            for(CReadHolder::string_iterator is = reads[1].sbegin(); is != reads[1].send(); ++is) {
                string read = *is;
                uint8_t color = m_graph_digger.CheckAndClipReadLite(read);
                if(color) {
                    clipped_reads[1].PushBack(read);
                    colors[1].push_back(color);
                }
            }
            reads[0].Swap(clipped_reads[0]);
            reads[1].Swap(clipped_reads[1]);
        }

        void InitKmerHash() {
            size_t kmer_estimate = 0;
            for(auto& gfa_graph : m_gfa_collection) {
                for(auto& seg : gfa_graph)
                    kmer_estimate += seg.m_seq.size();  
            }
            cerr << "Kmer estimate: " << kmer_estimate << endl;

            m_read_pos.Reset(5*kmer_estimate, m_ncores);
        
            list<function<void()>> jobs;            
            for(int thr = 0; thr < m_ncores; ++thr) {
                jobs.push_back(bind(&GraphCleaner::KmerHashJob, this));
            }            
            RunThreads(m_ncores, jobs);

            for(auto& gfa_graph : m_gfa_collection)
                gfa_graph.Sentinel() = 0;
        }

        void KmerHashJob() {
            for(auto& gfa_graph : m_gfa_collection) {
                if(!gfa_graph.Sentinel().Set(1, 0))
                    continue;

                // create entries for kmers in graphs   
                for(auto& seg : gfa_graph) {
                    for(auto& base : seg.m_seq) {
                        for(auto& node : base.m_left_kmers) {
                            string kmer_seq = m_graph.GetNodeSeq(node);
                            if(Entropy(kmer_seq.begin(), m_kmer_len) <= m_entropy_level)
                                continue;
                            auto kmer = m_graph.GetNodeKmer(node.DropStrand());
                            if(m_read_pos.Find(kmer) == nullptr) {
                                m_read_pos.FindOrInsert(kmer);
                                m_read_pos.FindOrInsert(revcomp(kmer, m_kmer_len));
                            }
                        }
                        for(auto& node : base.m_right_kmers) {
                            string kmer_seq = m_graph.GetNodeSeq(node);
                            if(Entropy(kmer_seq.begin(), m_kmer_len) <= m_entropy_level)
                                continue;
                            auto kmer = m_graph.GetNodeKmer(node.DropStrand());
                            if(m_read_pos.Find(kmer) == nullptr) {
                                m_read_pos.FindOrInsert(kmer);
                                m_read_pos.FindOrInsert(revcomp(kmer, m_kmer_len));
                            }
                        }
                    }
                }
            }
        }
                
        size_t IndexReads() {
            size_t block = 1000;
            atomic<size_t> total_reads(0);
            list<function<void()>> jobs;
            auto ic = m_read_colors.begin();
            for(auto& job_input : m_raw_reads) 
                jobs.push_back(bind(&GraphCleaner::IndexReadsJob, this, ref(job_input), ref(*ic++), ref(total_reads)));
            m_selected_reads_index.Reset(block);
            m_kmers_in_assembly = 0;
            
            RunThreads(m_ncores, jobs);

            return total_reads;
        }
        void IndexReadsJob(array<CReadHolder,2>& clipped_reads, array<deque<uint8_t>,2>& colors, atomic<size_t>& included_reads) {
            auto Index = [this](int dir, int pos, size_t index) { return dir+(pos << 1) + (index << m_selected_reads_shift); };
            size_t total_kmers = 0;
            size_t total_reads = 0;
            size_t indexb = 0;
            uint8_t color = (1 << m_gfa_collection.front().ID());
            CReadHolder::string_iterator* blockp = nullptr;
            m_selected_reads_index.GetBlock(blockp, indexb);
            size_t index = indexb;
            for(int p = 0; p < 2; ++p) {
                auto is = clipped_reads[p].sbegin();
                auto ic = colors[p].begin();
                for( ; is != clipped_reads[p].send(); ++is, ++ic) {
                    if(!(*ic&color))
                        continue;
                    int rlen = is.ReadLen();
                    if(rlen > (int)m_max_allowed_read_length)
                        throw runtime_error("Read longer "+to_string(m_max_allowed_read_length));                
                    if(rlen >= m_kmer_len) {
                        bool included = false;
                        int pos = rlen-m_kmer_len;
                        for(CReadHolder::kmer_iterator ik = is.KmersForRead(m_kmer_len); pos >= 0; ++ik, --pos) { // iteration from last kmer to first  
                            TKmer kmer = *ik;
                            auto rslt = m_read_pos.Find(kmer);
                            if(rslt == nullptr)
                                continue;
                            
                            while(!get<2>(*rslt).Set(1, 0));    //grab info 
                            if(get<0>(*rslt).empty())
                                ++total_kmers;
                            get<0>(*rslt).push_back(Index(0, pos, index));
                            if(pos == rlen-m_kmer_len)
                                get<1>(*rslt).push_back(Index(0, pos, index));
                            get<2>(*rslt) = 0;                  // release info 

                            rslt = m_read_pos.Find(revcomp(kmer, m_kmer_len));
                            while(!get<2>(*rslt).Set(1, 0));    //grab info 
                            if(get<0>(*rslt).empty())
                                ++total_kmers;
                            get<0>(*rslt).push_back(Index(1, rlen-m_kmer_len-pos, index));
                            get<2>(*rslt) = 0;                  // release info 

                            included = true;
                        }
                        if(included) {
                            ++total_reads;
                            blockp[index-indexb] = is;
                            if(++index == indexb+m_selected_reads_index.BlockSize()) {
                                m_selected_reads_index.GetBlock(blockp, indexb);
                                index = indexb;
                            }                        
                        }
                    }
                }
            }
            m_kmers_in_assembly.m_atomic += total_kmers;
            included_reads += total_reads;
        }

        enum EReadSupport { eSupported, eNotSupported, eLowCoverage };
        pair<EReadSupport, int> LengthSupportedByReads(Path& path, bool toright) {
 
            int slen = path.Length();
           
            map<int, tuple<int, int>> counts; // length before fork, count for aligned, count for not aligned  
            if(toright) {
                int len = path.m_segments.front()->m_seq.size()-path.m_left;
                counts.emplace(len, make_tuple(0, 0));
                for(int i = 1; i < (int)path.m_segments.size()-1; ++i) {
                    auto& s = path.m_segments[i]->m_seq; 
                    len += s.size();
                    if(s.back().m_fork&eRightFork)
                        counts.emplace(len, make_tuple(0, 0));
                }
            } else {
                int len = path.m_right+1;
                counts.emplace(len, make_tuple(0, 0));
                for(int i = (int)path.m_segments.size()-2; i > 0; --i) {
                    auto& s = path.m_segments[i]->m_seq;
                    len += s.size();
                    if(s.front().m_fork&eLeftFork)
                        counts.emplace(len, make_tuple(0, 0));
                }
            }
            if(counts.size() < 2)
                return make_pair(eSupported, slen);

            string seq = path.Sequence();
            CReadHolder rh(false);
            rh.PushBack(seq);
            vector<uint64_t> seqb((2*slen+63)/64);
            rh.sbegin().TrueBSeq(0, 0, !toright, seqb.data());

            typedef tuple<int, uint64_t> THitInfo;
            vector<THitInfo> hitsp;

            int extra_len = min(m_not_aligned_len/2, m_kmer_len-2);
            string kseq;
            int shift = 0;
            if(toright) {
                int segl = path.m_segments.front()->m_seq.size();
                shift = min(extra_len, segl-1);
                for(int i = segl-shift-1; i < segl-1; ++i)
                    kseq.push_back(path.m_segments.front()->m_seq[i].m_nt);
                kseq += seq.substr(0, m_kmer_len-shift);
            } else {
                int segl = path.m_segments.back()->m_seq.size();
                shift = min(extra_len, segl-1);
                kseq = seq.substr(slen-m_kmer_len+shift);
                for(int i = 1; i < shift+1; ++i)
                    kseq.push_back(path.m_segments.back()->m_seq[i].m_nt);
                ReverseComplementSeq(kseq.begin(), kseq.end());
            }
            TKmer kmer(kseq);

            auto rslt = m_read_pos.Find(kmer);
            if(rslt != nullptr) {
                TReadPosInfo& hits = get<0>(*rslt);                
                hitsp.reserve(hits.size());
                for(uint64_t index : hits) {
                    bool reversed; 
                    int rpos;
                    CReadHolder::string_iterator is;
                    IndexToRead(index, reversed, rpos, is);
                    int rlen = is.ReadLen();
                    int possible_extend = rlen-rpos-shift;
                    hitsp.emplace_back(possible_extend, index);
                }
            }

            int chunk_for_sort = 25;
            auto sorted_so_far = hitsp.begin();
            for(auto it = hitsp.begin(); it != hitsp.end(); ++it) {
                if(it == sorted_so_far) {
                    sorted_so_far += min(chunk_for_sort, int(hitsp.end()-sorted_so_far));
                    std::partial_sort(it, sorted_so_far, hitsp.end(), [](const THitInfo& a, const THitInfo& b) { return get<0>(a) > get<0>(b); });
                }

                int possible_extend = get<0>(*it);
                if(possible_extend < counts.begin()->first)
                    break;

                uint64_t index = get<1>(*it);
                bool reversed; 
                int rpos;
                CReadHolder::string_iterator is;
                IndexToRead(index, reversed, rpos, is);
                int rlen = is.ReadLen();
                
                int left_clip = rlen-possible_extend;
                int right_clip = max(0, possible_extend-slen);
                if(reversed)
                    std::swap(left_clip, right_clip);
                rlen -= left_clip+right_clip;

                vector<uint64_t> readb((2*rlen+63)/64);
                is.TrueBSeq(left_clip, right_clip, reversed, readb.data());
                int extend = min(rlen, (int)CReadHolder::string_iterator::CommomSeqLen(readb.data(), seqb.data(), readb.size()));

                for(auto ic = counts.begin(); ic != counts.end() && ic->first <= extend; ++ic) {
                    if(ic->first == extend) {  // exact touch
                        if(possible_extend-extend > m_not_aligned_len)
                            ++get<1>(ic->second);
                    } else {                   //cross
                        auto inext = next(ic);
                        int lim = inext == counts.end() ? slen-1 : inext->first-1;
                        if(extend > min(lim, ic->first+extra_len))
                           ++get<0>(ic->second); 
                    }
                }
                while(!counts.empty() && get<0>(counts.begin()->second) >= m_aligned_count)
                    counts.erase(counts.begin());
                if(counts.empty())
                    return make_pair(eSupported, slen);
            }

            for(auto& count : counts) {
                if(get<1>(count.second) >= m_not_aligned_count)
                    return make_pair(eNotSupported, count.first);
            }

            return make_pair(eLowCoverage, counts.begin()->first);
        }
        pair<EReadSupport, int> LengthSupportedByPairs(Path& path, Graph& gfa_graph, bool toright) {
            int slen = path.Length();
           
            map<int, tuple<int, int>> counts; // length before fork, count for aligned, count for not aligned  
            if(toright) {
                int len = path.m_segments.front()->m_seq.size()-path.m_left;
                counts.emplace(len, make_tuple(0, 0));
                for(int i = 1; i < (int)path.m_segments.size()-1; ++i) {
                    auto& s = path.m_segments[i]->m_seq; 
                    len += s.size();
                    if(s.back().m_fork&eRightFork)
                        counts.emplace(len, make_tuple(0, 0));
                }
            } else {
                int len = path.m_right+1;
                counts.emplace(len, make_tuple(0, 0));
                for(int i = (int)path.m_segments.size()-2; i > 0; --i) {
                    auto& s = path.m_segments[i]->m_seq;
                    len += s.size();
                    if(s.front().m_fork&eLeftFork)
                        counts.emplace(len, make_tuple(0, 0));
                }
            }
            if(counts.size() < 2)
                return make_pair(eSupported, slen);

            string seq = path.Sequence();
            CReadHolder rh(false);
            rh.PushBack(path.Sequence());
            vector<uint64_t> seqb((2*slen+63)/64);
            rh.sbegin().TrueBSeq(0, 0, !toright, seqb.data());

            // kmer hash for finding second mate position   
            int klen_mate = 21; // must be <= 32
            unordered_map<uint64_t, list<int>> seq_pos;
            seq_pos.reserve(slen);
            {
                uint64_t kmer = seqb[0] & ((1ULL << 2*klen_mate) - 1);              // first kmer   
                int spos = 0;
                seq_pos[kmer].push_back(spos);
                for(int last_bit_pair = 2*klen_mate; last_bit_pair < 2*slen; last_bit_pair += 2) {
                    kmer >>= 2;                                                      // drop leftmost nt    
                    uint64_t lbp = (seqb[last_bit_pair/64] >> last_bit_pair%64) & 3; // extract new rightmost nt    
                    kmer |= (lbp << 2*(klen_mate-1));                                // shift new nt and put in place 
                    seq_pos[kmer].push_back(++spos);
                }
            }

            int extra_len = min(m_not_aligned_len/2, m_kmer_len-2);
            int ex_shift = 0;
            deque<pair<int, TKmer>> starting_kmers;
            {
                list<string> starting_paths;
                if(toright) {
                    int segl = path.m_segments.front()->m_seq.size();
                    ex_shift = min(extra_len, segl-1);
                    string rend = seq.substr(1, m_kmer_len-1-ex_shift);
                    list<Path> extra_paths = gfa_graph.Expand(path.m_segments.front(), path.m_segments.front()->m_seq.size()-1, m_kmer_len-2, 0);
                    for(auto& path : extra_paths) 
                        starting_paths.push_back(path.Sequence()+rend);
                } else {
                    int segl = path.m_segments.back()->m_seq.size();
                    ex_shift = min(extra_len, segl-1);
                    string rend = seq.substr(slen-m_kmer_len+ex_shift, m_kmer_len-1-ex_shift);
                    ReverseComplementSeq(rend.begin(), rend.end());
                    list<Path> extra_paths = gfa_graph.Expand(path.m_segments.back(), 0, 0, m_kmer_len-2);
                    for(auto& path : extra_paths) {
                        starting_paths.push_back(path.Sequence());
                        ReverseComplementSeq(starting_paths.back().begin(), starting_paths.back().end());
                        starting_paths.back() += rend;
                    }
                }
                for(string& sp : starting_paths) {
                    CReadHolder rh(false);
                    rh.PushBack(sp);
                    int shift = ex_shift;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len); ik != rh.kend(); ++ik, ++shift) // iteration from last kmer to first 
                        starting_kmers.emplace_back(shift, *ik);
                }
            }
            std::sort(starting_kmers.begin(), starting_kmers.end());
            starting_kmers.erase(std::unique(starting_kmers.begin(), starting_kmers.end()), starting_kmers.end());

            typedef tuple<int, int, bool, uint64_t> THitInfo;  // possible extend, first mate clip, second mate exist, hit index    
            deque<THitInfo> hitsp;
            for(auto& shiftk : starting_kmers) {
                int shift = shiftk.first;
                TKmer& kmer = shiftk.second;
                auto rslt = m_read_pos.Find(kmer);
                if(rslt == nullptr)
                    continue;
                TReadPosInfo& hits = (shift == ex_shift) ? get<0>(*rslt) : get<1>(*rslt);
               
                for(uint64_t index : hits) {
                    bool reversed; 
                    int rpos;
                    CReadHolder::string_iterator it;
                    IndexToRead(index, reversed, rpos, it);
                    if(!it.HasMate() || reversed)
                        continue;

                    rpos += shift;
                    int rlen = it.ReadLen();
                    int extend = rlen-rpos;

                    auto itm = it.GetMate();
                    int mlen = itm.ReadLen();
                    if(mlen < klen_mate) {
                        hitsp.emplace_back(extend, rpos, false, index);
                        continue;
                    }

                    uint64_t mkmera = 0;
                    itm.TrueBSeq(mlen-klen_mate, 0, true, &mkmera);
                                        
                    auto rslt = seq_pos.find(mkmera);
                    int mextend = 0;
                    if(rslt != seq_pos.end()) {
                        for(int mp : rslt->second) {
                            mextend = mp+mlen;
                            if(mextend > extend)
                                hitsp.emplace_back(mextend, rpos, true, index);                            
                        }
                    }
                    if(mextend <= extend)
                        hitsp.emplace_back(extend, rpos, false, index);
                }        
            }
                
            int chunk_for_sort = 25;
            auto sorted_so_far = hitsp.begin();
            for(auto it = hitsp.begin(); it != hitsp.end(); ++it) {
                if(it == sorted_so_far) {
                    sorted_so_far += min(chunk_for_sort, int(hitsp.end()-sorted_so_far));
                    std::partial_sort(it, sorted_so_far, hitsp.end(), [](const THitInfo& a, const THitInfo& b) { return get<0>(a) > get<0>(b); });
                }

                int possible_extend = get<0>(*it);
                if(possible_extend < counts.begin()->first)
                    break;
                uint64_t index = get<3>(*it);
                CReadHolder::string_iterator is = IndexToRead(index);

                {
                    int rlen = is.ReadLen();
                    int left_clip = get<1>(*it);
                    int right_clip = max(0, rlen-left_clip-slen);
                    int remaining_len = rlen-left_clip;
                    rlen -= left_clip+right_clip;
            
                    vector<uint64_t> readb((2*rlen+63)/64);
                    is.TrueBSeq(left_clip, right_clip, false, readb.data());
                    int extend = min(rlen, (int)CReadHolder::string_iterator::CommomSeqLen(readb.data(), seqb.data(), readb.size()));

                    for(auto ic = counts.begin(); ic != counts.end() && ic->first <= extend; ++ic) {
                        if(ic->first == extend) {  // exact touch
                            if(remaining_len-extend > m_not_aligned_len)
                                ++get<1>(ic->second);
                        } else {                   //cross
                            auto inext = next(ic);
                            int lim = inext == counts.end() ? slen-1 : inext->first-1;
                            if(extend > min(lim, ic->first+extra_len))
                                ++get<0>(ic->second); 
                        }
                    }
                    while(!counts.empty() && get<0>(counts.begin()->second) >= m_aligned_count)
                        counts.erase(counts.begin());
                    if(counts.empty())
                        return make_pair(eSupported, slen);

                    if(extend < rlen || !get<2>(*it))
                        continue;
                }
            
                {
                    auto itm = is.GetMate();
                    int mlen = itm.ReadLen();
                    int mleft = possible_extend-mlen;

                    int mextend = 0;
                    int first_chunk_len = 0;
                    int mshift = mleft%32;
                    if(mshift > 0) {
                        first_chunk_len = min(mlen, 32-mshift);
                        uint64_t mfirst_chunk = 0;
                        itm.TrueBSeq(mlen-first_chunk_len, 0, true, &mfirst_chunk);
                        uint64_t sfirst_chunk = seqb[mleft/32] >> 2*mshift;
                        mextend = min(first_chunk_len, (int)CReadHolder::string_iterator::CommomSeqLen(&mfirst_chunk, &sfirst_chunk, 1));
                    }                
                    if(mextend == first_chunk_len) {
                        int mright_clip = first_chunk_len;
                        int mleft_clip = max(0, possible_extend-slen);
                        mlen -= mleft_clip+mright_clip;
                        if(mlen > 0) {
                            vector<uint64_t> mateb((2*mlen+63)/64);
                            itm.TrueBSeq(mleft_clip, mright_clip, true, mateb.data());
                            mextend += min(mlen, (int)CReadHolder::string_iterator::CommomSeqLen(mateb.data(), seqb.data()+(mleft+first_chunk_len)/32, mateb.size()));
                        }
                    }                
                    mextend += mleft;                 

                    auto ic_loop = counts.upper_bound(mleft+1); // first fork which could be affected
                    for( ; ic_loop != counts.end() && ic_loop->first <= mextend; ) {
                        auto ic = ic_loop++;
                        if(ic->first == mextend) {  // exact touch
                            if(possible_extend-mextend > m_not_aligned_len)
                                ++get<1>(ic->second);
                        } else {                   //cross
                            auto inext = next(ic);
                            int lim = inext == counts.end() ? slen-1 : inext->first-1;
                            if(mextend > min(lim, ic->first+extra_len)) {
                                if(++get<0>(ic->second) >= m_aligned_count)
                                    counts.erase(ic);
                            }
                        }
                    }
                    if(counts.empty())
                        return make_pair(eSupported, slen);
                }                   
            }

            for(auto& count : counts) {
                if(get<1>(count.second) >= m_not_aligned_count)
                    return make_pair(eNotSupported, count.first);                
            }

            return make_pair(eLowCoverage, counts.begin()->first);       
        }
        void AlignReads() {
            vector<Collection> thread_rslts(m_ncores);
            list<function<void()>> jobs;
            for(int thr = 0; thr < m_ncores; ++thr)
                jobs.push_back(bind(&GraphCleaner::AlignReadsJob, this, ref(thread_rslts[thr])));
            RunThreads(m_ncores, jobs);

            m_gfa_collection.clear();
            for(auto& trslt : thread_rslts)
                m_gfa_collection.splice(m_gfa_collection.end(), trslt);
            
            EnumerateCollection(m_gfa_collection);
        }
        void AlignReadsJob(Collection& rslts) { 
            size_t cutoff_for_break = 15;

            for(auto& gfa_graph : m_gfa_collection) {
                if(!gfa_graph.Sentinel().Set(1, 0))
                    continue;

                CStopWatch timer;
                timer.Restart();
                auto& acc = gfa_graph.Target();
                int group = gfa_graph.front().m_group;
                int min_len = 0;
                if(m_targetsp != nullptr)
                    min_len = get<1>((*m_targetsp)[acc]);
            
                for(auto it = gfa_graph.begin(); it != gfa_graph.end(); ++it) {
                    auto& seg = *it;
                    int seglen = seg.m_seq.size();
                        
                    if(seg.RightConnectionsNum() == 2) {
                        int shift = 0;
                        if(seglen > 1) {
                            char lastc = seg.m_seq.back().m_nt;
                            for(auto rc : seg.m_right_connections) {
                                if(lastc == rc->m_seq.front().m_nt) {                                
                                    while(++shift < min(seglen,m_kmer_len-1)-1 && seg.m_seq[seglen-1-shift].m_nt == lastc);
                                    break;
                                }
                            }
                        }

                        list<Path> paths = gfa_graph.Expand(it, it->m_seq.size()-1, shift, m_kmer_len-1-shift);
                        int minpath = numeric_limits<int>::max();
                        for(auto& path : paths)
                            minpath = min(minpath, path.Length());
                        if(minpath < m_kmer_len && seglen-shift-1 >= m_kmer_len-minpath) {
                            shift += m_kmer_len-minpath;
                            paths = gfa_graph.Expand(it, it->m_seq.size()-1, shift, m_kmer_len-1-shift);
                        }

                        if(paths.size() == 2 && paths.front().Length() == m_kmer_len && paths.back().Length() == m_kmer_len) {
                            bool asamepair = false;
                            auto& apath = paths.front();
                            TKmer akmer(apath.Sequence());                     
                            auto rslt = m_read_pos.Find(akmer);
                            if(rslt != nullptr) {
                                TReadPosInfo& hits = get<0>(*rslt);                  
                                if(hits.size() == 2) {
                                    auto is = IndexToRead(hits[0]);
                                    asamepair = is.HasMate() && is.GetMate() == IndexToRead(hits[1]);
                                }
                            }   

                            bool bsamepair = false;
                            auto& bpath = paths.back();
                            TKmer bkmer(bpath.Sequence());                     
                            rslt = m_read_pos.Find(bkmer);
                            if(rslt != nullptr) {
                                TReadPosInfo& hits = get<0>(*rslt);       
                                if(hits.size() == 2) {
                                    auto is = IndexToRead(hits[0]);
                                    bsamepair = is.HasMate() && is.GetMate() == IndexToRead(hits[1]);
                                }
                            }

                            Path* pathp = nullptr;
                            if(asamepair && !bsamepair)
                                pathp = &apath;
                            else if(!asamepair && bsamepair)
                                pathp = &bpath;

                            if(pathp != nullptr) {
                                bool simple_chain = true;
                                for(int j = 1; j < (int)pathp->m_segments.size()-1 && simple_chain; ++j)
                                    simple_chain = (pathp->m_segments[j]->LeftConnectionsNum() == 1);
                                if(simple_chain) {
                                    if(pathp->m_segments.size() > 2) {
                                        for(int j = 1; j < (int)pathp->m_segments.size()-1; ++j) {
                                            auto jsegp = pathp->m_segments[j];
                                            gfa_graph.RemoveSegment(jsegp);
                                        }
                                    } else {
                                        pathp->m_segments.front()->m_right_connections.remove(pathp->m_segments.back());
                                        pathp->m_segments.back()->m_left_connections.remove(pathp->m_segments.front());
                                    }
                                }
                            } 
                        }
                    }
                    if(seg.LeftConnectionsNum() == 2) {
                        int shift = 0;
                        if(seglen > 1) {
                            char firstc = seg.m_seq.front().m_nt;
                            for(auto rc : seg.m_left_connections) {
                                if(firstc == rc->m_seq.back().m_nt) {                                
                                    while(++shift < min(seglen,m_kmer_len-1)-1 && seg.m_seq[shift].m_nt == firstc);
                                    break;
                                }
                            }
                        }

                        list<Path> paths = gfa_graph.Expand(it, 0, m_kmer_len-1-shift, shift);
                        int minpath = numeric_limits<int>::max();
                        for(auto& path : paths)
                            minpath = min(minpath, path.Length());
                        if(minpath < m_kmer_len && seglen-shift-1 >= m_kmer_len-minpath) {
                            shift += m_kmer_len-minpath;
                            paths = gfa_graph.Expand(it, 0, m_kmer_len-1-shift, shift);
                        }

                        if(paths.size() == 2 && paths.front().Length() == m_kmer_len && paths.back().Length() == m_kmer_len) {
                            bool asamepair = false;
                            auto& apath = paths.front();
                            TKmer akmer(apath.Sequence());                     
                            auto rslt = m_read_pos.Find(akmer);
                            if(rslt != nullptr) {
                                TReadPosInfo& hits = get<0>(*rslt);                  
                                if(hits.size() == 2) {
                                    auto is = IndexToRead(hits[0]);
                                    asamepair = is.HasMate() && is.GetMate() == IndexToRead(hits[1]);
                                }
                            }   

                            bool bsamepair = false;
                            auto& bpath = paths.back();
                            TKmer bkmer(bpath.Sequence());                     
                            rslt = m_read_pos.Find(bkmer);
                            if(rslt != nullptr) {
                                TReadPosInfo& hits = get<0>(*rslt);       
                                if(hits.size() == 2) {
                                    auto is = IndexToRead(hits[0]);
                                    bsamepair = is.HasMate() && is.GetMate() == IndexToRead(hits[1]);
                                }
                            }

                            Path* pathp = nullptr;
                            if(asamepair && !bsamepair)
                                pathp = &apath;
                            else if(!asamepair && bsamepair)
                                pathp = &bpath;

                            if(pathp != nullptr) {
                                bool simple_chain = true;
                                for(int j = 1; j < (int)pathp->m_segments.size()-1 && simple_chain; ++j)
                                    simple_chain = (pathp->m_segments[j]->RightConnectionsNum() == 1);
                                if(simple_chain) {
                                    if(pathp->m_segments.size() > 2) {
                                        for(int j = 1; j < (int)pathp->m_segments.size()-1; ++j) {
                                            auto jsegp = pathp->m_segments[j];
                                            gfa_graph.RemoveSegment(jsegp);
                                        }           
                                    } else {
                                        pathp->m_segments.front()->m_right_connections.remove(pathp->m_segments.back());
                                        pathp->m_segments.back()->m_left_connections.remove(pathp->m_segments.front());
                                    }
                                }
                            } 
                        }
                    }
                } 

                gfa_graph.CutToChunks();
                Graph graph_copy(gfa_graph);
                size_t initial_gsize = gfa_graph.Size();
                size_t last_gsize = initial_gsize;
                TCopyInfo copies;  // get<0> - not checked, get<1> - checked    
                copies.reserve(gfa_graph.Size()+1);
                {
                    lock_guard<mutex> guard(m_out_mutex);
                    cerr << "Started graph: " << acc << ":" << group << " " << gfa_graph.Size() << endl;
                }

                bool interrupted = false;
                if(!m_no_reads) {
                    for(auto it = gfa_graph.begin(); it != gfa_graph.end(); ) {
                        if((it->m_left_check && it->m_right_check) || it->m_marked_for_erase) {
                            ++it;
                            continue;
                        }

                        list<GFAIterator> segments;
                        tuple<list<GFAIterator>,list<GFAIterator>>* copyp = nullptr;
                        auto copy_ofi = gfa_graph.end();
                        if(it->m_copy_of != nullptr) {
                            copy_ofi = it->m_copy_ofi;
                            copyp = &copies[it->m_copy_ofi];
                        } else if(copies.count(it)) {
                            copy_ofi = it;
                            copyp = &copies[it];
                        }
                        if(copyp != nullptr) {
                            segments.splice(segments.end(), get<0>(*copyp));
                            segments.push_front(copy_ofi);
                        } else {
                            segments.push_back(it);
                        }

                        for(int toright = 0; toright < 2; ++toright) {
                            unordered_set<Path, Path::Hash> accepted_paths;           // supported or low coverage      
                            unordered_map<Path, int, Path::Hash> not_supported_paths;
                            for(GFAIterator segp : segments) {                            
                                GFASegment& seg = *segp;
                                if(seg.m_marked_for_erase)
                                    continue;
                                bool& check = toright ? seg.m_right_check : seg.m_left_check;
                                if(toright)
                                    check = !(seg.m_seq.back().m_fork&eLeftBranch);
                                else
                                    check = !(seg.m_seq.front().m_fork&eRightBranch);

                                while(!check) {
                                    check = true;
                                    int extend = m_read_length;
                                    list<Path> paths = gfa_graph.Expand(segp, toright ? segp->m_seq.size()-1 : 0, toright ? 0 : extend, toright ? extend : 0, m_maxp, true);
                                    for(auto& path : paths) {
                                        int path_len = path.Length();
                                        if(path_len < m_kmer_len || accepted_paths.count(path))
                                            continue;                                       
                                        if(!path.IntactPath()) {
                                            check = false;
                                            continue;
                                        }

                                        int max_extend = -1;
                                        auto known_not_supported = not_supported_paths.find(path);
                                        if(known_not_supported != not_supported_paths.end())
                                            max_extend = known_not_supported->second;
                                        if(max_extend < 0) {
                                            auto rslt = LengthSupportedByReads(path, toright);
                                            if(rslt.first == eNotSupported) {
                                                max_extend = rslt.second;
                                                not_supported_paths[path] = max_extend;
                                            } else {
                                                accepted_paths.insert(path);
                                            }
                                        }

                                        if(max_extend >= 0) {
                                            int clip = path_len-max_extend-1;
                                            if(clip > 0) {
                                                if(toright)
                                                    path.ClipRight(clip);
                                                else
                                                    path.ClipLeft(clip);
                                            } 
                                            gfa_graph.RemovePath(path, toright, copies);
                                            accepted_paths.insert(path);
                                        }
                                    }
                                }
                            }
                        }

                        ++it;    //cleaning loop inserts segments at the list end

                        if(copyp != nullptr) {
                            segments.pop_front();
                            auto& checked = get<1>(*copyp);
                            checked.splice(checked.end(), segments);
                        }

                        if(gfa_graph.Size() > 2*last_gsize) {
                            it = gfa_graph.ReduceGraph(min_len, copies, it);
                            last_gsize = gfa_graph.Size();
                        }

                        if(gfa_graph.Size() > cutoff_for_break*initial_gsize) {
                            interrupted = true;
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted graphA: " << acc << ":" << group << " " << gfa_graph.Size() << endl;
                            break;
                        }
                    }
                    gfa_graph.ReduceGraph(min_len, copies, gfa_graph.end());
                    last_gsize = gfa_graph.Size();
                }

                if(interrupted) {
                    gfa_graph = graph_copy;
                } else if(!m_no_pairs && m_insert_length > 0) {
                    graph_copy = gfa_graph;
                                                                    
                    for(auto& seg : gfa_graph) {
                        seg.m_right_check = false;
                        seg.m_left_check = false;
                    }
                    for(auto& copy : copies) 
                        get<0>(copy.second).splice(get<0>(copy.second).end(), get<1>(copy.second));

                    for(auto it = gfa_graph.begin(); it != gfa_graph.end(); ) {
                        if((it->m_left_check && it->m_right_check) || it->m_marked_for_erase) {
                            ++it;
                            continue;
                        }

                        list<GFAIterator> segments;
                        tuple<list<GFAIterator>,list<GFAIterator>>* copyp = nullptr;
                        auto copy_ofi = gfa_graph.end();
                        if(it->m_copy_of != nullptr) {
                            copy_ofi = it->m_copy_ofi;
                            copyp = &copies[it->m_copy_ofi];
                        } else if(copies.count(it)) {
                            copy_ofi = it;
                            copyp = &copies[it];
                        }
                        if(copyp != nullptr) {
                            segments.splice(segments.end(), get<0>(*copyp));
                            segments.push_front(copy_ofi);
                        } else {
                            segments.push_back(it);
                        }

                        for(int toright = 0; toright < 2; ++toright) {
                            unordered_set<Path, Path::Hash> accepted_paths;           // supported or low coverage      
                            unordered_map<Path, int, Path::Hash> not_supported_paths;
                            for(GFAIterator segp : segments) {                            
                                GFASegment& seg = *segp;
                                if(seg.m_marked_for_erase)
                                    continue;
                                bool& check = toright ? seg.m_right_check : seg.m_left_check;
                                if(toright)
                                    check = !(seg.m_seq.back().m_fork&eLeftBranch);
                                else
                                    check = !(seg.m_seq.front().m_fork&eRightBranch);
                                
                                int extend = m_insert_max;
                                while(!check) {
                                    check = true;
                                    list<Path> paths = gfa_graph.Expand(segp, toright ? segp->m_seq.size()-1 : 0, toright ? 0 : extend, toright ? extend : 0, m_maxp, true);
                                    for(auto& path : paths) {
                                        int path_len = path.Length();
                                        if(path_len < m_kmer_len || accepted_paths.count(path))                                            
                                            continue;                                        
                                        
                                        if(!path.IntactPath()) {
                                            check = false;
                                            continue;
                                        }
                                        
                                        int max_extend = -1;
                                        auto known_not_supported = not_supported_paths.find(path);
                                        if(known_not_supported != not_supported_paths.end())
                                            max_extend = known_not_supported->second;
                                        if(max_extend < 0) {
                                            auto rslt = LengthSupportedByPairs(path, gfa_graph, toright);
                                            if(rslt.first == eNotSupported) { 
                                                max_extend = rslt.second;
                                                not_supported_paths[path] = max_extend;
                                            } else {
                                                accepted_paths.insert(path);
                                            }
                                        }
                                        if(max_extend >= 0) {
                                            int clip = path_len-max_extend-1;
                                            if(clip > 0) {
                                                if(toright)
                                                    path.ClipRight(clip);
                                                else
                                                    path.ClipLeft(clip);
                                            } 
                                            gfa_graph.RemovePath(path, toright, copies);
                                            accepted_paths.insert(path);
                                        }
                                    }
                                }
                            }
                        }

                        ++it;    //cleaning loop inserts segments at the list end

                        if(copyp != nullptr) {
                            segments.pop_front();
                            auto& checked = get<1>(*copyp);
                            checked.splice(checked.end(), segments);
                        }

                        if(gfa_graph.Size() > 2*last_gsize) {
                            it = gfa_graph.ReduceGraph(min_len, copies, it);
                            last_gsize = gfa_graph.Size();
                        }
                        
                        if(gfa_graph.Size() > cutoff_for_break*initial_gsize) {
                            lock_guard<mutex> guard(m_out_mutex);
                            cerr << "Interrupted graphB: " << acc << ":" << group << " " << gfa_graph.Size() << endl;
                            swap(gfa_graph, graph_copy);
                            break;
                        }
                    }
                    gfa_graph.ReduceGraph(min_len, copies, gfa_graph.end());
                }

                gfa_graph.RemoveRedundantPaths(m_maxp);
 
                size_t gsize = gfa_graph.Size();
                gfa_graph.MergeForks();                       // m_copy_of is not valid after this
                /*
                if(gfa_graph.RemoveHair(m_graph, m_fraction))
                    gfa_graph.MergeForks();
                */
                gfa_graph.AssignGroupNumber(); 
                Collection splitted = gfa_graph.SplitGroups();
                size_t graphs = splitted.size();
                rslts.splice(rslts.end(), splitted);
                
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Finished graph: " << acc << ":" << group << " " << gsize << " " << double(gsize)/initial_gsize << " " << graphs << " in " << timer.Elapsed();
            }
        }
        

    };


}; // namespace
#endif /* _GFA_ */

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

#ifndef _GuidedGraph_
#define _GuidedGraph_

using namespace std;
namespace DeBruijn {

class CGuidedGraph {
public:
    struct SeqSegment;
    typedef list<SeqSegment>::iterator GrIterator;
    struct GrIteratorHash { size_t operator()(GrIterator it) const { return std::hash<void*>()(&(*it)); } };
    typedef tuple<GrIterator, int> TSegmentP; // segment, position 
    struct SegmentPHash { size_t operator()(const TSegmentP& sp) const { return GrIteratorHash()(get<0>(sp))^std::hash<int>()(get<1>(sp)); } };
    typedef forward_list<TSegmentP> TSegmPList;

    struct SeqSegment {
        SeqSegment(const vector<SPathBase>& seq, int not_aligned_left, int not_aligned_right) : m_seq(seq), m_not_aligned_left(not_aligned_left), m_not_aligned_right(not_aligned_right) {}
        SeqSegment(const SegSeq& seq, int not_aligned_left, int not_aligned_right) : m_seq(seq), m_not_aligned_left(not_aligned_left), m_not_aligned_right(not_aligned_right) {}

        SegSeq m_seq;
        map<int, TSegmPList> m_left_connections;
        map<int, TSegmPList> m_right_connections;
        int m_not_aligned_left = 0;      
        int m_not_aligned_right = 0;      
    };

    typedef pair<Node,int> TAnchor;   // kmer, right end of kmer on target
    struct AnchorHash { size_t operator()(const TAnchor& anc) const { return Node::Hash()(get<0>(anc))^std::hash<int>()(get<1>(anc)); } };  

    TSegmentP End() { return TSegmentP(m_seq_container.end(), 0); }
    
    CGuidedGraph(int kmer_len, int secondary_kmer_len) : m_kmer_len(kmer_len), m_secondary_kmer_len(secondary_kmer_len) {}

    void StartNewAssembly(const string& init_kmer, int not_aligned_left, int not_aligned_right) {
        m_last_segments.clear();
        m_last_chain.clear();
        m_notaligned.clear();
        m_notaligned_reverse_map.clear();
        vector<SPathBase> seg;
        for(char c : init_kmer)
            seg.emplace_back(c);
        AddSegment(seg, -m_kmer_len, not_aligned_left, not_aligned_right);
    }

    void RemoveNotAlignedSegments(double anchor_frac) {
        m_last_segments.clear();
        m_last_chain.clear();
        m_notaligned.clear();
        m_notaligned_reverse_map.clear();

        stack<list<SeqSegment>::iterator> not_aligned_rights;
        stack<list<SeqSegment>::iterator> not_aligned_lefts;
        for(auto iseg = m_seq_container.begin(); iseg != m_seq_container.end(); ++iseg) {
            if(iseg->m_not_aligned_right > (int)iseg->m_seq.size()+anchor_frac*m_kmer_len && iseg->m_right_connections.empty())
                not_aligned_rights.push(iseg);
            if(iseg->m_not_aligned_left > (int)iseg->m_seq.size()+anchor_frac*m_kmer_len && iseg->m_left_connections.empty())
                not_aligned_lefts.push(iseg);
        }

        while(!not_aligned_rights.empty()) {
            auto iseg = not_aligned_rights.top();
            not_aligned_rights.pop();

            for(auto& ipos_lst : iseg->m_left_connections) {
                int ipos = ipos_lst.first;
                for(auto jp : ipos_lst.second) {
                    auto jptr = get<0>(jp);
                    int jpos = get<1>(jp);
                    auto& lst = jptr->m_right_connections[jpos];
                    lst.remove(TSegmentP(iseg, ipos));
                    if(lst.empty())
                        jptr->m_right_connections.erase(jpos);                       
                    if(jptr->m_not_aligned_right > (int)jptr->m_seq.size()+anchor_frac*m_kmer_len && jptr->m_right_connections.empty())
                        not_aligned_rights.push(jptr);                    
                }
            }

            m_seq_container.erase(iseg);
        }

        while(!not_aligned_lefts.empty()) {
            auto iseg = not_aligned_lefts.top();
            not_aligned_lefts.pop();

            for(auto& ipos_lst : iseg->m_right_connections) {
                int ipos = ipos_lst.first;
                for(auto jp : ipos_lst.second) {
                    auto jptr = get<0>(jp);
                    int jpos = get<1>(jp);
                    auto& lst = jptr->m_left_connections[jpos];
                    lst.remove(TSegmentP(iseg, ipos));
                    if(lst.empty())
                        jptr->m_left_connections.erase(jpos);                       
                    if(jptr->m_not_aligned_left > (int)jptr->m_seq.size()+anchor_frac*m_kmer_len && jptr->m_left_connections.empty())
                        not_aligned_lefts.push(jptr);                    
                }
            }

            m_seq_container.erase(iseg);
        }
    }
    void RewindLeftBranch(int starting_shift, double anchor_frac) {
        while(m_last_segments.size() > 1 && get<2>(*(m_last_segments.end()-2)) >= starting_shift) {
            GrIterator iseg = get<0>(m_last_segments.back());
            m_last_chain.erase(iseg);
            m_last_segments.pop_back();

            if(iseg->m_not_aligned_left > (int)iseg->m_seq.size()+anchor_frac*m_kmer_len && iseg->m_left_connections.empty()) {
                for(auto& ipos_lst : iseg->m_right_connections) {
                    int ipos = ipos_lst.first;
                    for(auto jp : ipos_lst.second) {
                        auto jptr = get<0>(jp);
                        int jpos = get<1>(jp);
                        auto& lst = jptr->m_left_connections[jpos];
                        lst.remove(TSegmentP(iseg, ipos));
                        if(lst.empty())
                            jptr->m_left_connections.erase(jpos);                       
                    }
                }

                for(auto& rinfo : m_notaligned_reverse_map[iseg]) {
                    auto mapi = get<0>(rinfo);
                    auto& tp = get<1>(rinfo);
                    mapi->second.remove(tp);
                    if(mapi->second.empty())
                        m_notaligned.erase(mapi);
                }
                m_notaligned_reverse_map.erase(iseg);

                m_seq_container.erase(iseg);
            }
        }
        if(m_last_segments.size() > 1) {
            int prev_len = get<2>(*(m_last_segments.end()-2));
            int len_in_seg = starting_shift-prev_len;                    // > 0
            get<1>(m_last_segments.back()) = len_in_seg;
            get<2>(m_last_segments.back()) = starting_shift;
            auto segp = get<0>(m_last_segments.back());
            m_last_chain[segp] = make_tuple(len_in_seg, starting_shift);
        }
    }
    void RewindRightBranch(int starting_shift, double anchor_frac) {
        while(m_last_segments.size() > 1 && get<2>(*(m_last_segments.end()-2)) >= starting_shift) {
            GrIterator iseg = get<0>(m_last_segments.back());
            m_last_chain.erase(iseg);
            m_last_segments.pop_back();

            if(iseg->m_not_aligned_right > (int)iseg->m_seq.size()+anchor_frac*m_kmer_len && iseg->m_right_connections.empty()) {
                for(auto& ipos_lst : iseg->m_left_connections) {
                    int ipos = ipos_lst.first;
                    for(auto jp : ipos_lst.second) {
                        auto jptr = get<0>(jp);
                        int jpos = get<1>(jp);
                        auto& lst = jptr->m_right_connections[jpos];
                        lst.remove(TSegmentP(iseg, ipos));
                        if(lst.empty())
                            jptr->m_right_connections.erase(jpos);                       
                    }
                }

                for(auto& rinfo : m_notaligned_reverse_map[iseg]) {
                    auto mapi = get<0>(rinfo);
                    auto& tp = get<1>(rinfo);
                    mapi->second.remove(tp);
                    if(mapi->second.empty())
                        m_notaligned.erase(mapi);
                }
                m_notaligned_reverse_map.erase(iseg);

                m_seq_container.erase(iseg);
            }
        }
        if(m_last_segments.size() > 1) {
            int prev_len = get<2>(*(m_last_segments.end()-2));
            int len_in_seg = starting_shift-prev_len;                    // > 0
            get<1>(m_last_segments.back()) = len_in_seg;
            get<2>(m_last_segments.back()) = starting_shift;
            auto segp = get<0>(m_last_segments.back());
            m_last_chain[segp] = make_tuple(len_in_seg, starting_shift);
        }
    }
    void CleanBranch() {
        m_last_segments.resize(1);
        m_last_chain.clear();
    }
    TSegmentP KnownRightAnchor(const TAnchor& anchor) { 
        auto rslt = m_right_anchors.find(anchor);

        if(rslt != m_right_anchors.end()) {
            TSegmentP& ancp = rslt->second;
            auto it = m_last_chain.find(get<0>(ancp));
            if(it == m_last_chain.end() || get<0>(it->second)-1 < get<1>(ancp))  // check for loop
                return ancp;
        }
        
        return End();
    }
    TSegmentP KnownLeftAnchor(const TAnchor& anchor) {
        auto rslt = m_left_anchors.find(anchor);

        if(rslt != m_left_anchors.end()) {
            TSegmentP& ancp = rslt->second;
            int plen = get<0>(ancp)->m_seq.size();
            auto it = m_last_chain.find(get<0>(ancp));
            if(it == m_last_chain.end() || get<0>(it->second) < plen-get<1>(ancp))  // check for loop 
                return ancp;
        }

        return End(); 
    }
    TSegmentP KnownRightNotAligned(Node node, int chunk_len, int starting_shift, const vector<SPathBaseKmer>& assembled_seq) {
        if(chunk_len == (int)assembled_seq.size())      // initial seq
            return End();
        auto rslt = m_notaligned.find(node);
        if(rslt == m_notaligned.end())
            return End();

        double fr = 0.2;
        for(TSegmentP& hook : rslt->second) {
            GrIterator segp = get<0>(hook);
            auto it = m_last_chain.find(segp);
            if(it != m_last_chain.end()) { 
                int fork_pos = get<0>(it->second)-1;
                int match_pos = get<1>(hook);
                int template_len = match_pos-fork_pos;
                int chain_len = starting_shift+chunk_len-get<1>(it->second);
                int diff = abs(template_len-chain_len);
                for(int p = m_kmer_len+1; p <= min(template_len, chain_len) && diff <= fr*chain_len; ++p) {
                    if(assembled_seq[assembled_seq.size()-p].m_nt != segp->m_seq[match_pos+1-p].m_nt)
                        ++diff;
                }
                if(diff <= fr*chain_len)
                    return hook; 
                else
                    return End();  // check only last matching branch
            }
        }               
        
        return End(); 
    }
    TSegmentP KnownLeftNotAligned(Node node, int chunk_len, int starting_shift, const vector<SPathBaseKmer>& assembled_seq) {
        if(chunk_len == (int)assembled_seq.size())      // initial seq
            return End();
        auto rslt = m_notaligned.find(node);
        if(rslt == m_notaligned.end())
            return End();

        double fr = 0.2;
        for(TSegmentP& hook : rslt->second) {
            GrIterator segp = get<0>(hook);
            auto it = m_last_chain.find(segp);
            if(it != m_last_chain.end()) { 
                int fork_pos = segp->m_seq.size()-get<0>(it->second);
                int match_pos = get<1>(hook);
                int template_len = fork_pos-match_pos+1;
                int chain_len = starting_shift+chunk_len-get<1>(it->second);
                int diff = abs(template_len-chain_len);
                for(int p = m_kmer_len+1; p <= min(template_len, chain_len) && diff <= fr*chain_len; ++p) {
                    if(Complement(assembled_seq[assembled_seq.size()-p].m_nt) != segp->m_seq[match_pos+p-1].m_nt)
                        ++diff;
                }
                if(diff <= fr*chain_len)
                    return hook;                
                else
                    return End();  // check only last matching branch
            }
        }               

        return End(); 
    }

    void AddLeftSegment(SPathChunk& chunk, const map<TAnchor, int>& lanchors, const map<Node, int>& lnotaligned, const TSegmentP& hook) {
        vector<SPathBase>& segment = chunk.m_seq;
        if(segment.empty())
            return;
        reverse(segment.begin(), segment.end());
        for(auto& base : segment) {
            base.m_nt = Complement(base.m_nt);
        }

        int starting_shift = chunk.m_starting_shift;
        int not_aligned = chunk.m_not_aligned;
        int shift = get<1>(m_last_segments.back()); // > 0; plen-shift is last position included from prev
        GrIterator prev_segmentp = get<0>(m_last_segments.back());
        auto& prev_segment = *prev_segmentp;

        GrIterator new_segmentp = AddSegment(segment, starting_shift, not_aligned, 0);
        SeqSegment& new_segment = *new_segmentp;
        int plen = prev_segment.m_seq.size();
        int nlen = new_segment.m_seq.size();
        prev_segment.m_left_connections[plen-shift].push_front(TSegmentP(new_segmentp, nlen-1));     // plen-1-shift+1   
        new_segment.m_right_connections[nlen-1].push_front(TSegmentP(prev_segmentp, plen-shift));
        
        // capture fork info
        TSegmentP lfork = StepRight(TSegmentP(new_segmentp, 0));
        TSegmentP rfork(new_segmentp, 0);
        TSegmentP rfork_secondary(new_segmentp, 0);
        for(int i = 0; i < m_kmer_len-1; ++i)
            rfork = StepRight(rfork);
        for(int i = 0; i < m_kmer_len-1; ++i)
            rfork_secondary = StepRight(rfork_secondary);
        for(int p = 0; p < nlen; ++p) {
            if(segment[p].m_fork & eRightFork)
                get<0>(lfork)->m_seq[get<1>(lfork)].m_fork |= eLeftFork;
            if(segment[p].m_fork & eLeftFork) {
                if(segment[p].m_fork & eSecondaryKmer)
                    get<0>(rfork_secondary)->m_seq[get<1>(rfork_secondary)].m_fork |= eRightFork;
                else
                    get<0>(rfork)->m_seq[get<1>(rfork)].m_fork |= eRightFork;
            }

            lfork = StepRight(lfork);
            rfork = StepRight(rfork);
            rfork_secondary = StepRight(rfork_secondary);
        }        

        if(get<0>(hook) != m_seq_container.end()) {
            GrIterator next_segmentp = get<0>(hook);
            int kmer_end = get<1>(hook);
                               
            new_segment.m_left_connections[0] = next_segmentp->m_left_connections[kmer_end];                          // inherit existing connections from next_segment (create empty if not existed) 
            if(kmer_end > 0)                                                                                          // anchor is not the end of the next segment
                new_segment.m_left_connections[0].push_front(TSegmentP(next_segmentp, kmer_end-1));                   // make connection to the rest of the next_segment
            for(auto& lconnection : new_segment.m_left_connections[0])
                get<0>(lconnection)->m_right_connections[get<1>(lconnection)].push_front(TSegmentP(new_segmentp, 0)); // make left connection in all newly connected segments (incuding next_segment)
        }

        for(auto& la : lanchors)
            m_left_anchors.emplace(la.first, TSegmentP(new_segmentp, nlen-1-la.second));            
        for(auto& lnal : lnotaligned) {
            auto itr = m_notaligned.emplace(lnal.first, TSegmPList()).first;
            itr->second.emplace_front(new_segmentp, nlen-1-lnal.second);
            m_notaligned_reverse_map[new_segmentp].emplace_front(itr, itr->second.front());
        }
    }
    void AddRightSegment(SPathChunk& chunk, const map<TAnchor, int>& ranchors, const map<Node, int>& rnotaligned, const TSegmentP& hook) {
        vector<SPathBase>& segment = chunk.m_seq;
        if(segment.empty())
            return;

        int starting_shift = chunk.m_starting_shift;
        int not_aligned = chunk.m_not_aligned;
        int shift = get<1>(m_last_segments.back()); // > 0; shift-1 is last position included from prev
        GrIterator prev_segmentp = get<0>(m_last_segments.back());
        auto& prev_segment = *prev_segmentp;

        GrIterator new_segmentp = AddSegment(segment, starting_shift, 0, not_aligned);
        SeqSegment& new_segment = *new_segmentp;
        prev_segment.m_right_connections[shift-1].push_front(TSegmentP(new_segmentp, 0)); // shift-1 position connects to new segment start
        new_segment.m_left_connections[0].push_front(TSegmentP(prev_segmentp, shift-1));  // new segment start connects to shift-1 position
        int nlen = new_segment.m_seq.size();

        // capture fork info
        TSegmentP rfork = StepLeft(TSegmentP(new_segmentp, nlen-1));
        TSegmentP lfork(new_segmentp, nlen-1);
        TSegmentP lfork_secondary(new_segmentp, nlen-1);
        for(int i = 0; i < m_kmer_len-1; ++i)
            lfork = StepLeft(lfork);
        for(int i = 0; i < m_secondary_kmer_len-1; ++i)
            lfork_secondary = StepLeft(lfork_secondary);
        for(int p = nlen-1; p >= 0; --p) {
            if(segment[p].m_fork & eSecondaryKmer)
                get<0>(lfork_secondary)->m_seq[get<1>(lfork_secondary)].m_fork |= (segment[p].m_fork & eLeftFork);
            else
                get<0>(lfork)->m_seq[get<1>(lfork)].m_fork |= (segment[p].m_fork & eLeftFork);
            get<0>(rfork)->m_seq[get<1>(rfork)].m_fork |= (segment[p].m_fork & eRightFork);
            lfork = StepLeft(lfork);
            lfork_secondary = StepLeft(lfork_secondary);
            rfork = StepLeft(rfork);
        }

        if(get<0>(hook) != m_seq_container.end()) {
            GrIterator next_segmentp = get<0>(hook);
            int kmer_end = get<1>(hook);
            new_segment.m_right_connections[nlen-1] = next_segmentp->m_right_connections[kmer_end];                         // inherit existing connections from next_segment (create empty if not existed)             
            if(kmer_end < (int)next_segmentp->m_seq.size()-1)                                                               // anchor is not the end of the next segment    
                new_segment.m_right_connections[nlen-1].push_front(TSegmentP(next_segmentp, kmer_end+1));                   // make connection to the rest of the next_segment
            for(auto& rconnection: new_segment.m_right_connections[nlen-1])
                get<0>(rconnection)->m_left_connections[get<1>(rconnection)].push_front(TSegmentP(new_segmentp, nlen-1));  // make left connection in all newly connected segments (incuding next_segment)
        }

        for(auto& ra : ranchors)
            m_right_anchors.emplace(ra.first, TSegmentP(new_segmentp, ra.second));  // don't update if anchor exists          
        for(auto& rnal : rnotaligned) {
            auto itr = m_notaligned.emplace(rnal.first, TSegmPList()).first;
            itr->second.emplace_front(new_segmentp, rnal.second);
            m_notaligned_reverse_map[new_segmentp].emplace_front(itr, itr->second.front());
        }
    }                        
    void GetGFAGraph(GFAGraph& gfa_segments) {
        list<SeqSegment> segments;
        unordered_map<TSegmentP, TSegmentP, SegmentPHash> link_map; 

        //break original segments
        for(GrIterator orig_segp = m_seq_container.begin(); orig_segp != m_seq_container.end(); ++orig_segp) {
            auto& orig_seg = *orig_segp;
            set<int> break_points; // break between i and i+1   
            for(auto& lc : orig_seg.m_left_connections) {
                if(lc.first > 0)
                    break_points.insert(lc.first-1);
            }
            for(auto& rc : orig_seg.m_right_connections) {
                if(rc.first < (int)orig_seg.m_seq.size()-1)
                    break_points.insert(rc.first);   
            }
            list<pair<int,int>> ranges;
            int left = 0;
            for(int right : break_points) {
                ranges.emplace_back(left, right);
                left = right+1;
            }
            ranges.emplace_back(left, orig_seg.m_seq.size()-1); 

            for(auto& range : ranges) {
                int left = range.first;
                int right = range.second;
                int len = right-left+1;
                int not_aligned_left = 0;
                if(left == 0)
                    not_aligned_left = orig_seg.m_not_aligned_left;
                int not_aligned_right = 0;
                if(right == (int)orig_seg.m_seq.size()-1)
                    not_aligned_right = orig_seg.m_not_aligned_right;
                segments.push_back(SeqSegment(orig_seg.m_seq.substr(left, len), not_aligned_left, not_aligned_right));
                SeqSegment& seg = segments.back();
                GrIterator segp = prev(segments.end());

                // at this point links are to 'old' segments        
                if(orig_seg.m_left_connections.count(left))
                    seg.m_left_connections[0] = orig_seg.m_left_connections[left];
                if(orig_seg.m_right_connections.count(right))
                    seg.m_right_connections[len-1] = orig_seg.m_right_connections[right];                    
                    
                if(left > 0) { // add new links for segment break       
                    SeqSegment& prev_seg = *next(segments.rbegin()); // garanteed to exist      
                    prev_seg.m_right_connections[prev_seg.m_seq.size()-1].emplace_front(orig_segp, left);
                    seg.m_left_connections[0].emplace_front(orig_segp, left-1);
                }                    
                    
                //populate link map     
                link_map[TSegmentP(orig_segp, left)] = TSegmentP(segp, 0);
                link_map[TSegmentP(orig_segp, right)] = TSegmentP(segp, seg.m_seq.size()-1);
            }                                    
        }

        // create gfa segments  
        gfa_segments.resize(segments.size());
        unordered_map<GrIterator, typename list<GFASegment>::iterator, GrIteratorHash> seg_to_gfa;
        auto it = gfa_segments.begin();
        for(GrIterator segp = segments.begin(); segp != segments.end(); ++segp)
            seg_to_gfa[segp] = it++; 

        it = gfa_segments.begin();
        for(SeqSegment& seg : segments) {
            GFASegment& gfa_seg = *it++;;
            swap(seg.m_seq, gfa_seg.m_seq);
            // transfer not aligned 
            gfa_seg.m_left_len = seg.m_not_aligned_left;
            gfa_seg.m_right_len = seg.m_not_aligned_right;
            // remap links      
            if(!seg.m_left_connections.empty()) {
                for(auto& lc : seg.m_left_connections.begin()->second)
                    gfa_seg.m_left_connections.push_front(seg_to_gfa[get<0>(link_map[lc])]);
            }
            if(!seg.m_right_connections.empty()) {
                for(auto& rc : seg.m_right_connections.begin()->second)
                    gfa_seg.m_right_connections.push_front(seg_to_gfa[get<0>(link_map[rc])]);
            }
        }

        // propagate not aligned to max depth
        bool keep_doing = true;
        while(keep_doing) {
            keep_doing = false;
            for(auto& seg : gfa_segments) {
                int len = seg.m_seq.size();
                if(seg.m_left_len > len) {
                    for(auto& rc : seg.m_right_connections) {
                        rc->m_left_len = seg.m_left_len-len;
                        if(rc->m_left_len > (int)rc->m_seq.size())
                            keep_doing = true;
                    }
                    seg.m_left_len = len;
                }
                if(seg.m_right_len > len) {
                    for(auto& lc : seg.m_left_connections) {
                        lc->m_right_len = seg.m_right_len - len;
                        if(lc->m_right_len > (int)lc->m_seq.size())
                            keep_doing = true;
                    }
                    seg.m_right_len = len;
                }
            }
        }
        // remove not aligned if there are 'aligned' connections
        keep_doing = true;
        while(keep_doing) {
            keep_doing = false;
            for(auto& seg : gfa_segments) {
                if(seg.m_left_len > 0) {
                    for(auto& lc : seg.m_left_connections) {
                        if(lc->m_left_len < (int)lc->m_seq.size()) {
                            seg.m_left_len = 0;
                            keep_doing = true;
                            break;
                        }
                    }
                }
                if(seg.m_right_len > 0) {
                    for(auto& rc : seg.m_right_connections) {
                        if(rc->m_right_len < (int)rc->m_seq.size()) {
                            seg.m_right_len = 0;
                            keep_doing = true;
                            break;
                        }                        
                    }
                }
            }
        }

        // clip not aligned
        for(auto it_loop = gfa_segments.begin(); it_loop != gfa_segments.end(); ) {
            auto it = it_loop++;
            auto& seg = *it;
            int len = seg.m_seq.size();

            if(seg.m_left_len >= len || seg.m_right_len >= len) {
                gfa_segments.RemoveLinksToSegment(it);
                gfa_segments.erase(it);
            } else if(seg.m_left_len > 0) {
                seg.m_left_connections.clear();
                seg.m_seq = seg.m_seq.substr(seg.m_left_len);
            } else if(seg.m_right_len > 0) {
                seg.m_right_connections.clear();
                seg.m_seq = seg.m_seq.substr(0, len-seg.m_right_len);
            }
        }

        // reset chain lengths
        for(auto& seg : gfa_segments) {
            seg.m_left_len = 0;
            seg.m_right_len = 0;
        }

        gfa_segments.Size() = gfa_segments.size();
    }

    size_t Total() const {
        size_t total = 0;
        for(auto& segm : m_seq_container) 
            total += segm.m_seq.size();
        return total;
    }
    size_t Size() const { return m_seq_container.size(); }

private:
    GrIterator AddSegment(const vector<SPathBase>& segment, int starting_shift, int not_aligned_left, int not_aligned_right) {
        m_seq_container.push_front(SeqSegment(segment, not_aligned_left, not_aligned_right));
        GrIterator new_segmentp = m_seq_container.begin();
        int seg_len = new_segmentp->m_seq.size();
        m_last_segments.emplace_back(new_segmentp, seg_len, starting_shift+seg_len);
        return new_segmentp;
    }

    // step functions follow simple tree; they DON'T make a thorough check
    TSegmentP StepRight(const TSegmentP& pos) {
        GrIterator segp = get<0>(pos);
        int p = get<1>(pos);
        if(++p < (int)segp->m_seq.size()) {
            return TSegmentP(segp, p);
        } else {
            auto last = segp->m_right_connections.rbegin();
            if(last == segp->m_right_connections.rend())
                return End();
            else
                return last->second.front();
        }
    }
    TSegmentP StepLeft(const TSegmentP& pos) {
        GrIterator segp = get<0>(pos);
        int p = get<1>(pos);
        if(--p >= 0) {
            return TSegmentP(segp, p);
        } else {
            auto first = segp->m_left_connections.begin();
            if(first == segp->m_left_connections.end())
                return End();
            else
                return first->second.front();
        }
    }

    unordered_map<TAnchor, TSegmentP, AnchorHash> m_left_anchors;
    unordered_map<TAnchor, TSegmentP, AnchorHash> m_right_anchors;
    typedef unordered_map<Node, TSegmPList, typename Node::Hash> TNotAlignedMap;
    TNotAlignedMap m_notaligned;
    unordered_map<GrIterator, forward_list<tuple<TNotAlignedMap::iterator, TSegmentP>>, GrIteratorHash> m_notaligned_reverse_map; 
    list<SeqSegment> m_seq_container;
    deque<tuple<GrIterator, int, int>> m_last_segments;                        // number bases included from segment, total number of bases     
    unordered_map<GrIterator, tuple<int,int>, GrIteratorHash> m_last_chain;  // number bases included from segment, total number of bases
    int m_kmer_len;
    int m_secondary_kmer_len;
};
} // namespace
#endif /* _GuidedGraph_ */

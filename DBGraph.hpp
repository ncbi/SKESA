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

#ifndef _DeBruijn_Graph_
#define _DeBruijn_Graph_

#include <iostream>
#include <bitset>
#include "counter.hpp"
#include "concurrenthash.hpp"

// This file contains classes which facilitate basic operation of storing reads, counting kmers,
// and creating and traversing a de Bruijn graph

using namespace std;
namespace DeBruijn {

    // Implementation of de Bruijn graph based on TKmerCount which stores kmer (smaller in the bit encoding of self and its reverse
    // complement), its count, fraction of times the stored kmer was seen as self, and information for presence/absence in graph
    // for each of the eight possible extensions to which this kmer can be connected
    // Allows basic traversing operations such as find kmer and its abundance (count) or find successors for a kmer
    // We use a node-centric definition of de Bruijn graph in which nodes of the graph are kmers
    class CDBGraph {
    public:

        // Construct graph from counted kmers and histogram
        // is_stranded indicates if count include reliable direction information (PlusFraction() and MinusFraction() could be used)
        CDBGraph(const TKmerCount& kmers, const TBins& bins, bool is_stranded) : m_graph_kmers(kmers.KmerLen()), m_bins(bins), m_is_stranded(is_stranded) {
            m_graph_kmers.PushBackElementsFrom(kmers);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
            m_visited.resize(GraphSize(), 0);
        }

        // Construct graph from temporary containers
        CDBGraph(TKmerCount&& kmers, TBins&& bins, bool is_stranded) :  m_graph_kmers(kmers.KmerLen()), m_is_stranded(is_stranded) {
            m_graph_kmers.Swap(kmers);
            m_bins.swap(bins);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
            m_visited.resize(GraphSize(), 0);
        }

        // Load from a file
        CDBGraph(istream& in) {
            string tag;
            if(!getline(in, tag) || tag != "Sorted Graph")
                throw runtime_error("Wrong format of graph file");           
            m_graph_kmers.Load(in);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);

            int bin_num;
            if(!in.read(reinterpret_cast<char*>(&bin_num), sizeof bin_num))
                throw runtime_error("Error in CDBGraph read");
            for(int i = 0; i < bin_num; ++i) {
                pair<int, size_t> bin;
                if(!in.read(reinterpret_cast<char*>(&bin), sizeof bin))
                    throw runtime_error("Error in CDBGraph read");
                m_bins.push_back(bin);
            }

            if(!in.read(reinterpret_cast<char*>(&m_is_stranded), sizeof m_is_stranded))
                throw runtime_error("Error in CDBGraph read");
            m_visited.resize(GraphSize(), 0);
        }

        // Save in a file
        void Save(ostream& out) const {
            out << "Sorted Graph\n";
            m_graph_kmers.Save(out);
            int bin_num = m_bins.size();
            out.write(reinterpret_cast<const char*>(&bin_num), sizeof bin_num);         
            out.write(reinterpret_cast<const char*>(&m_bins[0]), bin_num*(sizeof m_bins[0])); 
            out.write(reinterpret_cast<const char*>(&m_is_stranded), sizeof m_is_stranded);
            if(!out)
                throw runtime_error("Error in CDBGraph write"); 
        }

        class Node {
        public:
            explicit Node(size_t node = 0) : m_node(node) {}
            bool isValid() const { return m_node > 0; }
            bool isPlus() const { return (m_node%2 == 0); }
            bool isMinus() const { return (m_node%2 != 0); }
            Node ReverseComplement() const {
                if(m_node != 0)
                    return Node(m_node%2 == 0 ? m_node+1 : m_node-1);
                else
                    return Node(0);
            }
            Node DropStrand() const { return Node(2*(m_node/2)); }
            bool operator == (const Node& other) const { return m_node == other.m_node; }
            bool operator != (const Node& other) const { return m_node != other.m_node; }
            bool operator < (const Node& other) const { return m_node < other.m_node; }
            bool operator > (const Node& other) const { return m_node > other.m_node; }
            struct Hash { size_t operator()(const Node& node) const { return std::hash<u_int64_t>()(node.m_node); } };
        private:
            friend class CDBGraph;
            size_t Index() const { return m_node/2-1; }
            size_t m_node;
        };
        class Iterator : public Node {
        public:
            Iterator& operator++() {
                m_node += 2;
                return *this;
            }
        private:
            friend class CDBGraph;
            explicit Iterator(size_t node) : Node(node) {}
        };
        Iterator Begin() const { return Iterator(GraphSize() > 0 ? 2 : 0); }
        Iterator End() const { return Iterator(GraphSize() > 0 ? 2*(GraphSize()+1) : 0); }
        vector<Iterator> Chunks(int desired_num) {
            vector<Iterator> chunks;
            size_t step = GraphSize()/desired_num+1;
            for(size_t index = 0; index < GraphSize(); ++index) {
                if(index%step == 0)
                    chunks.push_back(Iterator(2*(index+1)));
            }
            if(!chunks.empty())
                chunks.push_back(End());

            return chunks;
        }


        // These two functions map kmers to integer indexes which could be used to retrieve kmer properties
        // 0 is returned for kmers not present in the graph
        // positive even numbers are for stored kmers
        // positive odd numbers are for reverse complement of stored kmers 
        Node GetNode(const TKmer& kmer) const {   // finds kmer in graph
            TKmer rkmer = revcomp(kmer, KmerLen());
            if(kmer < rkmer) {
                size_t index = m_graph_kmers.Find(kmer);
                return Node(index == GraphSize() ? 0 : 2*(index+1));
            } else {
                size_t index = m_graph_kmers.Find(rkmer);
                return Node(index == GraphSize() ? 0 : 2*(index+1)+1);
            }
        }
        Node GetNode(const string& kmer_seq) const {   // finds kmer in graph
            if(kmer_seq.find_first_not_of("ACGT") != string::npos || (int)kmer_seq.size() != KmerLen())   // invalid kmer
                return Node(0);
            TKmer kmer(kmer_seq);
            return GetNode(kmer);
        }

        // for all access with Node there is NO check that node is in range !!!!!!!!
        int Abundance(const Node& node) const { // total count for a kmer
            if(!node.isValid())
                return 0;
            else
                return m_graph_kmers.GetKmerCount(node.Index()).second;  // automatically clips out branching information!
        }
        // 32 bit count; 8 bit branching; 8 bit not used yet; 16 bit +/-
        double MinusFraction(const Node& node) const {  // fraction of the times kmer was seen in - direction
            double plusf = PlusFraction(node);
            return min(plusf,1-plusf);
        }
        double PlusFraction(const Node& node) const {  // fraction of the times kmer was seen in + direction
            double plusf = double(m_graph_kmers.GetKmerCount(node.Index()).second >> 48)/numeric_limits<uint16_t>::max();
            if(node.isMinus())
                plusf = 1-plusf;
            return plusf;
        }
        TKmer GetNodeKmer(const Node& node) const {  // returns kmer as TKmer
            if(node.isPlus()) 
                return m_graph_kmers.GetKmerCount(node.Index()).first;
            else
                return revcomp(m_graph_kmers.GetKmerCount(node.Index()).first, KmerLen());
        }
        string GetNodeSeq(const Node& node) const { // returnd kmer as string
            return GetNodeKmer(node).toString(KmerLen());
        }
        const uint64_t* getPointer(const Node& node) { return m_graph_kmers.getPointer(node.Index()); }        

        // multithread safe way to set visited value; returns true if value was as expected before and has been successfully changed
        // 1 is used for permanent holding; 2 is used for temporary holding; 3 for multi contig
        bool SetVisited(const Node& node, uint8_t value=1, uint8_t expected=0) {
            return m_visited[node.Index()].Set(value, expected);
        }
        void SetTempHolding(const Node& node) { SetVisited(node, 2, 1); }
        void SetMultContig(const Node& node) { SetVisited(node, 3, 1); }

        bool ClearVisited(const Node& node) { // multithread safe way to clear visited value; returns true if value was set before
            return m_visited[node.Index()].Set(0, 1) || m_visited[node.Index()].Set(0, 2) || m_visited[node.Index()].Set(0, 3);
        }

        uint8_t IsVisited(const Node& node) const { // returns visited value
            return m_visited[node.Index()];
        }
        bool IsMultContig(const Node& node) const { return IsVisited(node) == 3; }

        void ClearHoldings() { // clears temporary holdings
            for(auto& v : m_visited) 
                if(v == 2) v = 0;            
        }

        void ClearAllVisited() { // clears all visited
            for(auto& v : m_visited) 
                v = 0;            
        }
        
        struct Successor {
            Successor(const Node& node, char c) : m_node(node), m_nt(c) {}
            Node m_node;
            char m_nt;
            bool operator == (const Successor& other) const { return m_node == other.m_node; }
            bool operator != (const Successor& other) const { return m_node != other.m_node; }
            bool operator < (const Successor& other) const { return m_node < other.m_node; }
        };
        
        // Returns successors of a node 
        // These are nodes representing kmers produced by extending the right end of the kmer for
        // this node by one base and removing the leftmost base of the kmer
        // Each successor stores the successor's node and the extra base
        // Finding predecessors is done by finding successors of reverse complement of the kmer for the node
        vector<Successor> GetNodeSuccessors(const Node& node) const {
            vector<Successor> successors;
            if(!node.isValid())
                return successors;

            uint8_t branch_info = (m_graph_kmers.GetCount(node.Index()) >> 32);
            bitset<4> branches(node.isMinus() ? (branch_info >> 4) : branch_info);
            if(branches.count()) {
                TKmer shifted_kmer = (GetNodeKmer(node) << 2) & m_max_kmer;
                for(int nt = 0; nt < 4; ++nt) {
                    if(branches[nt]) {
                        Node successor = GetNode(shifted_kmer + TKmer(KmerLen(), nt));
                        successors.push_back(Successor(successor, bin2NT[nt]));
                    }
                }
            }

            return successors;            
        }

        // Revese complement node
        static Node ReverseComplement(Node node) { return node.ReverseComplement(); }

        int KmerLen() const { return m_graph_kmers.KmerLen(); }             // returns kmer length
        size_t GraphSize() const { return m_graph_kmers.Size(); }           // returns total number of elements
        size_t ElementSize() const { return m_graph_kmers.ElementSize(); }  // element size in bytes
        size_t MemoryFootprint() const {                                    // reserved memory in bytes
            return m_graph_kmers.MemoryFootprint()+m_visited.capacity()+sizeof(TBins::value_type)*m_bins.capacity(); 
        }
        bool GraphIsStranded() const { return m_is_stranded; }              // indicates if graph contains stranded information

        // returns minimum position for stored histogram
        int HistogramMinimum() const {
            pair<int,int> r = HistogramRange(m_bins);
            if(r.first < 0)
                return 0;
            else 
                return m_bins[r.first].first;
        }

        // useus simple heuristic to evaluate the genome size
        size_t GenomeSize() const { return CalculateGenomeSize(m_bins); }
        // returns histogram
        const TBins& GetBins() const { return m_bins; }
        // average count of kmers in the histogram with the main peak
        double AverageCount() const { return GetAverageCount(m_bins); }

    private:

        TKmerCount m_graph_kmers;     // only the minimal kmers are stored  
        TKmer m_max_kmer;             // contains 1 in all kmer_len bit positions  
        TBins m_bins;
        vector<SAtomic<uint8_t>> m_visited;
        bool m_is_stranded;
    };


    class CDBHashGraph {
    public:
        // Construct graph from temporary containers
        CDBHashGraph(CKmerHashCount&& kmers, bool is_stranded) :  m_graph_kmers(kmers.KmerLen()), m_is_stranded(is_stranded) {
            m_graph_kmers.Swap(kmers);
            m_bins = m_graph_kmers.GetBins();
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
            m_graph_size = 0;
            for(auto& bin : m_bins)
                m_graph_size += bin.second;
        }

        // Load from a file
        CDBHashGraph(istream& in) {
            string tag;
            if(!getline(in, tag) || tag != "Hash Graph")
                throw runtime_error("Wrong format of graph file");           
            m_graph_kmers.Load(in);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);

            int bin_num;
            if(!in.read(reinterpret_cast<char*>(&bin_num), sizeof bin_num))
                throw runtime_error("Error in CDBHashGraph read");
            for(int i = 0; i < bin_num; ++i) {
                pair<int, size_t> bin;
                if(!in.read(reinterpret_cast<char*>(&bin), sizeof bin))
                    throw runtime_error("Error in CDBHashGraph read");
                m_bins.push_back(bin);
            }
            m_graph_size = 0;
            for(auto& bin : m_bins)
                m_graph_size += bin.second;
            if(!in.read(reinterpret_cast<char*>(&m_is_stranded), sizeof m_is_stranded))
                throw runtime_error("Error in CDBHashGraph read");
            ClearAllVisited();
        }

        // Save in a file
        void Save(ostream& out) const {
            out << "Hash Graph\n";
            m_graph_kmers.Save(out);
            int bin_num = m_bins.size();
            out.write(reinterpret_cast<const char*>(&bin_num), sizeof bin_num);         
            out.write(reinterpret_cast<const char*>(&m_bins[0]), bin_num*(sizeof m_bins[0])); 
            out.write(reinterpret_cast<const char*>(&m_is_stranded), sizeof m_is_stranded);
            if(!out)
                throw runtime_error("Error in CDBHashGraph write");
        }

        class Node : public CKmerHashCount::Index {
        public:
            enum Status : int8_t { eMinus = -1, eNotValid = 0, ePlus = 1 };
            Node() : Index(), m_status(eNotValid) {}
            Node(CKmerHashCount::Index index, Status status) : Index(index), m_status(status) {}
            Node(CKmerHashCount::Iterator iter) : Index(iter), m_status(ePlus) {}
            bool isValid() const { return m_status != eNotValid; }
            bool isPlus() const { return m_status > 0; }
            bool isMinus() const { return m_status < 0; }
            Node DropStrand() const { 
                Node node = *this;
                if(isMinus())
                    node.m_status = ePlus;
                return node;
            }
            Node ReverseComplement() const {
                Node node = *this;
                switch(m_status) {
                case eMinus : node.m_status = ePlus; return node;
                case eNotValid : return node;
                case ePlus: node.m_status = eMinus; return node;
                }
                return node;
            }
            bool operator==(const Node& other) const { return Index::operator==(other) && m_status == other.m_status; }
            bool operator!=(const Node& other) const { return !operator==(other); }            
            bool operator<(const Node& other) const {
                if(Index::operator==(other))
                    return m_status < other.m_status;
                else
                    return Index::operator<(other);
            }
            bool operator>(const Node& other) const {
                if(Index::operator==(other))
                    return m_status > other.m_status;
                else
                    return Index::operator>(other);
            }            
            struct Hash { size_t operator()(const Node& node) const { return Index::Hash()(node)^std::hash<int8_t>()(node.m_status); } };
        private:
            Status m_status;
        };

        typedef CKmerHashCount::Iterator Iterator;
        Iterator Begin() { return m_graph_kmers.Begin(); }
        Iterator End() { return m_graph_kmers.End(); }
        vector<Iterator> Chunks(int desired_num) { return m_graph_kmers.Chunks(desired_num); }
 
        Node GetNode(const TKmer& kmer) const {   // finds kmer in graph
            TKmer rkmer = revcomp(kmer, KmerLen());
            CKmerHashCount::Index end = m_graph_kmers.EndIndex();
            if(kmer < rkmer) {
                CKmerHashCount::Index index = const_cast<CKmerHashCount&>(m_graph_kmers).FindIndex(kmer);
                return Node(index, index == end ? Node::eNotValid : Node::ePlus);
            } else {
                CKmerHashCount::Index index = const_cast<CKmerHashCount&>(m_graph_kmers).FindIndex(rkmer);
                return Node(index, index == end ? Node::eNotValid : Node::eMinus);
            }
        }
        Node GetNode(const string& kmer_seq) {   // finds kmer in graph
            if(kmer_seq.find_first_not_of("ACGT") != string::npos || (int)kmer_seq.size() != KmerLen())   // invalid kmer
                return Node();
            TKmer kmer(kmer_seq);
            return GetNode(kmer);
        }
        // Revese complement node
        static Node ReverseComplement(const Node& node) { return node.ReverseComplement(); }

        // for all access with Node there is NO check that node is in range !!!!!!!!
        int Abundance(const Node& node) const { // total count for a kmer
            if(!node.isValid())
                return 0;
            else
                return node.GetMapped(m_graph_kmers)->m_data;  // automatically clips out branching information!
        }
        // 32 bit count; 8 bit branching; 8 bit visited control; 16 bit +/-
        double MinusFraction(const Node& node) const {  // fraction of the times kmer was seen in - direction
            double plusf = PlusFraction(node);
            return min(plusf,1-plusf);
        }
        double PlusFraction(const Node& node) const {  // fraction of the times kmer was seen in + direction
            double plusf = double(node.GetMapped(m_graph_kmers)->m_data >> 48)/numeric_limits<uint16_t>::max();
            if(node.isMinus())
                plusf = 1-plusf;
            return plusf;
        }
        TKmer GetNodeKmer(const Node& node) const {  // returns kmer as TKmer
            if(node.isPlus()) 
                return node.GetElement(m_graph_kmers).first;
            else
                return revcomp(node.GetElement(m_graph_kmers).first, KmerLen());
        }
        string GetNodeSeq(const Node& node) const { // returnd kmer as string
            return GetNodeKmer(node).toString(KmerLen());
        }
        const uint64_t* getPointer(const Node& node) const { return node.GetKeyPointer(m_graph_kmers); }        

        enum Visited : uint64_t {eNull = 0, eVisited = 0x10000000000, eTemp = 0x20000000000, eMulti = 0x40000000000, eAll = 0xFF0000000000 };
        // multithread safe way to set visited value; returns true if value was as expected before and has been successfully changed
        // 1 is used for permanent holding; 2 is used for temporary holding; 4 for multi contig
        bool SetVisited(const Node& node, Visited value = eVisited, Visited expected = eNull) {
            // we assume that other bits are const
            auto& count = node.GetMapped(m_graph_kmers)->m_data;
            uint64_t other_bits = (~eAll)&count.Load();
            return count.Set(other_bits|value, other_bits|expected);
        }
        void SetTempHolding(const Node& node) { SetVisited(node, eTemp, eVisited); }
        void SetMultContig(const Node& node) { SetVisited(node, eMulti, eVisited); }
        void ClearVisited(const Node& node) { node.GetMapped(m_graph_kmers)->m_data.m_atomic &= ~eAll; }
        uint64_t IsVisited(const Node& node) const { return eAll&node.GetMapped(m_graph_kmers)->m_data; }
        bool IsMultContig(const Node& node) const { return eMulti&node.GetMapped(m_graph_kmers)->m_data; }
        void ClearHoldings() { // clears temporary holdings
            for(auto it = m_graph_kmers.Begin(); it != m_graph_kmers.End(); ++it) {
                auto& count = it.GetMapped()->m_data;
                if(eTemp&count)
                    count.m_atomic &= ~eAll;
            }
        }
        void ClearAllVisited() { // clears all visited
            for(auto it = m_graph_kmers.Begin(); it != m_graph_kmers.End(); ++it)
               it.GetMapped()->m_data .m_atomic &= ~eAll;
        }
        
        struct Successor {
            Successor(const Node& node, char c) : m_node(node), m_nt(c) {}
            Node m_node;
            char m_nt;
            bool operator == (const Successor& other) const { return m_node == other.m_node; }
            bool operator != (const Successor& other) const { return m_node != other.m_node; }
            bool operator < (const Successor& other) const { return m_node < other.m_node; }
        };
            
        // Returns successors of a node 
        // These are nodes representing kmers produced by extending the right end of the kmer for
        // this node by one base and removing the leftmost base of the kmer
        // Each successor stores the successor's node and the extra base
        // Finding predecessors is done by finding successors of reverse complement of the kmer for the node
        vector<Successor> GetNodeSuccessors(const Node& node) const {
            vector<Successor> successors;
            if(!node.isValid())
                return successors;

            uint8_t branch_info = node.GetMapped(m_graph_kmers)->m_data >> 32;
            bitset<4> branches(node.isMinus() ? (branch_info >> 4) : branch_info);
            if(branches.count()) {
                TKmer shifted_kmer = (GetNodeKmer(node) << 2) & m_max_kmer;
                for(int nt = 0; nt < 4; ++nt) {
                    if(branches[nt]) {
                        Node successor = GetNode(shifted_kmer + TKmer(KmerLen(), nt));
                        successors.push_back(Successor(successor, bin2NT[nt]));
                    }
                }
            }

            return successors;            
        }    
        
        bool GraphIsStranded() const { return m_is_stranded; }              // indicates if graph contains stranded information
        int KmerLen() const { return m_graph_kmers.KmerLen(); }             // returns kmer length
        // returns minimum position for stored histogram
        int HistogramMinimum() const {
            pair<int,int> r = HistogramRange(m_bins);
            if(r.first < 0)
                return 0;
            else 
                return m_bins[r.first].first;
        }

        // useus simple heuristic to evaluate the genome size
        size_t GenomeSize() const { return CalculateGenomeSize(m_bins); }
        // returns histogram
        const TBins& GetBins() const { return m_bins; }                     // returns histogram
        // average count of kmers in the histogram with the main peak
        double AverageCount() const { return GetAverageCount(m_bins); }
        size_t GraphSize() const { return m_graph_size; }

    private:
        CKmerHashCount m_graph_kmers;
        TKmer m_max_kmer;             // contains 1 in all kmer_len bit positions  
        TBins m_bins;
        size_t m_graph_size;
        bool m_is_stranded;
    };    
    

}; // namespace
#endif /* _DeBruijn_Graph_ */

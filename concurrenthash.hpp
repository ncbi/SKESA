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

#ifndef _Concurrent_Hash_
#define _Concurrent_Hash_


#include "Integer.hpp"
#include "common_util.hpp"

// This file contains classes which facilitate basic operation of storing reads, counting kmers,
// and creating and traversing a de Bruijn graph

using namespace std;
namespace DeBruijn {

    template<int BlockSize> // in bytes
    class CConcurrentBlockedBloomFilter {
    public:
        enum EInsertResult {eNewKmer = 0, eAboveThresholdKmer = 1, eExistingKmer = 2};
        // table_size - number of counting elements in bloom filter
        // counter_bit_size - number of bith per counter (2, 4, 8)
        // hash_num - number of has functions (generated from two) 
        CConcurrentBlockedBloomFilter(size_t table_size, int counter_bit_size, int hash_num, int min_count) {
            Reset(table_size, counter_bit_size, hash_num, min_count);
        } 
        void Reset(size_t table_size, int counter_bit_size, int hash_num, int min_count) {
            assert(counter_bit_size <= 8);
            m_counter_bit_size = counter_bit_size;
            m_hash_num = hash_num;           
            m_max_element = (1 << m_counter_bit_size) - 1;
            m_min_count = min(min_count, m_max_element);
            m_elements_in_block = 8*BlockSize/m_counter_bit_size;
            m_blocks = ceil((double)table_size/m_elements_in_block);
            m_table_size = m_blocks*m_elements_in_block;
            m_count_table.clear();
            m_status.clear();
            m_count_table.resize(m_blocks);
            m_status.resize(m_blocks);
        }
        EInsertResult Insert(size_t hashp, size_t hashm) {
            size_t ind = hashp%m_blocks;
            auto& block = m_count_table[ind];
            if(Test(hashp, hashm) >= m_min_count)
                return eExistingKmer;

            while(!m_status[ind].Set(1));
            int count = Test(hashp, hashm);            
            if(count >= m_min_count) {
                m_status[ind] = 0;
                return eExistingKmer;
            }
            
            for(int h = 1; h < m_hash_num; ++h) {
                hashp += hashm;
                size_t pos = (hashp&(m_elements_in_block-1))*m_counter_bit_size;  // bit position of the counting element in block
                auto& cell = block.m_data[pos >> m_bits_in_cell_log];
                int shift = pos&(m_bits_in_cell-1);
                int cnt = (cell >> shift)&m_max_element;
                if(cnt <= count)
                    cell += ((TCell)1 << shift);
            }
            m_status[ind] = 0;

            if(count == 0)
                return eNewKmer;
            else if(count == m_min_count-1)
                return eAboveThresholdKmer;
            else
                return eExistingKmer;
        }
        int Test(size_t hashp, size_t hashm) const {
            auto& block = m_count_table[hashp%m_blocks];
            int count = m_max_element;
            for(int h = 1; h < m_hash_num; ++h) {
                hashp += hashm;
                size_t pos = (hashp&(m_elements_in_block-1))*m_counter_bit_size;  // bit position of the counting element in block
                auto& cell = block.m_data[pos >> m_bits_in_cell_log];
                int shift = pos&(m_bits_in_cell-1);
                int cnt = (cell >> shift)&m_max_element;
                if(cnt < count)
                    count = cnt;
            }

            return count;
        }
        int MaxElement() const { return m_max_element; }
        int HashNum() const { return m_hash_num; }
        size_t TableSize() const { return m_table_size; } // number of counters
        size_t TableFootPrint() const { return (sizeof(SBloomBlock)+1)*m_count_table.size(); } // bytes

    private:
        typedef uint64_t TCell;
        struct alignas(64) SBloomBlock {
            SBloomBlock() { memset(m_data.data(), 0, BlockSize); }
            array<TCell, BlockSize/sizeof(TCell)> m_data;
        };
        vector<SBloomBlock> m_count_table;
        vector<SAtomic<uint8_t>> m_status;

        size_t m_elements_in_block;
        size_t m_blocks;
        size_t m_table_size;
        int m_counter_bit_size;
        int m_hash_num;
        int m_min_count;
        int m_max_element;
        int m_bits_in_cell = 8*sizeof(TCell);
        int m_bits_in_cell_log = log(m_bits_in_cell)/log(2);
    };
    typedef CConcurrentBlockedBloomFilter<128> TConcurrentBlockedBloomFilter;


    // minimalistic multithread safe forward list
    // allows reading and inserts in the beginning
    // reading thread will not see new entries inserted after reading started
    template <typename E>
    class CForwardList {
    public:
        struct SNode {
            E m_data = E();
            SNode* m_next = nullptr;
        };

        template<bool is_const>
        class iterator : public std::iterator<forward_iterator_tag, E> {
        public:
            template <bool flag, class IsTrue, class IsFalse> struct choose;
            template <class IsTrue, class IsFalse> struct choose<true, IsTrue, IsFalse> {  typedef IsTrue type; };
            template <class IsTrue, class IsFalse> struct choose<false, IsTrue, IsFalse> {  typedef IsFalse type; };
            typedef typename choose<is_const, const E&, E&>::type reference;
            typedef typename choose<is_const, const E*, E*>::type pointer;
            typedef typename choose<is_const, const SNode*, SNode*>::type node_pointer;
            iterator(SNode* node = nullptr) : m_node(node) {};
            iterator& operator++() {
                m_node = m_node->m_next;
                return *this;
            }
            reference& operator*() { return m_node->m_data; }
            pointer operator->() { return &m_node->m_data; }
            node_pointer NodePointer() { return m_node; }
            bool operator!=(const iterator& i) const { return i.m_node != m_node; }

        private:
            friend class CForwardList;
            SNode* m_node;
        };
        iterator<false> begin() { return iterator<false>(m_head.load()); }
        iterator<false> end() { return iterator<false>(); }
        iterator<true> begin() const { return iterator<true>(m_head.load()); }
        iterator<true> end() const { return iterator<true>(); }

        CForwardList() { m_head.store(nullptr); }
        // not mutithread safe
        CForwardList& operator=(const CForwardList& other) {
            Clear();
            for(auto it = other.begin(); it != other.end(); ++it)
                PushFront(*it);
            return *this;
        }
        CForwardList(const CForwardList& other) { 
            m_head.store(nullptr);
            for(auto it = other.begin(); it != other.end(); ++it)
                PushFront(*it);
        }
       ~CForwardList() { Clear(); }
        E& front() { return m_head.load()->m_data; }
        SNode* Head() const { return m_head; }
        SNode* NewNode(const E& e) {
            SNode* p = new SNode;
            p->m_data = e;
            p->m_next = m_head;
            return p;
        }
        SNode* NewNode() {
            SNode* p = new SNode;
            p->m_next = m_head;
            return p;
        }
        bool TryPushFront(SNode* nodep) { return m_head.compare_exchange_strong(nodep->m_next, nodep); }
        E* Emplace() {
            SNode* p = NewNode();
            while (!m_head.compare_exchange_strong(p->m_next, p));
            return &(p->m_data);
        }
        void PushFront(const E& e) {
            SNode* p = NewNode(e);
            while (!m_head.compare_exchange_strong(p->m_next, p));
        }
        // not mutithread safe
        template <typename P> void remove_if(const P& prob) {
            while(m_head.load() != nullptr && prob(m_head.load()->m_data)) {
                SNode* p = m_head;
                m_head = p->m_next;
                delete p;
            }
            for(SNode* p = m_head; p != nullptr; ) {
                SNode* after = p->m_next;
                if(after != nullptr && prob(after->m_data)) {
                    p->m_next = after->m_next;
                    delete after;
                } else {
                    p = p->m_next;
                }
            }
        }
        void Save(ostream& os) const {
            size_t vsize = sizeof(E);
            os.write(reinterpret_cast<const char*>(&vsize), sizeof vsize);
            size_t elements = distance(begin(), end());
            os.write(reinterpret_cast<const char*>(&elements), sizeof elements);
            for(auto& elem : *this) 
                os.write(reinterpret_cast<const char*>(&elem), sizeof elem);
            if(!os)
                throw runtime_error("Error in CForwardList write"); 
        }
        void Load(istream& is) { 
            Clear();
            size_t vsize; 
            if(!is.read(reinterpret_cast<char*>(&vsize), sizeof vsize))
                throw runtime_error("Error in CForwardList read");
            if(vsize != sizeof(E))
                throw runtime_error("Wrong format for CForwardList load");
            size_t elements;
            if(!is.read(reinterpret_cast<char*>(&elements), sizeof elements))
                throw runtime_error("Error in CForwardList read");
            for( ; elements > 0; --elements) {
                E* p = Emplace();
                if(!is.read(reinterpret_cast<char*>(p), sizeof *p))
                    throw runtime_error("Error in CForwardList read"); 
            }
        }
        void Clear() {
            for(SNode* p = m_head; p != nullptr; ) {
                auto tmp = p->m_next;
                delete p;
                p = tmp;
            }
            m_head.store(nullptr);
        }
        void Init() { m_head.store(nullptr); }
    private:
        atomic<SNode*> m_head;
    };

    // minimalistic deque-type container allowing multithread initialization of a large hash_table
    template <typename E>
    class CDeque {
    public:
        typedef E value_type;
        CDeque(size_t size = 0) : m_chunks(1) {
            m_data.resize(m_chunks);
            Reset(size, m_chunks);
        }        
        ~CDeque() {
            if(m_chunks  > 1) {
                list<function<void()>> jobs;
                for(unsigned chunk = 0; chunk < m_chunks; ++chunk)
                    jobs.push_back(bind(&CDeque::ReleaseChunk, this, chunk));                
                RunThreads(m_chunks, jobs);                            
            }
        }        
        
        E& operator[](size_t index) { return m_data[index/m_chunk_size][index%m_chunk_size]; }        
        const E& operator[](size_t index) const { return m_data[index/m_chunk_size][index%m_chunk_size]; }        
        size_t Size() const { return m_size; }
        void Reset(size_t size, size_t chunks) {
            m_chunks = chunks;
            m_data.resize(m_chunks);
            m_size = size;
            m_chunk_size = (size+m_chunks-1)/m_chunks;
            if(m_chunks == 1) {
                ResetChunk(0, m_chunk_size);
            } else {
                list<function<void()>> jobs;
                for(unsigned chunk = 0; chunk < m_chunks; ++chunk) {
                    size_t chunk_size = min(m_chunk_size, size);
                    jobs.push_back(bind(&CDeque::ResetChunk, this, chunk, chunk_size));
                    size -= min(m_chunk_size, size);                    
                }
                RunThreads(m_chunks, jobs);                            
            }
        }
        void Swap(CDeque& other) {
            swap(m_chunks, other.m_chunks);
            swap(m_size, other.m_size);
            swap(m_chunk_size, other.m_chunk_size);
            swap(m_data, other.m_data);
        } 
        void Save(ostream& os) const {
            size_t vsize = sizeof(E);
            os.write(reinterpret_cast<const char*>(&vsize), sizeof vsize);
            os.write(reinterpret_cast<const char*>(&m_chunks), sizeof m_chunks);
            os.write(reinterpret_cast<const char*>(&m_size), sizeof m_size);
            os.write(reinterpret_cast<const char*>(&m_chunk_size), sizeof m_chunk_size);
            for(auto& chunk : m_data) {
                size_t num = chunk.size();
                os.write(reinterpret_cast<const char*>(&num), sizeof num);
                os.write(reinterpret_cast<const char*>(chunk.data()), num*vsize);
            }
            if(!os)
                throw runtime_error("Error in CDeque write"); 
        }
        void Load(istream& is) {
            size_t vsize; 
            if(!is.read(reinterpret_cast<char*>(&vsize), sizeof vsize))
                throw runtime_error("Error in CDeque read");
            if(vsize != sizeof(E))
                throw runtime_error("Wrong format for CDeque load");
            if(!is.read(reinterpret_cast<char*>(&m_chunks), sizeof m_chunks))
                throw runtime_error("Error in CDeque read");
            if(!is.read(reinterpret_cast<char*>(&m_size), sizeof m_size))
                throw runtime_error("Error in CDeque read");
            if(!is.read(reinterpret_cast<char*>(&m_chunk_size), sizeof m_chunk_size))
                throw runtime_error("Error in CDeque read");
            m_data.clear();
            m_data.resize(m_chunks);
            for(auto& chunk : m_data) {
                size_t num;
                if(!is.read(reinterpret_cast<char*>(&num), sizeof num))
                    throw runtime_error("Error in CDeque read");
                chunk.resize(num);
                if(!is.read(reinterpret_cast<char*>(chunk.data()), num*vsize))
                    throw runtime_error("Error in CDeque read");
            }
        }
    private:
        void ResetChunk(size_t chunk, size_t chunk_size) {
            m_data[chunk].clear();
            m_data[chunk].resize(chunk_size);
        }        
        void ReleaseChunk(size_t chunk) { vector<E>().swap(m_data[chunk]); }
        
        size_t m_chunks = 0;
        size_t m_size = 0;
        size_t m_chunk_size = 0;
        vector<vector<E>> m_data;
    };

    // BucketBlock <= 32
    // Moderate value of BucketBlock will improve memory cache use
    // Larger values will reduce the number of entries in the spillover lists but eventually will increase the search time
    // 0 (all entries in the lists) is permitted and could be used for reduction of the table size for large sizeof(V)
    template<class Key, class V, int BucketBlock>
    struct SHashBlock {
        static_assert(BucketBlock <= 32, "");
        typedef Key large_t;
        typedef V mapped_t;
        typedef pair<large_t,V> element_t;
        typedef CForwardList<element_t> list_t;
        typedef typename list_t::SNode snode_t;
        enum States : uint64_t {eAssigned = 1, eKeyExists = 2};
        enum { eBucketBlock = BucketBlock };
        SHashBlock() : m_status(0) {}
        SHashBlock(const SHashBlock& other) : m_data(other.m_data), m_extra(other.m_extra), m_status(other.m_status.load()) {} // used for table initialisation only
        SHashBlock& operator=(const SHashBlock& other) {
            m_data = other.m_data;
            m_extra = other.m_extra;
            m_status.store(other.m_status.load());
            return *this;
        }
            
        pair<int, typename list_t::SNode*> Find(const large_t& k, int hint) {
            if(BucketBlock > 0) {
                //try exact position first
                if(isEmpty(hint)) {
                    return make_pair(BucketBlock+1, nullptr);
                } else {
                    Wait(hint);
                    if(m_data[hint].first == k)
                        return make_pair(hint, nullptr);
                }

                //scan array        
                for(int shift = 0; shift < BucketBlock; ++shift) {
                    if(shift != hint) {
                        if(isEmpty(shift)) {
                            return make_pair(BucketBlock+1, nullptr);
                        } else {
                            Wait(shift);
                            if(m_data[shift].first == k)
                                return make_pair(shift, nullptr);
                        }
                    }
                }       
            }

            //scan spillover list   
            for(auto it = m_extra.begin(); it != m_extra.end(); ++it) {
                auto& cell = *it;
                if(cell.first == k)
                    return make_pair(BucketBlock, it.NodePointer());
            }            
            return make_pair(BucketBlock+1, nullptr);            
        }

        // 1. Try to put to exact position prescribed by hash   
        // 2. Put in the lowest available array element 
        // 3. Put in the spillover list 
        mapped_t* FindOrInsert(const large_t& k, int hint) {
            auto TryCell = [&](int shift) {
                auto& cell = m_data[shift];
                //try to grab   
                if(Lock(shift, k))
                    return &cell.second;
                
                //already assigned to some kmer 
                //wait if kmer is not stored yet    
                Wait(shift);
               
                if(cell.first == k) // kmer matches 
                    return &cell.second; 
                else
                    return (mapped_t*)nullptr; // other kmer        
            };

            if(BucketBlock > 0) {
                //try exact position first  
                auto rslt = TryCell(hint);
                if(rslt != nullptr)
                    return rslt;            

                //scan remaining array      
                for(int shift = 0; shift < BucketBlock; ++shift) {
                    if(shift != hint) {
                        auto rslt = TryCell(shift);
                        if(rslt != nullptr)
                            return rslt;
                    }
                }
            }

            //scan spillover list       
            auto existing_head = m_extra.Head();
            for(auto p = existing_head; p != nullptr; p = p->m_next) {
                if(p->m_data.first == k)
                    return &(p->m_data.second); 
            }

            typename list_t::SNode* nodep = new typename list_t::SNode;
            nodep->m_data.first = k;
            nodep->m_next = existing_head;
            while(!m_extra.TryPushFront(nodep)) {
                //check if a new elemet matches     
                for(auto p = nodep->m_next; p != existing_head; p = p->m_next) {
                    if(p->m_data.first == k) {
                        delete nodep;
                        return &(p->m_data.second); 
                    }
                }
                existing_head = nodep->m_next;
            }                

            return &(nodep->m_data.second);
        }

        element_t* IndexGet(int shift, void* lstp) {
            if(shift < BucketBlock) {    // array element   
                return &m_data[shift];
            } else {                     //list element 
                snode_t* ptr = reinterpret_cast<snode_t*>(lstp);
                return &(ptr->m_data);  
            }              
        }

        bool Lock(int shift, const large_t& kmer) {
            uint64_t assigned = eAssigned << 2*shift;
            uint64_t expected = m_status;

            do { if(expected&assigned) return false; } 
            while(!m_status.compare_exchange_strong(expected, expected|assigned));

            m_data[shift].first = kmer;
            m_status |= eKeyExists << 2*shift;
            return true;
        }
        void Wait(int shift) { 
            uint64_t keyexists = eKeyExists << 2*shift;
            while(!(m_status&keyexists)); 
        }
        bool isEmpty(int shift) const { 
            uint64_t assigned = eAssigned << 2*shift;
            return (!(m_status&assigned)); 
        }

        void Move(element_t& cell, int to) {
            m_data[to] = cell;
            m_status |= (eAssigned|eKeyExists) << 2*to;
            cell.second = V();
        }
        void Move(int from, int to) {
            Move(m_data[from], to);
            m_status &= ~((eAssigned|eKeyExists) << 2*from); // clear bits  
        }
        void Clear(int shift) {
            m_data[shift].second = V();
            m_status &= ~((uint64_t)(eAssigned|eKeyExists) << 2*shift); // clear bits   
        }

        //assumes that cel is not assigned; not mutithread safe
        mapped_t* InitCell(const large_t& kmer, int shift) { 
            if(shift < BucketBlock) {
                m_status |= (eAssigned|eKeyExists) << 2*shift;
                m_data[shift].first = kmer;
                return &(m_data[shift].second);
            }

            auto nodep = m_extra.NewNode();
            nodep->m_data.first = kmer;
            return &(nodep->m_data.second);
        }      

        array<element_t, BucketBlock> m_data;
        list_t m_extra;
        atomic<uint64_t> m_status;
    };

    //TODO needs testing
    template<typename Key, class MappedV, int BucketBlock, class Hash = std::hash<Key>>
    class CHashMap {
    public:
        CHashMap(size_t size) {
            size_t blocks = size/max(1,BucketBlock);
            if(size%max(1,BucketBlock))
                ++blocks;
            m_table_size = max(1,BucketBlock)*blocks;
            m_hash_table.Reset(blocks,1);
        }
        MappedV* Find(const Key& k) { 
            if(m_table_size == 0) {
                return nullptr;
            } else {
                size_t pos = Hash()(k)%m_table_size;
                auto& bucket = m_hash_table[pos/max(1,BucketBlock)];
                int hint = pos%max(1,BucketBlock);               

                auto rslt = bucket.Find(k, hint);
                if(rslt.first < BucketBlock)                 // found in array  
                    return &bucket.m_data[rslt.first].second;
                else if(rslt.first == BucketBlock)           // found in list   
                    return &rslt.second->m_data.second;
                else                                         // not found   
                    return nullptr; 
            }
        } 
        MappedV* FindOrInsert(const Key& k) {
            size_t index = Hash()(k)%m_table_size;
            size_t bucket_num = index/max(1,BucketBlock);
            int hint = index%max(1,BucketBlock);               
            return m_hash_table[bucket_num].FindOrInsert(k, hint);
        }
        size_t TableSize() const { return m_table_size; }

    private:
        CDeque<SHashBlock<Key, MappedV, BucketBlock>> m_hash_table;
        size_t m_table_size;
    };

    template <typename MappedV, int BucketBlock> 
    class CKmerHashMap {
    public:
        CKmerHashMap(int kmer_len = 0, size_t size = 0) : m_kmer_len(kmer_len) {
            if(m_kmer_len > 0)
                m_hash_table = CreateVariant<TKmerHashTable<MappedV>, THashBlockVec, MappedV>((m_kmer_len+31)/32);
            Reset(size, 1);
        }
        void Reset(size_t size, size_t chunks) {
            size_t blocks = size/max(1,BucketBlock);
            if(size%max(1,BucketBlock))
                ++blocks;
            m_table_size = max(1,BucketBlock)*blocks;
            apply_visitor(resize(blocks, chunks), m_hash_table);
        }

        class Index {
        public:
            Index(size_t ind = 0, void* ptr = nullptr) : m_index(ind), m_lstp(ptr) {}
            void Advance(CKmerHashMap& hash) { apply_visitor(CKmerHashMap::index_advance(*this), hash.m_hash_table); }
            pair<TKmer, MappedV*> GetElement(CKmerHashMap& hash) const { return apply_visitor(CKmerHashMap::index_get(*this), hash.m_hash_table); };
            pair<TKmer, const MappedV*> GetElement(const CKmerHashMap& hash) const { return apply_visitor(CKmerHashMap::index_get(*this), const_cast<CKmerHashMap&>(hash).m_hash_table); };
            MappedV* GetMapped(CKmerHashMap& hash) const { return apply_visitor(CKmerHashMap::index_get_mapped(*this), hash.m_hash_table); };
            const MappedV* GetMapped(const CKmerHashMap& hash) const { return apply_visitor(CKmerHashMap::index_get_mapped(*this), const_cast<CKmerHashMap&>(hash).m_hash_table); };
            const uint64_t* GetKeyPointer(const CKmerHashMap& hash) const { return apply_visitor(CKmerHashMap::index_get_keyp(*this), const_cast<CKmerHashMap&>(hash).m_hash_table); };
            bool operator==(const Index& other) const { return m_index == other.m_index && m_lstp == other.m_lstp; }
            bool operator!=(const Index& other) const { return !operator==(other); }
            bool operator<(const Index& other) const {
                if(m_index == other.m_index)
                    return m_lstp < other.m_lstp;
                else
                    return m_index < other.m_index;
            }
            bool operator>(const Index& other) const {
                if(m_index == other.m_index)
                    return m_lstp > other.m_lstp;
                else
                    return m_index > other.m_index;
            }
            struct Hash { size_t operator()(const Index& index) const { return std::hash<size_t>()(index.m_index)^std::hash<void*>()(index.m_lstp); } };
        protected:
            friend class CKmerHashMap;

            size_t m_index; // index; list considered a single entry
            void* m_lstp;   // pointer to list element
        };
        Index EndIndex() const { return Index((BucketBlock+1)*BucketsNum(), nullptr); }

        class Iterator : public Index {
        public:
            Iterator(const Index& index, CKmerHashMap* hp) : Index(index), hashp(hp) {}
            Iterator(size_t ind, void* ptr, CKmerHashMap* hp) : Index(ind, ptr), hashp(hp) {}
            Iterator& operator++() { 
                this->Advance(*hashp);
                return *this;
            }
            pair<TKmer, MappedV*> GetElement() { return Index::GetElement(*hashp); }
            MappedV* GetMapped() { return Index::GetMapped(*hashp); }
            const uint64_t* GetKeyPointer() { return Index::GetKeyPointer(*hashp); }
            size_t HashPosition() const { return Index::m_index; }
        private:
            CKmerHashMap* hashp;
        };
        Iterator Begin() { return Iterator(apply_visitor(hash_begin(0), m_hash_table), this); }
        Iterator End() { return Iterator((BucketBlock+1)*BucketsNum(), nullptr, this); }
        Iterator FirstForBucket(size_t bucket) { return Iterator(apply_visitor(hash_begin(bucket), m_hash_table), this); }
        vector<Iterator> Chunks(int desired_num) {
            vector<Iterator> chunks;
            if(BucketsNum() == 0)
                return chunks;

            size_t step = BucketsNum()/desired_num+1;
            for(size_t bucket = 0; bucket < BucketsNum(); ) {
                chunks.push_back(FirstForBucket(bucket));
                bucket = chunks.back().m_index/(BucketBlock+1)+step;
            }
            if(chunks.back() != End())
                chunks.push_back(End());

            return chunks;
        }

        //returns pointer to mapped value if exists, otherwise nullptr
        MappedV* Find(const TKmer& kmer) { 
            if(m_table_size == 0)
                return nullptr;
            else
                return apply_visitor(find(kmer), m_hash_table); 
        } 
        //returns index in hash table
        Index FindIndex(const TKmer& kmer) { 
            if(m_table_size == 0)
                return EndIndex();
            else
                return apply_visitor(find_index(kmer), m_hash_table); 
        }
        // if kmer already included returns pointer to mapped value
        // if not inserts a new entry and returns pointer to default value
        // caller MUST update the mapped value
        // assumes that any updates will be atomic
        MappedV* FindOrInsert(const TKmer& kmer) {
            size_t index = kmer.oahash()%m_table_size;
            return FindOrInsertInBucket(kmer, index);
        }
        MappedV* FindOrInsertInBucket(const TKmer& kmer, size_t index) {return apply_visitor(find_or_insert(kmer,index), m_hash_table); }

        MappedV* InitCell(const TKmer& kmer, size_t index) { return apply_visitor(init_cell(kmer,index), m_hash_table); }        

        void Swap(CKmerHashMap& other) {            
            apply_visitor(swap_with_other(), m_hash_table, other.m_hash_table);  
            swap(m_table_size, other.m_table_size);
            swap(m_kmer_len, other.m_kmer_len);
        }

        int KmerLen() const { return m_kmer_len; }
        size_t TableSize() const { return m_table_size; }
        size_t TableFootPrint() const { return apply_visitor(hash_footprint(), m_hash_table); }
        size_t BucketsNum() const { return m_table_size/max(1,BucketBlock); }
        void Info() const { apply_visitor(info(), m_hash_table); }
        void Save(ostream& os) const {
            os.write(reinterpret_cast<const char*>(&m_table_size), sizeof m_table_size);
            os.write(reinterpret_cast<const char*>(&m_kmer_len), sizeof m_kmer_len);
            apply_visitor(save(os), m_hash_table);
            if(!os)
              throw runtime_error("Error in CKmerHashMap write");  
        }
        void Load(istream& is) { 
            if(!is.read(reinterpret_cast<char*>(&m_table_size), sizeof m_table_size))
                throw runtime_error("Error in CKmerHashMap read");
            if(!is.read(reinterpret_cast<char*>(&m_kmer_len), sizeof m_kmer_len))
                throw runtime_error("Error in CKmerHashMap read");
            m_hash_table = CreateVariant<TKmerHashTable<MappedV>, THashBlockVec, MappedV>((m_kmer_len+31)/32);
            apply_visitor(load(is), m_hash_table);
        }

    protected:
        friend class Index;
        
        template<int N, class V> using THashBlockVec = CDeque<SHashBlock<LargeInt<N>,V,BucketBlock>>;
        template<class V> using TKmerHashTable = BoostVariant<THashBlockVec,V>;

        struct save : public boost::static_visitor<void> {
            save(ostream& out) : os(out) {}
            template <typename T> 
            void operator()(const T& v) const {
                v.Save(os);
                size_t list_num = 0;
                for(size_t i = 0; i < v.Size(); ++i) {
                    if(v[i].m_extra.Head() != nullptr)
                        ++list_num;
                }
                os.write(reinterpret_cast<const char*>(&list_num), sizeof list_num);
                for(size_t i = 0; i < v.Size(); ++i) {
                    if(v[i].m_extra.Head() != nullptr) {
                        os.write(reinterpret_cast<const char*>(&i), sizeof i);
                        v[i].m_extra.Save(os);
                    }
                }                    
            }
            ostream& os;
        };
        struct load : public boost::static_visitor<void> {
            load(istream& in) : is(in) {}
            template <typename T> 
            void operator()(T& v) const {
                v.Load(is);
                size_t list_num;
                if(!is.read(reinterpret_cast<char*>(&list_num), sizeof list_num))
                    throw runtime_error("Error in CKmerHashMap read");
                for( ; list_num > 0; --list_num) {
                    size_t i;
                    if(!is.read(reinterpret_cast<char*>(&i), sizeof i))
                        throw runtime_error("Error in CKmerHashMap read");
                    v[i].m_extra.Init();
                    v[i].m_extra.Load(is);
                }
            }
            istream& is;
        };

        struct swap_with_other : public boost::static_visitor<> {
            template <typename T> void operator() (T& a, T& b) const { a.Swap(b); }
            template <typename T, typename U> void operator() (T& a, U& b)  const { throw runtime_error("Can't swap different type containers"); }
        };

        struct index_get : public boost::static_visitor<pair<TKmer, MappedV*>> {
            index_get(const Index& ind) : index(ind) {}
            template <typename T> 
            pair<TKmer, MappedV*> operator()(T& v) const {
                auto elemp = v[index.m_index/(BucketBlock+1)].IndexGet(index.m_index%(BucketBlock+1), index.m_lstp); 
                return make_pair(TKmer(elemp->first), &(elemp->second));
            }
            const Index& index;
        };

        struct index_get_mapped : public boost::static_visitor<MappedV*> {
            index_get_mapped(const Index& ind) : index(ind) {}
            template <typename T> 
            MappedV* operator()(T& v) const {
                auto elemp = v[index.m_index/(BucketBlock+1)].IndexGet(index.m_index%(BucketBlock+1), index.m_lstp); 
                return &(elemp->second);
            }
            const Index& index;
        };

        struct index_get_keyp : public boost::static_visitor<const uint64_t*> {
            index_get_keyp(const Index& ind) : index(ind) {}
            template <typename T> 
            const uint64_t* operator()(T& v) const {
                auto elemp = v[index.m_index/(BucketBlock+1)].IndexGet(index.m_index%(BucketBlock+1), index.m_lstp); 
                return elemp->first.getPointer();
            }
            const Index& index;
        };

        template <typename T> 
        static Index next_available(T& v, size_t from) {
            for(size_t i = from; i < v.Size(); ++i) {
                auto& bucket = v[i];
                for(int shift = 0; shift < BucketBlock; ++shift) {
                    if(!bucket.isEmpty(shift)) 
                        return Index(i*(BucketBlock+1)+shift, nullptr);
                }
                if(bucket.m_extra.Head() != nullptr)
                    return Index(i*(BucketBlock+1)+BucketBlock, bucket.m_extra.Head());
            }
                
            return Index((BucketBlock+1)*v.Size(), nullptr);
        }
        
        struct index_advance : public boost::static_visitor<> {
            index_advance(Index& ind) : index(ind) {}
            template <typename T> 
            void operator()(T& v) const {
                typedef typename  T::value_type::list_t::SNode snode_t;

                size_t ind = index.m_index;
                size_t i = ind/(BucketBlock+1);
                auto& bucket = v[i];
                int shift = ind%(BucketBlock+1);
                if(shift < BucketBlock-1) {                    // not last array element - check all next elements
                    while(++shift < BucketBlock) {
                        if(!bucket.isEmpty(shift)) {
                            index.m_index = i*(BucketBlock+1)+shift;
                            return;
                        }
                    }
                } else if(shift == BucketBlock-1) {            // last array element - check spillover list
                    if(bucket.m_extra.Head() != nullptr) {
                        ++index.m_index;
                        index.m_lstp = bucket.m_extra.Head();
                        return;
                    }
                } else {                                       // spillover list - check next list element
                    snode_t* ptr = reinterpret_cast<snode_t*>(index.m_lstp);
                    if(ptr->m_next != nullptr) {
                        index.m_lstp = ptr->m_next;
                        return;
                    }
                }
                
                index = next_available(v, i+1); // look for next bucket
            }
            Index& index;
        };

        struct hash_begin : public boost::static_visitor<Index> {
            hash_begin(size_t fr) : from(fr) {}
            template <typename T> 
            Index operator()(T& v) const { return next_available(v, from); }
            size_t from;
        };

        struct hash_footprint : public boost::static_visitor<size_t> {
            template <typename T> size_t operator()(T& v) const { return sizeof(typename  T::value_type)*v.Size(); }
        };        
                
        //returns pointer to mapped value if exists, otherwise nullptr
        struct find : public boost::static_visitor<MappedV*> {
            find(const TKmer& k) : kmer(k) {}
            template <typename T> MappedV* operator()(T& v) const {
                auto& k = kmer.get<typename  T::value_type::large_t>();
                size_t pos = k.oahash()%(v.Size()*max(1,BucketBlock));
                auto& bucket = v[pos/max(1,BucketBlock)];
                int hint = pos%max(1,BucketBlock);               

                auto rslt = bucket.Find(k, hint);
                if(rslt.first < BucketBlock)                 // found in array
                    return &bucket.m_data[rslt.first].second;
                else if(rslt.first == BucketBlock)           // found in list
                    return &rslt.second->m_data.second;
                else                                         // not found
                    return nullptr; 
            } 
            const TKmer& kmer;
        };         
       
        //returns Index
        struct find_index : public boost::static_visitor<Index> {
            find_index(const TKmer& k) : kmer(k) {}
            template <typename T> Index operator()(T& v) const {
                typedef typename  T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();

                size_t pos = k.oahash()%(v.Size()*max(1,BucketBlock));
                size_t bucket_num = pos/max(1,BucketBlock);
                int hint = pos%max(1,BucketBlock);               

                auto rslt = v[bucket_num].Find(k, hint);
                if(rslt.first <= BucketBlock)                                         // found in array
                    return Index(bucket_num*(BucketBlock+1)+rslt.first, rslt.second);
                else                                                                  // not found
                    return Index((BucketBlock+1)*v.Size(), nullptr); 
            } 
            const TKmer& kmer;
        };        


        // if kmer already included returns pointer to mapped value
        // if not inserts a new entry and returns pointer to default value
        // caller MUST update the mapped value
        // assumes that any updated will be atomic
        struct find_or_insert : public boost::static_visitor<MappedV*> {
            find_or_insert(const TKmer& k, size_t i) : kmer(k), index(i) {}
            template <typename T> MappedV* operator()(T& v) const {
                typedef typename T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();
                size_t bucket_num = index/max(1,BucketBlock);
                int hint = index%max(1,BucketBlock);               
                return v[bucket_num].FindOrInsert(k, hint);
            }
            const TKmer& kmer;
            size_t index;
        };

        struct init_cell : public boost::static_visitor<MappedV*> {
            init_cell(const TKmer& k, size_t i) : kmer(k), index(i) {}
            template <typename T> MappedV* operator()(T& v) const {
                typedef typename T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();
                size_t bucket_num = index/(BucketBlock+1);
                int shift = index%(BucketBlock+1);
                return v[bucket_num].InitCell(k, shift);
            }
            const TKmer& kmer;
            size_t index;
        };


        struct info : public boost::static_visitor<> {
            template <typename T> void operator()(T& v) const {
                
                map<int,int> numbers;
                for(size_t i = 0; i < v.Size(); ++i) {
                    auto& bucket= v[i];
                    int num = distance(bucket.m_extra.begin(), bucket.m_extra.end());  
                    for(int shift = 0; shift < BucketBlock; ++shift) {
                        if(!bucket.isEmpty(shift))
                            ++num;
                    }                                        
                    ++numbers[num];
                }
                for(auto& rslt : numbers)
                    cerr << "Bucket:\t" << rslt.first << "\t" << rslt.second << endl;
            }            
        };

        struct resize : public boost::static_visitor<> { 
            resize(size_t s, size_t c) : size(s), chunks(c) {}
            template <typename T> void operator()(T& v) const { v.Reset(size, chunks); }            
            size_t size;
            size_t chunks;
        };

        TKmerHashTable<MappedV> m_hash_table;
        size_t m_table_size;
        int m_kmer_len;
    };

    struct SKmerCounter {
        SKmerCounter() : m_data(0) {}
        bool operator==(const SKmerCounter& kc) const { return kc.m_data == m_data; }
        uint32_t Increment(bool is_plus) { return (m_data.m_atomic += (is_plus ? 0x100000001 : 1)); }
        uint32_t Count() const { return m_data; } // clips plus part
        SAtomic<uint64_t> m_data;
    };
    
    class CKmerHashCount : public CKmerHashMap<SKmerCounter, 8> {
    public:
        CKmerHashCount(int kmer_len = 0, size_t size = 0) : CKmerHashMap(kmer_len, size) {}
        // returns true if kmer was new
        bool UpdateCount(const TKmer& kmer, bool is_plus) {
            size_t index = kmer.oahash()%m_table_size;
            return (FindOrInsertInBucket(kmer, index)->Increment(is_plus) == 1);
        } 
        size_t UpdateCounts(const CReadHolder::string_iterator& is, const TConcurrentBlockedBloomFilter& bloom, int min_count) {
            return apply_visitor(update_counts(is, bloom, min_count, m_kmer_len), m_hash_table);
        }
        // rehash bucket from other container
        void RehashOtherBuckets(CKmerHashCount& other, size_t bucket_from, size_t bucket_to) {
            apply_visitor(rehash_bucket(bucket_from, bucket_to, *this), m_hash_table, other.m_hash_table); 
        }
        //remove false positives
        size_t CleanBuckets(int min_count, size_t bucket_from, size_t bucket_to) {
            return apply_visitor(clean_buckets(min_count, bucket_from, bucket_to, TableSize()), m_hash_table);
        }
        TBins GetBins() {
            map<int,size_t> hist;
            for(auto index = Begin(); index != End(); ++index) {
                ++hist[index.GetMapped()->Count()];            
            }
            return TBins(hist.begin(), hist.end());
        }
    private:
        struct update_counts : public boost::static_visitor<size_t> {
            update_counts(const CReadHolder::string_iterator& i, const TConcurrentBlockedBloomFilter& bl, int mc, unsigned kl) : is(i), bloom(bl), min_count(mc), kmer_len(kl) {}
            template <typename T> size_t operator() (T& v) const {
                if(v.Size() == 0)
                    return 0;

                typedef typename T::value_type::large_t large_t;

                size_t read_len = is.ReadLen();
                if(read_len < kmer_len)
                    return 0;

                unsigned kmer_bytes = (2*kmer_len+7)/8;                 //number of whole bytes in kmer
                unsigned kmer_size = (2*kmer_len+63)/64;                //number of whole 8-byte words in kmer
                int partial_bits = (2*kmer_len)%64;                     //number of used bits in partial 8 byte word (if any)
                uint64_t mask = numeric_limits<uint64_t>::max();
                if(partial_bits > 0)
                    mask = (uint64_t(1) << partial_bits) - 1;
                size_t buf_size = (2*read_len+63)/64+1;
                uint64_t* read_buf = new uint64_t[buf_size]; //(enough + 1) 8-byte words for read (one extra because we'll copy kmers using whole bytes which can go beyond the sequence)

                size_t new_kmers = 0;
                large_t kmer(0);
                for(int shift = 0; shift < 4 && read_len-shift >= kmer_len; ++shift) {
                    memset(read_buf, 0, 8*buf_size);
                    is.BSeq(shift, read_buf);
                    for(unsigned k = 0; k <= read_len-shift-kmer_len; k += 4) { // every 4th kmer on the byte boundary
                        memcpy(kmer.getPointer(), (uint8_t*)read_buf+k/4, kmer_bytes);
                        kmer.getPointer()[kmer_size-1] &= mask;
                        large_t rkmer = revcomp(kmer, kmer_len);
                        large_t* min_kmerp = &rkmer; 
                        bool is_plus = false;
                        size_t hashp = rkmer.oahash();
                        size_t hashm = kmer.oahash();
                        if(kmer < rkmer) {
                            min_kmerp = &kmer;
                            is_plus = true;
                            swap(hashp, hashm);
                        }

                        int bucket_block = T::value_type::eBucketBlock;
                        size_t pos = hashp%(v.Size()*max(1,bucket_block));
                        size_t bucket_num = pos/max(1,bucket_block);
                        int hint = pos%max(1,bucket_block);
                        auto& bucket = v[bucket_num];
                        
                        auto rslt = bucket.Find(*min_kmerp, hint);
                        if(rslt.first < bucket_block) {                          // found in array
                            if(bucket.m_data[rslt.first].second.Increment(is_plus) == 1)
                                ++new_kmers;
                            continue;
                        } else if(rslt.first == bucket_block) {                  // found in list
                            if(rslt.second->m_data.second.Increment(is_plus) == 1)
                                ++new_kmers;
                            continue;
                        }                        
                        
                        if(min_count > 1 && bloom.Test(hashp, hashm) < min(min_count, bloom.MaxElement()))
                            continue;  

                        if(bucket.FindOrInsert(*min_kmerp, hint)->Increment(is_plus) == 1)
                            ++new_kmers;
                    }
                }

                delete[] read_buf;
                return new_kmers;
            }
            const CReadHolder::string_iterator& is;
            const TConcurrentBlockedBloomFilter& bloom;
            int min_count;
            unsigned kmer_len;
        };

        struct rehash_bucket : public boost::static_visitor<> {
            rehash_bucket(size_t bf, size_t bt, CKmerHashCount& h) : bucket_from(bf), bucket_to(bt), hash(h) {}
            template <typename T> void operator() (T& a, T& b) const {
                typedef typename T::value_type::element_t element_t;
                for(size_t indexb = bucket_from; indexb <= bucket_to; ++indexb) {
                    auto& bucket_b = b[indexb];

                    for(auto& cell : bucket_b.m_data) {
                        if(cell.second.Count() != 0) {
                            auto& kmer = cell.first;
                            size_t indexa = kmer.oahash()%hash.TableSize();
                            *hash.FindOrInsertInBucket(TKmer(kmer), indexa) = cell.second;
                        }
                    }
                    
                    for(element_t& cell : bucket_b.m_extra) {
                        auto& kmer = cell.first;
                        size_t indexa = kmer.oahash()%hash.TableSize();
                        *hash.FindOrInsertInBucket(TKmer(kmer), indexa) = cell.second;
                    }
                }
            }
            template <typename T, typename U> void operator() (T& a, U& b) const { throw runtime_error("Can't rehash from different type container"); }
            
            size_t bucket_from;
            size_t bucket_to;
            CKmerHashCount& hash;
        };

        struct clean_buckets : public boost::static_visitor<size_t> {
            clean_buckets(int mc, size_t bf, size_t bt, size_t tb) : min_count(mc), bucket_from(bf), bucket_to(bt), table_size(tb) {}            
            template <typename T> size_t operator()(T& v) const {
                typedef typename T::value_type::element_t element_t;

                size_t num = 0;
                for(size_t bind = bucket_from; bind <= bucket_to; ++bind) {
                    auto& bucket = v[bind];
                    int empty_cells = 0;

                    auto Reposition = [&bucket, this](element_t& cell, unsigned limit) {
                        size_t index = cell.first.oahash()%table_size;
                        size_t orig_shift = index%bucket.m_data.size();
                        if(orig_shift == limit)
                            return false;

                        if(bucket.m_data[orig_shift].second.Count() < min_count) {
                            if(limit < bucket.m_data.size())
                                bucket.Move(limit, orig_shift);
                            else
                                bucket.Move(cell, orig_shift);

                            return orig_shift > limit;
                        } 

                        for(unsigned shift = 0; shift < limit; ++shift) {
                            if(shift != orig_shift && bucket.m_data[shift].second.Count() < min_count) {
                                if(limit < bucket.m_data.size())
                                    bucket.Move(limit, shift);
                                else
                                    bucket.Move(cell, shift);

                                return false;
                            }                         
                        }

                        return false;
                    };

                                      
                    for(unsigned shift = 0; shift < bucket.m_data.size(); ++shift) {
                        auto& cell = bucket.m_data[shift];
                        auto count = cell.second.Count();
                        if(count < min_count) {
                            ++empty_cells;                            
                            if(count > 0)
                                bucket.Clear(shift);
                        } else {
                            if(Reposition(cell, shift))
                                ++empty_cells;          // moved down and created new empty cell (will be counted later)
                            else
                                ++num;                  // stayed or moved up
                        }
                    }
 
                    for(auto& cell : bucket.m_extra) {
                        if(cell.second.Count() >= min_count) {
                            ++num;
                            if(empty_cells > 0) {
                                Reposition(cell, bucket.m_data.size());
                                --empty_cells;
                            }
                        }
                    }
                                         
                    bucket.m_extra.remove_if([this](const element_t& elem) {return elem.second.Count() < min_count;});
                 }

                return num;
            }            
            unsigned min_count;
            size_t bucket_from;
            size_t bucket_to;
            size_t table_size;
        };
    };


    class CKmerHashCounter {
    public:
        CKmerHashCounter(const list<array<CReadHolder,2>>& reads, int kmer_len, int min_count, size_t estimated_kmer_num, bool is_stranded, int ncores, bool skip_bloom) : 
            m_kmer_len(kmer_len), m_min_count(min_count), m_is_stranded(is_stranded), m_ncores(ncores), m_skip_bloom(skip_bloom), m_hash_table(m_kmer_len), 
            m_estimated_table_size(0),  m_estimated_uniq_kmers(0), m_kmer_num(0), m_kmer_num_raw(0), m_kmer_count(0), m_rehash_status(false) {

            m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
            for(auto& rholder : reads) 
                m_start_position.push_back(make_pair(0, rholder[0].sbegin()));
            
            CStopWatch timer;

            TConcurrentBlockedBloomFilter bloom(0, 2, 1, m_min_count);
            if(!m_skip_bloom)  {
                m_estimated_uniq_kmers.store(estimated_kmer_num);
                while(true) {
                    timer.Restart();
                    int counter_bit_size = 2;
                    for( ; counter_bit_size < 8 &&  (1 << counter_bit_size)-1 < m_min_count; counter_bit_size *= 2);
                    double false_positive_rate = 0.03;
                    //                    size_t bloom_table_size = -1.5*(double)m_estimated_uniq_kmers.load()*log(false_positive_rate)/log(2.)/log(2.); // 50% extra because blocked
                    size_t bloom_table_size = -(double)m_estimated_uniq_kmers.load()*log(false_positive_rate)/log(2.)/log(2.); // 50% extra because blocked
                    int hash_num = ceil(-log(false_positive_rate)/log(2.));
                    bloom.Reset(bloom_table_size, counter_bit_size, hash_num, m_min_count);

                    cerr << "\nBloom table size: " << bloom.TableSize() << "(" << 0.1*(bloom.TableFootPrint()/100000) << "MB)" << " Counter bit size: " << counter_bit_size << " Hash num: " << hash_num << endl;  

                    m_estimated_table_size.store(0);
                    m_estimated_uniq_kmers.store(0);

                    list<function<void()>> jobs;
                    for(auto& job_input : reads) {
                        if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty                   
                            jobs.push_back(bind(&CKmerHashCounter::InsertInBloomJob, this, ref(job_input), ref(bloom)));
                        }
                    }
                    RunThreads(m_ncores, jobs);            
            
                    if(m_min_count == 1)
                        m_estimated_table_size.store(m_estimated_uniq_kmers.load());
                    double kmers = m_estimated_uniq_kmers.load();
                    false_positive_rate = pow(1.-exp(-hash_num*kmers/bloom_table_size), hash_num);
                    cerr << "Estimated kmers above threshold: " << m_estimated_table_size.load() << " Estimated uniq kmers: " << m_estimated_uniq_kmers.load() << " Estimated bloom false positive rate " << false_positive_rate << endl;
                    cerr << "Bloom filter in " << timer.Elapsed();
                    
                    if(false_positive_rate < 0.15)
                        break;

                    cerr << "\nBloom filter false positive rate is too high - increasing the bloom filter size and recalculating" << endl;
                }
            } else {
                m_estimated_table_size.store(estimated_kmer_num);
            }

            timer.Restart();
            m_hash_table.Reset(1.5*m_estimated_table_size.load(), m_ncores);            
            while(m_hash_table.TableSize() > 0) {
                {
                    CStopWatch timer;
                    timer.Restart();
                    list<function<void()>> jobs;
                    auto start_pos = m_start_position.begin();
                    for(auto& job_input : reads) {
                        if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty               
                            jobs.push_back(bind(&CKmerHashCounter::CountKmersJob, this, ref(job_input), ref(*start_pos), ref(bloom)));
                        }
                        ++start_pos;
                    }
                    RunThreads(m_ncores, jobs);
               }
                
                if(!m_rehash_status.load())
                    break;

                //Rehash
                {
                    CStopWatch timer;
                    timer.Restart();
                    size_t new_size = m_hash_table.TableSize()*m_increase_factor;
                    cerr << "Rehash new size: " << new_size << endl;
                    CKmerHashCount hash_table_tmp(m_kmer_len);
                    hash_table_tmp.Reset(new_size, m_ncores);
                    swap(m_hash_table, hash_table_tmp);
                    m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
                    m_rehash_status.store(false);
                
                    list<function<void()>> jobs;
                    size_t step = ceil((double)hash_table_tmp.BucketsNum()/m_ncores);
                    for(int thr = 0; thr < m_ncores; ++thr) {
                        size_t from = step*thr;
                        size_t to = min(hash_table_tmp.BucketsNum()-1,from+step-1);
                        if(to >= from)
                            jobs.push_back(bind(&CKmerHashCounter::RehashJob, this, ref(hash_table_tmp), from, to));
                    }
                    RunThreads(m_ncores, jobs);
                    cerr << "Rehashing in " << timer.Elapsed();
                }
            }

            cerr << "Create hash in " << timer.Elapsed();
            timer.Restart();            
            //Remove false positives
            RemoveLowCountKmers(m_min_count);
            cerr << "Clean hash in " << timer.Elapsed();
            timer.Restart();            

            cerr << "Initial kmers: " << m_kmer_num_raw.load() << " Kmers above threshold: " << m_kmer_num.load() << " Total kmers: " << m_kmer_count.load() << " Hash table size: " <<  m_hash_table.TableSize() << "(" << 0.1*(m_hash_table.TableFootPrint()/100000) << "MB)" << endl;

        }
        void RemoveLowCountKmers(int min_count) {
            m_min_count = min_count;
            if(m_hash_table.TableSize() > 0) {
                list<function<void()>> jobs;
                size_t step = ceil((double)m_hash_table.BucketsNum()/m_ncores);
                for(int thr = 0; thr < m_ncores; ++thr) {
                    size_t from = step*thr;
                    size_t to = min(m_hash_table.BucketsNum()-1,from+step-1);
                    if(to >= from)
                        jobs.push_back(bind(&CKmerHashCounter::CleanJob, this, from, to));
                }
                RunThreads(m_ncores, jobs);
            }
        }                                           
        // prepares kmer counts to be used in  de Bruijn graph
        void GetBranches() {
            CStopWatch timer;
            timer.Restart();

            list<function<void()>> jobs;
            size_t step = ceil((double)m_hash_table.BucketsNum()/m_ncores);
            for(int thr = 0; thr < m_ncores; ++thr) {
                size_t from = step*thr;
                size_t to = min(m_hash_table.BucketsNum()-1,from+step-1);
                if(to >= from)
                    jobs.push_back(bind(&CKmerHashCounter::GetBranchesJob, this, from, to));
            }
            RunThreads(m_ncores, jobs);

            cerr << "Kmers branching in " << timer.Elapsed();
        }
        void Info() const {
            m_hash_table.Info();
        }
        CKmerHashCount& Kmers() { return m_hash_table; }
        size_t KmerNum() const { return m_kmer_num; }
        

    private:

        void GetBranchesJob(size_t bucket_from, size_t bucket_to) {
            TKmer max_kmer(string(m_kmer_len, bin2NT[3]));
            for(size_t bucket = bucket_from; bucket <= bucket_to; ++bucket) {
                for(auto index = m_hash_table.FirstForBucket(bucket); index != m_hash_table.FirstForBucket(bucket+1); ++index) {
                    uint64_t branches = 0;
                    pair<TKmer, SKmerCounter*> kmer_count = index.GetElement();
                    //direct        
                    TKmer shifted_kmer = (kmer_count.first << 2) & max_kmer;
                    //inverse       
                    TKmer shifted_rkmer = (revcomp(kmer_count.first, m_kmer_len) << 2) & max_kmer;
                    for(int nt = 0; nt < 4; ++nt) {
                        TKmer k = shifted_kmer + TKmer(m_kmer_len, nt);
                        SKmerCounter* nbrp = m_hash_table.Find(min(k, revcomp(k, m_kmer_len)));
                        // New kmer is a neighbor if it exists in reads and is not same as current kmer
                        if(nbrp != nullptr && nbrp != kmer_count.second)
                            branches |= (1 << nt);

                        k = shifted_rkmer + TKmer(m_kmer_len, nt);
                        nbrp = m_hash_table.Find(min(k, revcomp(k, m_kmer_len)));
                        if(nbrp != nullptr && nbrp != kmer_count.second)
                            branches |= (1 << (nt+4));
                    }

                    uint64_t count = kmer_count.second->m_data;
                    uint32_t total_count = count;
                    uint32_t plus_count = (count >> 32);
                    uint64_t plusf = uint16_t(double(plus_count)/total_count*numeric_limits<uint16_t>::max()+0.5);
                    kmer_count.second->m_data = (plusf << 48)+(branches << 32)+total_count;
                }
            }
        }
        
        void CleanJob(size_t bucket_from, size_t bucket_to) {
             m_kmer_num += m_hash_table.CleanBuckets(m_min_count, bucket_from, bucket_to);
        }

        class CBloomInserter : public TKmer {
        public:
            CBloomInserter(int kmer_len) : TKmer(kmer_len, 0), m_kmer_len(kmer_len) {}
            pair<size_t,size_t> InsertInBloom(const CReadHolder::string_iterator& is, TConcurrentBlockedBloomFilter& bloom) {
                return apply_visitor(insert_in_bloom(is, bloom, m_kmer_len), v);
            }

        private:
            unsigned m_kmer_len;

            struct insert_in_bloom : public boost::static_visitor<pair<size_t,size_t>> {
                insert_in_bloom(const CReadHolder::string_iterator& i, TConcurrentBlockedBloomFilter& bl, unsigned kl) : is(i), bloom(bl), kmer_len(kl) {}
                template <typename large_t> pair<size_t,size_t> operator() (large_t& kmer) const {
                    size_t above_threshold_kmers = 0;
                    size_t uniq_kmers = 0;

                    size_t read_len = is.ReadLen();
                    if(read_len < kmer_len)
                        return make_pair(above_threshold_kmers, uniq_kmers);

                    unsigned kmer_bytes = (2*kmer_len+7)/8;                 //number of whole bytes in kmer
                    unsigned kmer_size = (2*kmer_len+63)/64;                //number of whole 8-byte words in kmer
                    int partial_bits = (2*kmer_len)%64;                     //number of used bits in partial 8 byte word (if any)
                    uint64_t mask = numeric_limits<uint64_t>::max();
                    if(partial_bits > 0)
                        mask = (uint64_t(1) << partial_bits) - 1;
                    size_t buf_size = (2*read_len+63)/64+1;
                    uint64_t* read_buf = new uint64_t[buf_size]; //(enough + 1) 8-byte words for read (one extra because we'll copy kmers using whole bytes which can go beyond the sequence)

                    for(int shift = 0; shift < 4 && read_len-shift >= kmer_len; ++shift) {
                        memset(read_buf, 0, 8*buf_size);
                        is.BSeq(shift, read_buf);
                        for(unsigned k = 0; k <= read_len-shift-kmer_len; k += 4) { // every 4th kmer on the byte boundary
                            memcpy(kmer.getPointer(), (uint8_t*)read_buf+k/4, kmer_bytes);
                            kmer.getPointer()[kmer_size-1] &= mask;
                            large_t rkmer = revcomp(kmer, kmer_len);

                            size_t hashp = rkmer.oahash();
                            size_t hashm = kmer.oahash();
                            if(kmer < rkmer)
                                swap(hashp, hashm);  

                            switch(bloom.Insert(hashp, hashm)) {
                            case TConcurrentBlockedBloomFilter::eNewKmer : ++uniq_kmers; continue;
                            case TConcurrentBlockedBloomFilter::eAboveThresholdKmer : ++above_threshold_kmers; continue;
                            default : continue;
                            }
                        }
                    }
                    delete[] read_buf;
                    return make_pair(above_threshold_kmers, uniq_kmers);
                }                
                const CReadHolder::string_iterator& is;
                TConcurrentBlockedBloomFilter& bloom;
                unsigned kmer_len;
            };
        };

        void InsertInBloomJob(const array<CReadHolder,2>& rholder, TConcurrentBlockedBloomFilter& bloom) {
            size_t above_threshold_kmers = 0;
            size_t uniq_kmers = 0;
 
            CBloomInserter bloom_inserter(m_kmer_len);
            for(int p = 0; p < 2; ++p) {
               for(CReadHolder::string_iterator is = rholder[p].sbegin(); is != rholder[p].send(); ++is) {
                   auto rslt = bloom_inserter.InsertInBloom(is, bloom);
                   above_threshold_kmers += rslt.first;
                   uniq_kmers += rslt.second;
               }
            }

            m_estimated_table_size += above_threshold_kmers;
            m_estimated_uniq_kmers += uniq_kmers;
        }
        void RehashJob(CKmerHashCount& other_hash_table, size_t bucket_from, size_t bucket_to) {
            m_hash_table.RehashOtherBuckets(other_hash_table, bucket_from, bucket_to);
        }

        void CountKmersJob(const array<CReadHolder,2>& rholder, pair<int,CReadHolder::string_iterator>& start_pos, const TConcurrentBlockedBloomFilter& bloom) {
            size_t kmer_num = 0;
            size_t kmer_count = 0;
            for(int p = start_pos.first; p < 2; ++p) {
                CReadHolder::string_iterator from = start_pos.second;
                if(p != start_pos.first)
                    from = rholder[p].sbegin();
                for(CReadHolder::string_iterator is = from; is != rholder[p].send(); ++is) {
                    size_t read_len = is.ReadLen();
                    if(read_len >= (unsigned)m_kmer_len)
                        kmer_count += read_len-m_kmer_len+1;
                    else
                        continue;

                    kmer_num += m_hash_table.UpdateCounts(is, bloom, m_skip_bloom ? 0 : m_min_count);

                    if(kmer_num >= m_kmer_step) {
                        m_kmer_num_raw += kmer_num;
                        m_kmer_count += kmer_count;
                        kmer_num = 0;
                        kmer_count = 0;
                        if(m_kmer_num_raw.load() > m_hash_table.TableSize()*m_max_load_factor)
                            m_rehash_status.store(true); 
                        if(m_rehash_status.load()) {
                            start_pos.first = p;
                            ++is;
                            start_pos.second = is;
                            return;
                        }
                    }                    
                }
            }
            m_kmer_num_raw += kmer_num;            
            m_kmer_count += kmer_count;
            start_pos.first = 2;
        }

        int m_kmer_len;
        int m_min_count;
        bool m_is_stranded;
        int m_ncores;
        bool m_skip_bloom;

        CKmerHashCount m_hash_table;
        atomic<size_t> m_estimated_table_size;
        atomic<size_t> m_estimated_uniq_kmers;
        atomic<size_t> m_kmer_num;
        atomic<size_t> m_kmer_num_raw;
        atomic<size_t> m_kmer_count;
        atomic<bool> m_rehash_status;
        size_t m_kmer_step;
        double m_max_load_factor = 1;
        int m_increase_factor = 2;
        list<pair<int,CReadHolder::string_iterator>> m_start_position;
    };


}; // namespace
#endif /*_Concurrent_Hash_*/

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

#ifndef _common_util_
#define _common_util_
#include <atomic>
#include <future>
#include <thread>
#include <boost/timer/timer.hpp>
#include <cmath>

using namespace std;
namespace DeBruijn {

    // Wraps around atomic<> to make it possible to use in containers
    // IMPORTANT: don't concurrently create or modify containers of SAtomic!!!!!
    template <typename T>
    struct SAtomic {
        typedef T Type;
        SAtomic(T t = 0) { m_atomic.store(t); }
        SAtomic(const atomic<T> &a) { m_atomic.store(a.load()); }
        SAtomic(const SAtomic &other) { m_atomic.store(other.m_atomic.load()); }
        SAtomic& operator=(const SAtomic &other) { 
            m_atomic.store(other.m_atomic.load());
            return *this;
        }        
        SAtomic& operator=(T t) { 
            m_atomic.store(t); 
            return *this;
        } 
        bool Set(T value, T expected = 0) { return m_atomic.compare_exchange_strong(expected, value); }
        operator T() const { return m_atomic.load(); }
        T Load() const { return m_atomic.load(); }

        atomic<T> m_atomic;
    };

    class CStopWatch : public boost::timer::cpu_timer {
    public:
        void Restart() { start(); }
        string Elapsed() const { return format(); }
        void Stop() { stop (); }
        void Resume() { resume(); }
    };


    // runs ncores threads until all jobs are exhausted
    void RunThreads(int ncores, list<function<void()>>& jobs) {
        typedef list<future<void>> ThreadsStatus;
        ThreadsStatus active_threads_status;

        //        int total_jobs = jobs.size();
        //        cerr << "Remaining " << total_jobs << " jobs from " << total_jobs << endl;

        //create ncores threads
        for(int i = 0; i < ncores && !jobs.empty(); ++i) {
            active_threads_status.push_front(async(launch::async, jobs.front()));
            jobs.pop_front();
        }

        //for each finished thread create a new one until done  
        chrono::milliseconds span (1);
        while(!active_threads_status.empty()) {
            for(auto iloop = active_threads_status.begin(); iloop != active_threads_status.end(); ) {
                auto done = iloop++;
                if(done->wait_for(span) == future_status::timeout)   // not ready    
                    continue;                

                done->get();
                active_threads_status.erase(done);
                if(!jobs.empty()) {
                    active_threads_status.push_front(async(launch::async, jobs.front()));
                    jobs.pop_front();
                }
                //                cerr << "Remaining jobs " << jobs.size()+active_threads_status.size() << " from " << total_jobs << endl;    
            }
        }
        //        cerr << endl; 
    }

    // Stores DNA sequences using 4 letter alphabet
    // The sequences and kmers could be accessed sequentially using iterator-type classes
    // 
    class CReadHolder {
    public:
        CReadHolder(bool contains_paired) :  m_total_seq(0), m_contains_paired(contains_paired) {};

        // inserts read at the end
        template <typename Container>
        void PushBack(const Container& read) {
            int shift = (m_total_seq*2)%64;
            int read_len = 0;
            for(auto it = read.rbegin(); it != read.rend(); ++it) {   // put backward for kmer compatibility
                if(shift == 0)
                    m_storage.push_back(0);
                m_storage.back() += ((find(bin2NT.begin(), bin2NT.end(),  *it) - bin2NT.begin()) << shift);
                shift = (shift+2)%64;
                ++read_len;
            }
            m_read_length.push_back(read_len);
            m_total_seq += read_len;
        }

        template <typename RandomIterator>
        void PushBack(RandomIterator begin, uint32_t len) {
            int shift = (m_total_seq*2)%64;
            for(RandomIterator it = begin+len-1; ; --it) {
                if(shift == 0)
                    m_storage.push_back(0);
                m_storage.back() += ((find(bin2NT.begin(), bin2NT.end(),  *it) - bin2NT.begin()) << shift);
                shift = (shift+2)%64;

                if(it == begin)
                    break;
            }
            m_read_length.push_back(len);
            m_total_seq += len;
        }

        // insert sequence from other container
        class string_iterator;
        void PushBack(const string_iterator& is) {
            size_t read_len = is.ReadLen();
            m_read_length.push_back(read_len);
            size_t destination_first_bit = 2*m_total_seq;
            m_total_seq += read_len;
            m_storage.resize((2*m_total_seq+63)/64);

            const CReadHolder& other_holder = *is.m_readholderp;
            size_t bit_from = is.m_position;
            size_t bit_to = bit_from+2*read_len;
            other_holder.CopyBits(bit_from, bit_to, m_storage, destination_first_bit, m_storage.size());
        }

        // swaps contents with other
        void Swap(CReadHolder& other) {
            swap(m_storage, other.m_storage);
            swap(m_read_length, other.m_read_length);
            swap(m_total_seq, other.m_total_seq);
        }

        // deletes all sequences and releases  memory
        void Clear() { CReadHolder(m_contains_paired).Swap(*this); }

        // Total nucleotide count of the sequnce
        size_t TotalSeq() const { return m_total_seq; }

        // Maximal length of included sequences
        size_t MaxLength() const { 
            if(m_read_length.empty())
                return 0;
            else
                return *max_element(m_read_length.begin(), m_read_length.end()); 
        }

        // the number of kmers of give length that could be generated
        size_t KmerNum(unsigned kmer_len) const {
            size_t num = 0;
            if(m_read_length.empty())
                return num;
            for(auto l : m_read_length) {
                if(l >= kmer_len)
                    num += l-kmer_len+1;
            }
            return num;
        }

        // total number of sequences
        size_t ReadNum() const { return m_read_length.size(); }

        size_t MemoryFootprint() const { return 8*m_storage.capacity()+4*m_read_length.capacity(); }  // memory in bytes
        void Reserve(size_t seq, size_t num = 0) {
            m_storage.reserve(seq/32+1);
            if(num > 0)
                m_read_length.reserve(num);
        }

        // shortest sequence length at xx% of total length
        size_t NXX(double xx) const {
            vector<uint32_t> read_length(m_read_length.begin(), m_read_length.end());
            sort(read_length.begin(), read_length.end());
            size_t nxx = 0;
            size_t len = 0;
            for(int j = (int)read_length.size()-1; j >= 0 && len < xx*m_total_seq; --j) {
                nxx = read_length[j];
                len += read_length[j];
            }
            
            return nxx;
        }
        // shortest sequence length at 50% of total length
        size_t N50() const { return NXX(0.5); }

        // iterator-type clas to access kmers
        class kmer_iterator;
        kmer_iterator kend() const { return kmer_iterator(0, *this, 2*m_total_seq); }
        kmer_iterator kbegin(int kmer_len) const { return kmer_iterator(kmer_len, *this); }

        class kmer_iterator {
        public:

            // dereference operator; returns value!
            TKmer operator*() const {
                TKmer kmer(m_kmer_len, 0);
                uint64_t* guts = kmer.getPointer();
                size_t bit_from = m_position;
                size_t bit_to = bit_from+2*m_kmer_len;
                m_readholderp->CopyBits(bit_from, bit_to, guts, 0, (2*m_kmer_len+63)/64);
                
                return kmer;
            }

            // iterator advance
            kmer_iterator& operator++() {
                if(m_position == 2*(m_readholderp->m_total_seq-m_kmer_len)) {
                    m_position = 2*m_readholderp->m_total_seq;
                    return *this;
                }

                m_position += 2;
                if(++m_position_in_read == m_readholderp->m_read_length[m_read]-m_kmer_len+1) {
                    m_position += 2*(m_kmer_len-1);
                    ++m_read;
                    m_position_in_read = 0;
                    SkipShortReads();
                } 

                return *this;
            }
            // doesn't check read boundaries - should be used only if landing in the SAME read
            kmer_iterator& operator+=(int l) {
                m_position += 2*l;
                m_position_in_read += l;

                return *this;
            }

            friend bool operator==(kmer_iterator const& li, kmer_iterator const& ri) { return li.m_position == ri.m_position && li.m_readholderp == ri.m_readholderp; }
            friend bool operator!=(kmer_iterator const& li, kmer_iterator const& ri) { return li.m_position != ri.m_position || li.m_readholderp != ri.m_readholderp; }
            friend class CReadHolder;

        private:
            kmer_iterator(int kmer_len, const CReadHolder& rholder, size_t position = 0, size_t position_in_read = 0, size_t read = 0) : m_readholderp(&rholder), m_read(read), m_position(position), m_kmer_len(kmer_len), m_position_in_read(position_in_read) {
                SkipShortReads();
            }

            void SkipShortReads() {
                while(m_position < 2*m_readholderp->m_total_seq && m_read < m_readholderp->m_read_length.size() && m_readholderp->m_read_length[m_read] < m_kmer_len)
                    m_position += 2*m_readholderp->m_read_length[m_read++];                               
            }
            const CReadHolder* m_readholderp;
            size_t m_read;               // read number
            size_t m_position;           // BIT num in concatenated sequence
            uint32_t m_kmer_len;
            uint32_t m_position_in_read; // SYMBOL in read
        };

        // iterator-type clas to access reads
        string_iterator send() const { return string_iterator(*this, 2*m_total_seq, m_read_length.size()); }
        string_iterator sbegin() const { return string_iterator(*this); }

        enum {eSingle = 0, eFirstMate = 1, eSecondMate = 2};
        class string_iterator {
        public:
            string_iterator() : m_readholderp(nullptr), m_position(0), m_read(0) {}

            string operator*() const {
                int read_length = m_readholderp->m_read_length[m_read];
                string read;
                read.reserve(read_length);
                size_t position = m_position+2*(read_length-1);
                for(int i = 0; i < read_length; ++i) {
                    read.push_back(bin2NT[(m_readholderp->m_storage[position/64] >> position%64) & 3]);
                    position -= 2;
                }
                return read;
            }
            // returns inversed binary sequence (not complemented) 
            // assumes that destination is extended properly and filled with 0s    
            void BSeq(int shift, uint64_t* destination) const {
                size_t position = m_position+2*shift;
                size_t len = 2*(ReadLen()-shift);
                m_readholderp->CopyBits(position, position+len, destination, 0, (len+63)/64);
            }
            // returns clipped binary sequence in correct order
            // assumes that destination is extended properly and filled with 0s  
            // left/right refer to the original sequence  
            void TrueBSeq(size_t left_clip, size_t right_clip, bool reverse_complement, uint64_t* destination) const {
                auto Reverse = [](uint64_t& word) {
                    word = ((word & 0x3333333333333333) << 2)  | ((word >> 2)  & 0x3333333333333333); // swap adjacent pairs
                    word = ((word & 0x0F0F0F0F0F0F0F0F) << 4)  | ((word >> 4)  & 0x0F0F0F0F0F0F0F0F); // swap nibbles
                    word = ((word & 0x00FF00FF00FF00FF) << 8)  | ((word >> 8)  & 0x00FF00FF00FF00FF); // swap bytes
                    word = ((word & 0x0000FFFF0000FFFF) << 16) | ((word >> 16) & 0x0000FFFF0000FFFF); // swap 16 bit chunks
                    word = ((word & 0x00000000FFFFFFFF) << 32) | ((word >> 32) & 0x00000000FFFFFFFF); // swap 32 bit chunks                                        
                }; 

                size_t position = m_position+2*right_clip;  // sequence stored reversed
                size_t len = 2*(ReadLen()-right_clip-left_clip);
                size_t destination_size = (len+63)/64;
                
                if(reverse_complement) {
                    m_readholderp->CopyBits(position, position+len, destination, 0, destination_size);                   // already reversed; not complemented
                    for(size_t p = 0; p < destination_size; ++p)  // complement (will also convert trailing As into Ts)
                        destination[p] ^= 0xAAAAAAAAAAAAAAAA;
                    int partial_bits = len%64;
                    if(partial_bits > 0)                          // remove trailing Ts
                        destination[destination_size-1] &= (1ULL << partial_bits) - 1;
                } else {
                    int shift_to_right_end = 64*destination_size-len;
                    m_readholderp->CopyBits(position, position+len, destination, shift_to_right_end, destination_size);  // reversed and shifted to the end of the destination
                    for(size_t p = 0; p < destination_size/2; ++p) {
                        swap(destination[p], destination[destination_size-1-p]);
                        Reverse(destination[p]);
                        Reverse(destination[destination_size-1-p]);
                    }
                    if(destination_size%2)
                        Reverse(destination[destination_size/2]);
                }
            }
            // returns number of equal nucleotides (2bit) from the beginning
            // could be longer than actual sequence length if sequence is not multiple of 32
            static size_t CommomSeqLen(const uint64_t* seq1p, const uint64_t* seq2p, size_t word_len) {
                auto last = seq1p+word_len;
                auto mism = mismatch(seq1p, last, seq2p);
                size_t extend = 32*(mism.first-seq1p);
                if(mism.first != last)
                    extend += (ffsll(*mism.first ^ *mism.second)-1)/2; // after ^ all matches are 0s; ffs returns 1-based position of the first bit set to 1
            
                return extend;
            }
            string_iterator& operator++() {
                if(m_read == m_readholderp->m_read_length.size())
                    return *this;
                m_position +=  2*m_readholderp->m_read_length[m_read++]; 
                return *this;
            }
            size_t ReadLen() const { return m_readholderp->m_read_length[m_read]; }
            kmer_iterator KmersForRead(int kmer_len) const {
                if(kmer_len <= (int)m_readholderp->m_read_length[m_read])
                    return kmer_iterator(kmer_len, *m_readholderp, m_position, 0, m_read); 
                else
                    return m_readholderp->kend();
            }
            size_t Hash() const { 
                hash<const CReadHolder*> h1;
                hash<size_t> h2;
                return h1(m_readholderp)^h2(m_position); 
            }
            struct SHash { size_t operator()(const string_iterator& is) const { return is.Hash(); } };
            bool HasMate() const { return m_readholderp->m_contains_paired; }
            int PairType() const {
                if(!m_readholderp->m_contains_paired)
                    return eSingle;
                else if(m_read%2) // odd
                    return eSecondMate;
                else              // even
                    return eFirstMate;
            }
            string_iterator GetMate() const { // undefined behavior if not paired container
                if(m_read%2) // odd
                    return string_iterator(*m_readholderp, m_position-2*m_readholderp->m_read_length[m_read-1], m_read-1);
                else         // even
                    return string_iterator(*m_readholderp, m_position+2*m_readholderp->m_read_length[m_read], m_read+1);
            }

            friend bool operator==(const string_iterator& li, const string_iterator& ri) { return li.m_read == ri.m_read && li.m_readholderp == ri.m_readholderp; }
            friend bool operator!=(const string_iterator& li, const string_iterator& ri) { return li.m_read != ri.m_read || li.m_readholderp != ri.m_readholderp; }
            friend class CReadHolder;
        

        private:
            string_iterator(const CReadHolder& rholder, size_t position = 0, size_t read = 0) : m_readholderp(&rholder), m_position(position), m_read(read) {}

            const CReadHolder* m_readholderp;
            size_t m_position;
            size_t m_read;
        };

    private:
        // efficiently copies sequence to destination without converting it to string
        // assumes that destination is extended properly and filled with 0; destination_size - number of 'used' 8-byte words in destination after copy
        template <typename Dest>
        void CopyBits(size_t bit_from, size_t bit_to, Dest& destination, size_t destination_bit_from, size_t destination_size) const {
            if(bit_to <= bit_from)
                return;

            size_t word = bit_from/64;
            size_t last_word = (bit_to-1)/64;
            unsigned shift = bit_from%64;
            size_t destination_word = destination_bit_from/64; 
            unsigned destination_shift = destination_bit_from%64;
            if(shift > 0) {                                                               // first word partial
                uint64_t chunk = (m_storage[word++] >> shift);
                if(destination_shift > 0) {                                               // destination word partial
                    destination[destination_word] += (chunk << destination_shift);
                    if(shift <= destination_shift)                                        // we used all remaining destination word
                        ++destination_word;
                    if(shift < destination_shift && destination_word < destination_size)  // first word spills out
                        destination[destination_word] += (chunk >> (64-destination_shift));
                } else {                                                                  // desination word is not partial - it is bigger than chunk
                    destination[destination_word] = chunk;
                }
                destination_shift = (destination_shift+64-shift)%64;
            }
            for( ; word <= last_word; ++word, ++destination_word) {
                if(destination_shift > 0) {
                    destination[destination_word] += (m_storage[word] << destination_shift);
                    if(destination_word+1 < destination_size)
                        destination[destination_word+1] += (m_storage[word] >> (64-destination_shift));
                } else {
                    destination[destination_word] = m_storage[word];
                }
            }
            int partial_bits = (destination_bit_from+bit_to-bit_from)%64;
            if(partial_bits > 0) {
                uint64_t mask = (1ULL << partial_bits) - 1;
                destination[destination_size-1] &= mask;
            }
        }


        vector<uint64_t> m_storage;
        vector<uint32_t> m_read_length;
        size_t m_total_seq;
        bool m_contains_paired;
    };

    typedef vector<pair<int,size_t>> TBins; // pair of position,count

    // simple heuristic to find a valley/peak in a histogram
    int FindValleyAndPeak(const TBins& bins, int rlimit) {
        int SLOPE_LEN = 5;
        int peak = min(rlimit,(int)bins.size()-SLOPE_LEN-1);
        while(peak >= SLOPE_LEN) {
            bool maxim = true;
            for(int i = 1; i <= SLOPE_LEN && maxim; ++i) 
                maxim = bins[peak+i].second < bins[peak].second;
            for(int i = 1; i <= SLOPE_LEN && maxim; ++i) 
                maxim = bins[peak-i].second < bins[peak].second;
            if(maxim)
                break;
            --peak;
        }

        if(peak < SLOPE_LEN)
            return -1;

        int valley = 0;
        for(int i = 1; i <= peak; ++i) {
            if(bins[i].second < bins[valley].second)
                valley = i;
        }
        if(valley == peak)
            return -1;

        for(int i = valley; i < (int)bins.size(); ++i) {
            if(bins[i].second > bins[peak].second)
                peak = i;
        }
        
        if(bins[valley].second < 0.7*bins[peak].second)
            return valley;
        else            
            return -1;
    }

    // a simple heuristic to find main range in a histogram
    pair<int,int> HistogramRange(const TBins& bins) {  // returns <valley,rlimit>; valley == -1 if not found
        unsigned MIN_NUM = 100;
        size_t gsize = 0;
        for(auto& bin : bins) {
            if(bin.second >= MIN_NUM) 
                gsize += bin.first*bin.second;
        }

        // step back over repeats and plasmids that are not likely to be more than 20 percent of the genome
        int rl = 0;
        size_t gs = 0;
        for(auto& bin : bins) {
            gs += bin.first*bin.second;
            if(rl < (int)bins.size()-1)
                ++rl;
            if(gs > 0.8*gsize)
                break;
        }

        // find histogram portion with biggest volume and estimate genome size as number of kmers in the portion
        int valley = -1;
        int rlimit = rl;
        size_t genome = 0;
        size_t genome_vol = 0;

        while(true) {
            int v = FindValleyAndPeak(bins, rl);

            size_t g = 0;
            size_t g_vol = 0;
            for(int i = max(0, v); i <= rl; ++i)
            {
                g_vol += (bins[i].first*bins[i].second);
                g += bins[i].second;
            }

            if((v >= 0 && g > genome) || g_vol > genome_vol) {
                valley = v;
                rlimit = rl;
                genome = g;
                genome_vol = g_vol;
                //                cerr << valley << " " << rlimit << " " << genome << endl;
            }

            if(v < 0)
                break;
            rl = v;
        }
        
        return make_pair(valley, rlimit);
    }

    double GetAverageCount(const TBins& bins) {
        pair<int,int> grange =  HistogramRange(bins);
        if(grange.first < 0)
            grange.first = 0;

        size_t genome = 0;
        size_t kmers = 0;
        for(int i = grange.first; i <= grange.second; ++i) {
            genome += bins[i].second;
            kmers += bins[i].first*bins[i].second;
        }

        if(genome > 0)
            return double(kmers)/genome;
        else
            return 0.;
    }

    size_t CalculateGenomeSize(const TBins& bins) {
        pair<int,int> grange =  HistogramRange(bins);
        if(grange.first < 0)
            grange.first = 0;
        size_t genome = 0;
        for(int i = grange.first; i <= grange.second; ++i)
            genome += bins[i].second;            
        
        return genome;
    }
    
    template <typename V> class CKmerMap {
    // A hash with kmer as a key
    // Implemented using a boost::variant of unordered_map<<LargeInt<N>,V> with maximal N = 16 which allows kmer size up to 512

    public:
        typedef V MappedType;
        typedef TKmerMapN<V> Type;

        CKmerMap(int kmer_len = 0) : m_kmer_len(kmer_len) {
            if(m_kmer_len > 0)
                m_container = CreateVariant<TKmerMapN<V>, TLargeIntMap, V>((m_kmer_len+31)/32);
        }
        size_t Size() const { return apply_visitor(container_size(), m_container); }  // number of elements in the container
        void Reserve(size_t rsrv) { apply_visitor(reserve(rsrv), m_container); }      // reserves hash table for rsrv elements
        void Clear() { apply_visitor(clear(), m_container); }                         // clear hash table 
        V& operator[] (const TKmer& kmer) { 
            if(m_kmer_len == 0) 
                throw runtime_error("Can't insert in uninitialized container");
            return apply_visitor(mapper(kmer), m_container);
        }
        V* Find(const TKmer& kmer) { return apply_visitor(find(kmer), m_container); } // returns nullptr if not found
        int KmerLen() const { return m_kmer_len; }

        template <typename Prob> 
        void GetInfo(Prob& prob) { apply_visitor(get_info<Prob>(prob), m_container); } // scans the containier and calls prob(k, v) for each mapped element

    private:

        template <typename Prob>
        struct get_info : public boost::static_visitor<> {
            get_info(Prob& p) : prob(p) {}
            template <typename T> void operator()(T& v) const {
                for(auto& val : v)
                    prob(TKmer(val.first), val.second);
            }

            Prob& prob;
        };
        struct container_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.size();} };
        struct clear : public boost::static_visitor<> { template <typename T> void operator()(const T& v) const { v.clear();} };
        struct reserve : public boost::static_visitor<> { 
            reserve(size_t r) : rsrv(r) {}
            template <typename T> void operator() (T& v) const { v.reserve(rsrv); }
            size_t rsrv;
        };
        struct mapper : public boost::static_visitor<V&> { 
            mapper(const TKmer& k) : kmer(k) {}
            template <typename T> V& operator()(T& v) const { 
                typedef typename  T::key_type large_t;
                return v[kmer.get<large_t>()];
            } 
            const TKmer& kmer;
        };
        struct find : public boost::static_visitor<V*> {
            find(const TKmer& k) : kmer(k) {}
            template <typename T> V* operator()(T& v) const { 
                typedef typename  T::key_type large_t;
                typename T::iterator it = v.find(kmer.get<large_t>());
                if(it != v.end())
                    return &(it->second);
                else
                    return 0;
            } 
            const TKmer& kmer;
        };

        Type m_container;
        int m_kmer_len;
    };
    template <typename V>
    using  TKmerMap = CKmerMap<V>; // for compatibility with previous code

} // namespace
#endif /* _common_util_ */

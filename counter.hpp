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

#ifndef _KmerCounter_
#define _KmerCounter_

#include "Integer.hpp"
#include "common_util.hpp"

namespace DeBruijn {

    class CKmerCount {
    // Class for kmer counting and searching implemented using a boost::variant of vector<pair<LargeInt<N>,size_t>>
    // Currently, maximum N defined in config.hpp is 16 that allows kmers of length at most 512 to be stored.
    // Only smaller (in the bit encoding) of a kmer and its reverse complement is stored
    //
    // When vector is sorted, binary search on the first element of the pair that represents a kmer can be used for retrieving
    // the information stored for the kmer in the second element of the pair.
    // First 32 bits of the second element stores total count for kmer (self and reverse complement)
    // Remaining 32 bits store count for kmer for self only during the counting operations but are modified to additionally store
    // branching information when used inside CDBGraph

    public:
        typedef TKmerCountN Type;

        CKmerCount(int kmer_len = 0) : m_kmer_len(kmer_len) {
            if(m_kmer_len > 0)
                m_container = CreateVariant<TKmerCountN, TLargeIntVec>((m_kmer_len+31)/32);
        }
        size_t Size() const { return apply_visitor(container_size(), m_container); }         // number of elements in the container
        void Reserve(size_t rsrv) { apply_visitor(reserve(rsrv), m_container); }             // reserves memory for rsrv elements
        void Clear() { apply_visitor(clear(), m_container); }                                // clears container (doesn't release memory)
        size_t Capacity() const { return apply_visitor(container_capacity(), m_container); } // tells how many elements could be stored in reserved memory
        size_t ElementSize() const { return apply_visitor(element_size(), m_container); }    // size of one vector element in bytes
        size_t MemoryFootprint() const { return Capacity()*ElementSize(); }                  // reserved memory in bytes
        void PushBack(const TKmer& kmer, size_t count) {                                     // push back one element
            if(m_kmer_len == 0)
                throw runtime_error("Can't insert in uninitialized container");
            apply_visitor(push_back(kmer, count), m_container); 
        }
        void PushBackElementsFrom(const CKmerCount& other) {         // push back elements from other container
            if(m_kmer_len == 0)
                throw runtime_error("Can't insert in uninitialized container");
            apply_visitor(push_back_elements(), m_container, other.m_container); 
        }
        size_t Find(const TKmer& kmer) const { return apply_visitor(find_kmer(kmer), m_container); } // finds index for a kmer (returns Size() if not found)
        void UpdateCount(size_t count, size_t index) { apply_visitor(update_count(count, index), m_container); }          // updates count at the index position
        size_t GetCount(size_t index) const { return apply_visitor(get_count(index), m_container); }                      // gets count at the index position
        pair<TKmer,size_t> GetKmerCount(size_t index) const { return apply_visitor(get_kmer_count(index), m_container); } // gets kmer and count at the index position
        const uint64_t* getPointer(size_t index) { return apply_visitor(get_pointer(index), m_container); }               // gets access to binary kmer sequence
        int KmerLen() const { return m_kmer_len; }
        void Sort() { apply_visitor(container_sort(), m_container); }
        void SortAndExtractUniq(int min_count, CKmerCount& uniq) {  // sorts container, aggregates counts, copies elements with count >= min_count into uniq
            uniq = CKmerCount(m_kmer_len); // init
            Sort();
            apply_visitor(extract_uniq(min_count), m_container, uniq.m_container);
        }
        void SortAndUniq(int min_count) { // sorts container, aggregate counts, keeps elements with count >= min_count
            Sort();
            apply_visitor(uniq(min_count), m_container);
        }
        void RemoveLowCountKmers(int min_count) { apply_visitor(remove_low_count(min_count), m_container); }
        void MergeTwoSorted(const CKmerCount& other) { // merges with other assuming both sorted
            if(m_kmer_len != other.KmerLen())
                throw runtime_error("Can't merge kmers of different lengths");
            apply_visitor(merge_sorted(), m_container, other.m_container);
        }
        void Swap(CKmerCount& other) { // swaps with other
            swap(m_kmer_len, other.m_kmer_len);
            apply_visitor(swap_with_other(), m_container, other.m_container);    
        }
        void Save(ostream& out) const { 
            out.write(reinterpret_cast<const char*>(&m_kmer_len), sizeof(m_kmer_len));
            apply_visitor(save(out), m_container);
            if(!out)
                throw runtime_error("Error in counter write"); 
        }
        void Load(istream& in) {
            if(!in.read(reinterpret_cast<char*>(&m_kmer_len), sizeof(m_kmer_len)))
                throw runtime_error("Error in counter read");
            m_container = CreateVariant<TKmerCountN, TLargeIntVec>((m_kmer_len+31)/32);
            apply_visitor(load(in), m_container);
        }

    private:

        struct find_kmer : public boost::static_visitor<size_t> { 
            find_kmer(const TKmer& k) : kmer(k) {}
            template <typename T> size_t operator()(const T& v) const { 
                typedef typename T::value_type pair_t;
                typedef typename pair_t::first_type large_t;
                auto it = lower_bound(v.begin(), v.end(), kmer.get<large_t>(), [](const pair_t& element, const large_t& target){ return element.first < target; });
                if(it == v.end() || it->first != kmer.get<large_t>())
                    return v.size();
                else
                    return it-v.begin();
            } 
            const TKmer& kmer;
        };
        struct reserve : public boost::static_visitor<> { 
            reserve(size_t r) : rsrv(r) {}
            template <typename T> void operator() (T& v) const { v.reserve(rsrv); }
            size_t rsrv;
        };
        struct container_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.size();} };
        struct container_capacity : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.capacity();} };
        struct element_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return sizeof(typename  T::value_type);} };
        struct clear : public boost::static_visitor<> { template <typename T> void operator()(T& v) const { v.clear();} };
        struct push_back : public boost::static_visitor<> { 
            push_back(const TKmer& k, size_t c) : kmer(k), count(c) {}
            template <typename T> void operator() (T& v) const {
                typedef typename  T::value_type::first_type large_t;
                v.push_back(make_pair(kmer.get<large_t>(), count)); 
            }
            const TKmer& kmer;
            size_t count;
        };
        struct push_back_elements : public boost::static_visitor<> {
            template <typename T> void operator() (T& a, const T& b) const { a.insert(a.end(), b.begin(), b.end()); }
            template <typename T, typename U> void operator() (T& a, const U& b) const { throw runtime_error("Can't copy from different type container"); }
        };
        struct merge_sorted : public boost::static_visitor<> {
            template <typename T> void operator() (T& a, const T& b) const { 
                T merged;
                merged.reserve(a.size()+b.size());
                merge(a.begin(), a.end(), b.begin(), b.end(), back_inserter(merged));
                merged.swap(a); 
            }
            template <typename T, typename U> void operator() (T& a, const U& b) const { throw runtime_error("Can't merge different type containers"); }
        };
        struct update_count : public boost::static_visitor<> {
            update_count(size_t c, size_t i) : count(c), index(i) {}
            template <typename T> void operator() (T& v) const { v[index].second = count; }  
            size_t count;
            size_t index;
        };
        struct get_count : public boost::static_visitor<size_t> {
            get_count(size_t i) : index(i) {}
            template <typename T> size_t operator() (T& v) const { return v[index].second; }        
            size_t index;
        };
        struct get_kmer_count : public boost::static_visitor<pair<TKmer,size_t>> {
            get_kmer_count(size_t i) : index(i) {}
            template <typename T> pair<TKmer,size_t> operator() (T& v) const { return make_pair(TKmer(v[index].first), v[index].second); }        
            size_t index;
        };
        struct get_pointer : public boost::static_visitor<const uint64_t*> {
            get_pointer(size_t i) : index(i) {}
            template <typename T> const uint64_t* operator() (T& v) const { return v[index].first.getPointer(); }        
            size_t index;
        };
        struct container_sort : public boost::static_visitor<> { template <typename T> void operator() (T& v) const { sort(v.begin(), v.end()); }};
        struct swap_with_other : public boost::static_visitor<> { 
            template <typename T> void operator() (T& a, T& b) const { a.swap(b); }
            template <typename T, typename U> void operator() (T& a, U& b)  const { throw runtime_error("Can't swap different type containers"); }
        };
        
        struct remove_low_count : public boost::static_visitor<> {
            remove_low_count(int mc) : min_count(mc) {}
            template <typename T> void operator() (T& v) const {
                v.erase(remove_if(v.begin(), v.end(), [this](const typename T::value_type& pair) { return (uint32_t)pair.second < this->min_count; }), v.end());
            }

            unsigned min_count;
        };

        struct uniq : public boost::static_visitor<> {
            uniq(int mc) : min_count(mc) {}
            template <typename T> void operator() (T& v) const {
                typedef typename T::iterator iter_t;
                iter_t nextp = v.begin();
                for(iter_t ip = v.begin(); ip != v.end(); ) {
                    iter_t workp = ip;
                    while(++ip != v.end() && workp->first == ip->first)
                        workp->second += ip->second;   // accumulate all 8 bytes; we assume that count will not spill into higher half
                    if((uint32_t)workp->second >= min_count)
                        *nextp++ = *workp;
                }
                v.erase(nextp, v.end());
            }

            unsigned min_count;
        };
        struct extract_uniq : public boost::static_visitor<> {
            extract_uniq(int mc) : min_count(mc) {}
            template <typename T> void operator() (T& a, T& b) const {
                if(a.empty()) return;
                size_t num = 1;
                uint32_t count = a[0].second;  // count only 4 bytes!!!!!!
                for(size_t i = 1; i < a.size(); ++i) {
                    if(a[i-1].first < a[i].first) {
                        if(count >= min_count)
                            ++num;
                        count = a[i].second;
                    } else {
                        count += a[i].second;
                    }
                }
                if(count < min_count)
                    --num;
                b.reserve(num+1);
                b.push_back(a[0]);
                for(size_t i = 1; i < a.size(); ++i) {
                    if(b.back().first < a[i].first) {
                        if((uint32_t)b.back().second < min_count)
                            b.pop_back();
                        b.push_back(a[i]);
                    } else {
                        b.back().second += a[i].second;  // accumulate all 8 bytes; we assume that count will not spill into higher half
                    }
                }
                if((uint32_t)b.back().second < min_count)
                    b.pop_back();
            }
            template <typename T, typename U> void operator() (T& a, U& b) const { throw runtime_error("Can't extract into different type container"); }

            unsigned min_count;
        };
        struct save : public boost::static_visitor<> {
            save(ostream& out) : os(out) {}
            template <typename T> void operator() (T& v) const {
                size_t num = v.size();
                os.write(reinterpret_cast<const char*>(&num), sizeof num); 
                if(num > 0)
                    os.write(reinterpret_cast<const char*>(&v[0]), num*sizeof(v[0]));
            }
            ostream& os;
        };
        struct load : public boost::static_visitor<> {
            load(istream& in) : is(in) {}
            template <typename T> void operator() (T& v) const {
                size_t num;
                if(!is.read(reinterpret_cast<char*>(&num), sizeof num))
                    throw runtime_error("Error in counter read");
                if(num > 0) {
                    v.resize(num);
                    if(!is.read(reinterpret_cast<char*>(&v[0]), num*sizeof(v[0])))
                        throw runtime_error("Error in counter read");
                }
            }
            istream& is;
        };

        Type m_container;
        int m_kmer_len;
    };
    typedef CKmerCount TKmerCount; // for compatibility with previous code


    // CKmerCounter counts kmers in reads using multiple threads and stores them in TKmerCount
    // It also finds neighbors (in GetBranches) if a user wants to use this class to build a CDBGraph (de Bruijn graph)
    // As Kmer counting could be memory expensive, CKmerCounter accepts an upper limit for the memory available and will 
    // subdivide the task, if needed.
    // If the number of subtasks exceeds 10, it will throw an exception asking for more memory.

    class CKmerCounter {
    public:

        // reads - raw reads (ncores or more elements in the list)
        // kmer_len - size of kmer
        // min_count - minimal count for accepted kmers
        // is_stranded - flag indicating whether kmers are from input reads where strand is informative or from connected paired
        //               reads generated internally by the program where strand is not a meaningful observation
        // mem_available - allowed memory in bytes
        // ncores - number of cores
        CKmerCounter(const list<array<CReadHolder,2>>& reads, int kmer_len, int min_count, bool is_stranded, int64_t mem_available, int ncores) : 
            m_kmer_len(kmer_len), m_min_count(min_count), m_is_stranded(is_stranded), m_ncores(ncores), m_reads(reads) {

            cerr << endl << "Kmer len: " << m_kmer_len << endl;
            CStopWatch timer;
            timer.Restart();

            int64_t raw_kmer_num = 0;
            for(const auto& reads : m_reads)
                raw_kmer_num += reads[0].KmerNum(m_kmer_len)+reads[1].KmerNum(m_kmer_len);

            int64_t GB = 1000000000;
            int kmer_size = TKmerCount(m_kmer_len).ElementSize();
            int64_t mem_needed = 1.2*raw_kmer_num*kmer_size;

            int max_cycles = 10;  // maximum cycles allowed
            int64_t mbuf = 2*GB;  // memory buffer for allocation uncertainity
            if(mem_needed >= max_cycles*(mem_available-mbuf)) {
                throw runtime_error("Memory provided is insufficient to do runs in 10 cycles for the read coverage. We find that 16 Gb for 20x coverage of a 5 Mb genome is usually sufficient");
            }
            int cycles = ceil(double(mem_needed)/(mem_available-mbuf));

            cerr << "Raw kmers: " << raw_kmer_num  << " Memory needed (GB): " << double(mem_needed)/GB << " Memory available (GB): " << double(mem_available-mbuf)/GB << " " << cycles << " cycle(s) will be performed" << endl;
        
            int njobs = 8*m_reads.size();   // many buckets reduce short-lived memory overhead spike in SortAndMergeJob    
            int kmer_buckets = cycles*njobs; 
    
            for(int cycl = 0; cycl < cycles; ++cycl) {
                pair<int,int> bucket_range(cycl*njobs, (cycl+1)*njobs-1);
                list<vector<TKmerCount>> raw_kmers;

                list<function<void()>> jobs;
                for(auto& job_input : m_reads) {
                    if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty       
                        raw_kmers.push_back(vector<TKmerCount>());
                        jobs.push_back(bind(&CKmerCounter::SpawnKmersJob, this, ref(job_input), kmer_buckets, bucket_range, ref(raw_kmers.back())));
                    }
                }
                RunThreads(ncores, jobs);

                // size_t total = 0;
                // for(auto& v : raw_kmers) {
                //     for(auto& tc : v) 
                //         total += tc.MemoryFootprint();
                // }
                SortAndMergeKmers(raw_kmers);
            }
    
            size_t utotal = 0;
            for(auto& c : m_uniq_kmers)
                utotal += c.Size();

            cerr << "Distinct kmers: " << utotal << endl;    
            cerr << "Kmer count in " << timer.Elapsed();
            
            MergeSortedKmers();
            if(m_uniq_kmers.empty())
                m_uniq_kmers.push_back(TKmerCount(m_kmer_len));                        
        }
        virtual ~CKmerCounter() {}

        // reference to counted kmers
        TKmerCount& Kmers() { return m_uniq_kmers.front(); }
        const TKmerCount& Kmers() const { return m_uniq_kmers.front(); }

        // average count of kmers in the histogram with the main peak
        double AverageCount() const {
            map<int,size_t> bins;
            for(size_t index = 0; index < Kmers().Size(); ++index) {
                ++bins[Kmers().GetCount(index)];                  // count clipped to integer automatically
            }
            TBins hist(bins.begin(), bins.end());
            return GetAverageCount(hist);
        }

        // prepares kmer counts to be used in  CDBGraph (de Bruijn graph)
        // runs multiple instances of GetBranchesJob
        void GetBranches() {
            CStopWatch timer;
            timer.Restart();
            if(Kmers().Size() > 0) {
                vector<uint8_t> branches(Kmers().Size());
                size_t bucket_size = Kmers().Size()/m_ncores+1;
                list<function<void()>> jobs;
                for(int i = 0; i < m_ncores; ++i) {
                    pair<size_t,size_t> range(bucket_size*i,min(bucket_size*(i+1)-1,Kmers().Size()-1));
                    if(range.second >= range.first)
                        jobs.push_back(bind(&CKmerCounter::GetBranchesJob, this, range, ref(branches)));
                }
                RunThreads(m_ncores, jobs);

                for(size_t index = 0; index < Kmers().Size(); ++index) {
                    size_t b = branches[index];
                    size_t count = Kmers().GetCount(index);
                    uint32_t total_count = count;
                    uint32_t plus_count = (count >> 32);
                    size_t plusf = uint16_t(double(plus_count)/total_count*numeric_limits<uint16_t>::max()+0.5);
                    Kmers().UpdateCount((plusf << 48)+(b << 32)+total_count, index);  // we put strand info and branching in the high half of the count!!!!!                   
                }
            }

            cerr << "Kmers branching in " << timer.Elapsed();
        }

        bool IsStranded() const { return m_is_stranded; }              // indicates if contains stranded information

    private:

        // one-thread worker producing kmers and putting them in multiple non-overlapping buckets
        // rholder - input reads 
        // buckets - total number of buckets
        // bucket_range - range of buckets used by this worker
        // kmers - output kmers
        void SpawnKmersJob(const array<CReadHolder,2>& rholder, int buckets, pair<int,int> bucket_range,  vector<TKmerCount>& kmers) {
            size_t total = rholder[0].KmerNum(m_kmer_len)+rholder[1].KmerNum(m_kmer_len);
            size_t reserve = 1.1*total/buckets;
            int active_buckets = bucket_range.second-bucket_range.first+1;
            kmers.resize(active_buckets, TKmerCount(m_kmer_len));
            for(auto& k : kmers)
                k.Reserve(reserve);

            for(int p = 0; p < 2; ++p) {
                for(CReadHolder::kmer_iterator itk = rholder[p].kbegin(m_kmer_len); itk != rholder[p].kend(); ++itk) {
                    TKmer kmer = *itk;
                    TKmer rkmer = revcomp(kmer, m_kmer_len);
                    size_t count = 1;
                    TKmer* min_kmerp = &rkmer;
                    if(kmer < rkmer) {
                        min_kmerp = &kmer;
                        count += (size_t(1) << 32);
                    }
                    int bucket = min_kmerp->oahash()%buckets;
                    if(bucket < bucket_range.first || bucket > bucket_range.second)
                        continue;
                    // good to go   
                    int ind = bucket - bucket_range.first;
                    if(kmers[ind].Size() == kmers[ind].Capacity()) { //expensive plan B for the case of failed hash uniformity          
                        //            cerr << "Warning: Hash uniformity problem" << endl;           
                        TKmerCount bigger(m_kmer_len);
                        bigger.Reserve(kmers[ind].Size()*1.2);
                        bigger.PushBackElementsFrom(kmers[ind]);
                        bigger.Swap(kmers[ind]);
                    }
                    kmers[ind].PushBack(*min_kmerp, count);            
                }
            }
        }

        //SortAndMergeJob briefly doubles the input memory - should be executed in small chunks!!!!!!   
        // one-thread worker which accepts all containers for a given bucket and merges, sorts and counts them
        // group - list of containers
        // ukmers - counted kmers
        typedef list<TKmerCount*> TContainerPList;
        void SortAndMergeJob(TContainerPList group, TKmerCount& ukmers) {
            TKmerCount all_kmers(group.front()->KmerLen());
            if(group.size() == 1) {
                all_kmers.Swap(*group.front());
            } else {
                size_t total = 0;
                for(auto p : group)
                    total += p->Size();
                   
                all_kmers.Reserve(total); // doubles the input memory!!!!       
            
                for(auto p : group) {
                    all_kmers.PushBackElementsFrom(*p);
                    TKmerCount(p->KmerLen()).Swap(*p);
                }
            }

            all_kmers.SortAndExtractUniq(m_min_count, ukmers);         
        }

        // runs multiple instances of SortAndMergeJob and stores results in m_uniq_kmers
        // raw_kmers - input kmers
        void SortAndMergeKmers(list<vector<TKmerCount>>& raw_kmers) {

            list<function<void()>> jobs;
            int bucken_num = raw_kmers.front().size();

            for(int bucket = 0; bucket < bucken_num; ++bucket) {
                TContainerPList job_input;
                for(auto& vec : raw_kmers)
                    job_input.push_back(&vec[bucket]);
                m_uniq_kmers.push_back(TKmerCount());
                jobs.push_back(bind(&CKmerCounter::SortAndMergeJob, this, job_input, ref(m_uniq_kmers.back())));
            }
            RunThreads(m_ncores, jobs);
        }

        // one-thread worker which merges two sorted buckets
        static void MergeSortedJob(TKmerCount& akmers, TKmerCount& bkmers) {
            akmers.MergeTwoSorted(bkmers);
            TKmerCount(bkmers.KmerLen()).Swap(bkmers); // release bkmers memory
        }

        // runs multiple instances of MergeSortedJob
        // at the end m_uniq_kmers has only one element with final kmers
        void MergeSortedKmers() {
            CStopWatch timer;
            timer.Restart();
            while(m_uniq_kmers.size() > 1) {
                list<function<void()>> jobs;
                for(list<TKmerCount>::iterator first = m_uniq_kmers.begin(); first != m_uniq_kmers.end(); ++first) {
                    list<TKmerCount>::iterator second = first;
                    if(++second != m_uniq_kmers.end()) {
                        jobs.push_back(bind(MergeSortedJob, ref(*first), ref(*second)));
                        first = second;
                    }
                }
                RunThreads(m_ncores, jobs);
                for(auto iloop = m_uniq_kmers.begin(); iloop != m_uniq_kmers.end(); ) {
                    auto it = iloop++;
                    if(it->Size() == 0)
                        m_uniq_kmers.erase(it);
                }
            }
            cerr << "Uniq kmers merging in " << timer.Elapsed();
        }

        // one-thread worker which calculates the branching information (neighbors) for a range of kmers
        // range - from,to indexes for kmers
        // branches - vector of branching information (one bit is used for each of the eight possible neighbors)  
        void GetBranchesJob(pair<size_t,size_t> range, vector<uint8_t>& branches) {
            TKmer max_kmer(string(m_kmer_len, bin2NT[3]));
            for(size_t index = range.first; index <= range.second; ++index) {
                pair<TKmer,size_t> kmer_count = Kmers().GetKmerCount(index);
                //direct        
                TKmer shifted_kmer = (kmer_count.first << 2) & max_kmer;
                //inverse       
                TKmer shifted_rkmer = (revcomp(kmer_count.first, m_kmer_len) << 2) & max_kmer;
                for(int nt = 0; nt < 4; ++nt) {
                    TKmer k = shifted_kmer + TKmer(m_kmer_len, nt);
                    size_t new_index = Kmers().Find(min(k, revcomp(k, m_kmer_len)));
                    // New kmer is a neighbor if it exists in reads and is not same as current kmer
                    if(new_index != Kmers().Size() && new_index != index) 
                        branches[index] |= (1 << nt);

                    k = shifted_rkmer + TKmer(m_kmer_len, nt);
                    new_index = Kmers().Find(min(k, revcomp(k, m_kmer_len)));
                    if(new_index != Kmers().Size() && new_index != index) 
                        branches[index] |= (1 << (nt+4));
                }
            }
        }

        int m_kmer_len;
        int m_min_count;
        bool m_is_stranded;
        int m_ncores;
        const list<array<CReadHolder,2>>& m_reads;
        list<TKmerCount> m_uniq_kmers;                       // storage for kmer buckets; at the end will have one element which is the result     
    };

}; // namespace
#endif /* _KmerCounter_ */

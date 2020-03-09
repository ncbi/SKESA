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

#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/seek.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "gfa.hpp"
#include "readsgetter.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

namespace DeBruijn {

    void ConnectContigsJob(map<string, string>& contigs, map<string,Node>& lkmers, map<string,Node>& rkmers, 
                           DBGraph& graph, double fraction, int ext_len, 
                           vector<pair<string, SAtomic<uint8_t>>>& sentinels, TSpiderCollection& spiders) {
        array<string, 2> cends = {"3p", "5p"};
        for(auto& sentinel : sentinels) {
            if(!sentinel.second.Set(1, 0))
                continue;
            const string& acc = sentinel.first;
            for(int e = 0; e < 2; ++e) {
                Node node = (e == 0) ? rkmers[acc] : lkmers[acc]; 
                spiders.emplace_back(contigs, graph, fraction, acc+":"+cends[e]);
                spiders.back().ConnectOneEnd(node, ext_len);
                if(spiders.back().empty())
                    spiders.pop_back();
                else
                    spiders.back().DetectCycles();
            }        
        } 
    }

}; // namespace


int main(int argc, const char* argv[])
{

    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]")
        ("estimated_kmers", value<int>()->default_value(100), "Estimated number of unique kmers for bloom filter (millions) [integer]")
        ("skip_bloom_filter", "Don't do bloom filter; use --estimated_kmers as the hash table size [flag]")
        ;

    options_description input("Input/output options");
    input.add_options()
#ifndef NO_NGS
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
#endif
        ("reads", value<vector<string>>(), "Input fasta/fastq file(s) (could be used multiple times for different runs, could be gzipped) [string]")
        ("use_paired_ends", "Indicates that a single (not comma separated) fasta/fastq file contains paired reads [flag]")
        ("contigs", value<string>(), "Input file with contigs [string]")
        ("gfa", value<string>(), "GFA graph output (stdout if not specified) [string]")
        ("dbg", value<string>(), "Input de Bruijn graph (optional) [string]")
        ("csv", value<string>(), "CSV file ouput (optional) [string]")
        ("contig_color", value<string>()->default_value("Purple"), "Color for contigs [string]")
        ("connector_color", value<string>()->default_value("Lime Green"), "Color for connectors [string]")
        ;

    options_description assembly("Assembly options");
    assembly.add_options()
        ("kmer", value<int>()->default_value(41), "Kmer length for assembly [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("vector_percent", value<double>()->default_value(0.05, "0.05"), "Count for vectors as a fraction of the read number (1. disables) [float (0,1]]")
        ("fraction", value<double>()->default_value(0.1, "0.1"), "Threshold for extension")
        ("entropy", value<double>()->default_value(0.51, "0.51"), "Minimal entropy for a seed kmer [float]")
        ("ext_len", value<int>()->default_value(2000), "Maximal length for extension")
        ;

    options_description filter("Graph cleaning options");
    filter.add_options()
        ("not_aligned_len", value<int>()->default_value(10), "Not aligned read length for break count")
        ("not_aligned_count", value<int>()->default_value(3), "Number of not aligned reads to make a break")
        ("aligned_count", value<int>()->default_value(2), "Number of aligned reads to confirm a connection")
        ("max_path", value<int>()->default_value(1000), "Maximal number of paths allowed in 1 step of filtering")
        ("no_filter_by_reads", "Don't use full length reads for variants filtering [flag]")
        ("no_filter_by_pairs", "Don't use mate pairs for variants filtering [flag]")
        ;

    options_description all("");
    all.add(general).add(input).add(assembly).add(filter);

    try {
        variables_map argm;                                // boost arguments
        store(parse_command_line(argc, argv, all), argm);
        notify(argm);    

        if(argm.count("help")) {
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cout << all << "\n";
            return 0;
        }

        if(argm.count("version")) {
            cout << "gfa_connector 1.1.0" << endl;
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl << endl;
#endif
            return 0;
        }

        ofstream gfa_out;
        if(argm.count("gfa")) {
            gfa_out.open(argm["gfa"].as<string>());
            if(!gfa_out.is_open()) {
                cerr << "Can't open file " << argm["gfa_out"].as<string>() << endl;
                exit(1);
            }
        }

        ofstream csv_out;
        if(argm.count("csv")) {
            csv_out.open(argm["csv"].as<string>());
            if(!csv_out.is_open()) {
                cerr << "Can't open file " << argm["csv_out"].as<string>() << endl;
                exit(1);
            }
        }

        int not_aligned_len = argm["not_aligned_len"].as<int>(); 
        int not_aligned_count = argm["not_aligned_count"].as<int>();
        int aligned_count = argm["aligned_count"].as<int>();
        int maxp = argm["max_path"].as<int>();
        bool no_reads = argm.count("no_filter_by_reads");
        bool no_pairs = argm.count("no_filter_by_pairs");
        bool need_reads = !argm.count("dbg") || !no_reads || !no_pairs;
        if(need_reads && !argm.count("reads") 
#ifndef NO_NGS
           && !argm.count("sra_run")
#endif
                                                                   ) {
            cerr << "Provide some input reads" << endl;
            cerr << all << "\n";
            return 1;
        }

        vector<string> sra_list;
        vector<string> reads_list;
#ifndef NO_NGS
        if(argm.count("sra_run")) {
            sra_list = argm["sra_run"].as<vector<string>>();
            unsigned num = sra_list.size();
            sort(sra_list.begin(), sra_list.end());
            sra_list.erase(unique(sra_list.begin(),sra_list.end()), sra_list.end());
            if(sra_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from SRA run list" << endl; 
        }
#endif
        if(argm.count("reads")) {
            reads_list = argm["reads"].as<vector<string>>();
            unsigned num = reads_list.size();
            sort(reads_list.begin(), reads_list.end());
            reads_list.erase(unique(reads_list.begin(),reads_list.end()), reads_list.end());
            if(reads_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from reads file list" << endl; 
        }

        bool usepairedends = argm.count("use_paired_ends");

        int ncores = thread::hardware_concurrency();
        if(argm["cores"].as<int>()) {
            int nc = argm["cores"].as<int>();
            if(nc < 0) {
                cerr << "Value of --cores must be >= 0" << endl;
                exit(1);
            } else if(nc > ncores) {
                cerr << "WARNING: number of cores was reduced to the hardware limit of " << ncores << " cores" << endl;
            } else if(nc > 0) {
                ncores = nc;
            }
        }

        double fraction = argm["fraction"].as<double>();
        if(fraction >= 1.) {
            cerr << "Value of --fraction must be < 1 (more than 0.25 is not recommended)" << endl;
            exit(1);
        }
        if(fraction < 0.) {
            cerr << "Value of --fraction must be >= 0" << endl;
            exit(1);
        }
        double entropy_level = argm["entropy"].as<double>();
        int min_count = argm["min_count"].as<int>();
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }
        double vector_percent = argm["vector_percent"].as<double>();
        if(vector_percent > 1.) {
            cerr << "Value of --vector_percent  must be <= 1" << endl;
            exit(1);
        }
        if(vector_percent <= 0.) {
            cerr << "Value of --vector_percent  must be > 0" << endl;
            exit(1);
        }

        int estimated_kmer_num =  argm["estimated_kmers"].as<int>();
        if(estimated_kmer_num <= 0) {
            cerr << "Value of --estimated_kmers must be > 0" << endl;
            exit(1);
        }
        bool skip_bloom_filter = argm.count("skip_bloom_filter");

        list<array<CReadHolder,2>> reads;
        if(need_reads) {
            CReadsGetter readsgetter(sra_list, reads_list, ncores, usepairedends);
            if(vector_percent < 1.) {
                readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmer_num, skip_bloom_filter);
                readsgetter.PrintAdapters();
            } else {
                cerr << "Adapters clip is disabled" << endl;
            }
            reads.splice(reads.end(), readsgetter.Reads());
        }

        int kmer_len = argm["kmer"].as<int>();
        if(kmer_len%2 ==0) {
            cerr << "Kmer must be an odd number" << endl;
            return 1;
        }

        unique_ptr<DBGraph> graphp;
        if(argm.count("dbg")) {
            ifstream file(argm["dbg"].as<string>(), ios::binary | ios::in);
            if(!file.is_open()) {
                cerr << "Can't open file " << argm["dbg"].as<string>() << endl;
                return 1;             
            }
            graphp.reset(new DBGraph(file));
            kmer_len = graphp->KmerLen();
            cerr << "Loaded hash graph for kmer: " << kmer_len << endl;
        } else {
            int64_t M = 1000000;
            CKmerHashCounter  kmer_counter(reads, kmer_len, min_count, M*estimated_kmer_num, true, ncores, skip_bloom_filter);
            if(kmer_counter.KmerNum() == 0)
                throw runtime_error("Not enough quality reads which are at least "+to_string(kmer_len)+"bp long");
            kmer_counter.GetBranches();
            graphp.reset(new CDBHashGraph(move(kmer_counter.Kmers()), true));
        }

        map<string, string> contigs;  // acc, seq
        map<string,Node> rkmers; 
        map<string,Node> lkmers;      // reversed
        map<Node,int> end_kmers_count; 
        if(argm.count("contigs")) {
            ifstream file(argm["contigs"].as<string>());
            if(!file.is_open()) {
                cerr << "Can't open file " << argm["contigs"].as<string>() << endl;
                return 1;             
            }
            char c;
            if(!(file >> c) || c != '>')
                throw runtime_error("Invalid fasta file format for contigs");
            string record;
            while(getline(file, record, '>')) {
                while(!file.eof() && record.back() != '\n') {
                    string part;
                    getline(file, part, '>');
                    record += '>'+part;
                }
                size_t first_ret = min(record.size(),record.find('\n'));
                if(first_ret == string::npos)
                    throw runtime_error("Invalid fasta file format for contigs");
                string acc = record.substr(0, first_ret);
                acc = acc.substr(0, acc.find_first_of(" \t"));
                string seq = record.substr(first_ret+1);
                seq.erase(remove(seq.begin(),seq.end(),'\n'),seq.end());
                for(char& c : seq) c = toupper(c);
                if(seq.find_first_not_of("ACGTYRWSKMDVHBXN-") != string::npos)
                    throw runtime_error("Invalid sequence in fasta file for contigs");
                if((int)seq.size() < 2*kmer_len+1)
                    continue;
                //                seq = seq.substr(kmer_len, seq.size()-2*kmer_len);
                int len = seq.size();                
                contigs[acc] = seq;
                {
                    string rkmer = seq.substr(len-kmer_len);
                    Node node = graphp->GetNode(rkmer);
                    if(!node.isValid()) {
                        cerr << "Rkmer: " << rkmer << " not in graph" << endl;                    
                    } else {
                        rkmers[acc] = node;
                        ++end_kmers_count[node];
                    }
                }
                {
                    string lkmer = seq.substr(0, kmer_len);
                    Node node = graphp->GetNode(lkmer).ReverseComplement();
                    if(!node.isValid()) {
                        cerr << "Lkmer: " << lkmer << " not in graph" << endl;
                        continue;
                    } else {
                        lkmers[acc] = node;
                        ++end_kmers_count[node];
                    }
                }
            }
            if(file.bad())
                throw runtime_error("Error in reading contigs");

            cerr << "Contigs: " << contigs.size() << endl;
        } else {
            cerr << "Provide contigs" << endl;
            cerr << all << "\n";
            return 1;
        }

        CStopWatch timer;
        timer.Restart(); 

        int ext_len = argm["ext_len"].as<int>();                
        TSpiderCollection spiders; 

        {
            vector<pair<string, SAtomic<uint8_t>>> sentinels;
            for(auto& contig : contigs)
                sentinels.emplace_back(contig.first, 0);
 
            list<TSpiderCollection> job_rslts;
            list<function<void()>> jobs;            
            for(int thr = 0; thr < ncores; ++thr) {
                job_rslts.emplace_back();
                jobs.push_back(bind(ConnectContigsJob, ref(contigs), ref(lkmers), ref(rkmers), ref(*graphp), fraction, ext_len, ref(sentinels), ref(job_rslts.back())));
            }
            RunThreads(ncores, jobs);

            for(auto& rslt : job_rslts)
                spiders.splice(spiders.end(), rslt);
            SortCollection(spiders);
        }

        cerr << "Assembling in " << timer.Elapsed();
        timer.Restart(); 
                 
        GraphCleaner<TSpiderCollection> cleaner(spiders, *graphp, fraction, entropy_level, not_aligned_len, not_aligned_count, aligned_count, maxp, no_reads, no_pairs, reads, ncores);

        {
            TSpiderCollection cleaned;
            for(auto& spider : spiders) {
                while(!spider.empty()) {
                    spider.UpdateEndKmers();
                    spider.MergeForks();
                    if(!spider.RemoveLooseEnds())
                        break;
                }
                if(!spider.empty()) {
                    spider.AssignGroupNumber();
                    TSpiderCollection sp = spider.SplitGroups();
                    cleaned.splice(cleaned.end(), sp);
                } else {
                    cerr << "ErasedA " << spider.Target() << endl;
                }
            }

            for(auto it_loop = cleaned.begin(); it_loop != cleaned.end(); ) {
                auto it = it_loop++;
                if(it->Connections() < 2) {
                    cerr << "ErasedB " << it->Target() << endl;
                    cleaned.erase(it);
                } else {
                    it->GenerateKmers(*graphp);
                }
            }

            RemoveSpiderSubGraphs(cleaned);            
            swap(spiders, cleaned);
        } 
                
        for(bool keep_doing = true; keep_doing; ) {
            keep_doing = false;
            for(auto it = spiders.begin(); it != spiders.end(); ++it) {
                for(auto jt_loop = next(it); jt_loop != spiders.end(); ) {
                    auto jt = jt_loop++;
                    if(Spider::EndsIntersect(it->EndKmers(), jt->EndKmers())) {  
                        cerr << it->Target() << ":" << it->front().m_group << " " << it->KSignature().size() << " absorbs " << jt->Target() << ":" << jt->front().m_group << " " << jt->KSignature().size() << endl;
                        keep_doing = true;
                        it->Absorb(*jt);
                        it->GenerateKmers(*graphp);
                        spiders.erase(jt);
                    }
                }
            }
        }
        cerr << "Cleaning in " << timer.Elapsed();
        timer.Restart();

        struct SpiderSegmentP {
            SpiderSegmentP(TSpiderCollection::iterator spideri, GFAIterator segmi) : m_spideri(spideri), m_segmi(segmi) {}
            bool operator==(const SpiderSegmentP& other) const { return m_spideri == other.m_spideri && m_segmi == other.m_segmi; }
            struct Hash { size_t operator()(const SpiderSegmentP& p) const { return hash<void*>()(&(*p.m_spideri))^hash<void*>()(&(*p.m_segmi));} };

            TSpiderCollection::iterator m_spideri;
            GFAIterator m_segmi;
        };
        struct SpiderSegmentPDir : public SpiderSegmentP {
            SpiderSegmentPDir (TSpiderCollection::iterator spideri, GFAIterator segmi, bool rend) : SpiderSegmentP(spideri, segmi), m_right_end(rend) {}
            bool operator==(const SpiderSegmentPDir& other) const { return m_spideri == other.m_spideri && m_segmi == other.m_segmi && m_right_end == other.m_right_end; }
            struct Hash { size_t operator()(const SpiderSegmentPDir& p) const { return SpiderSegmentP::Hash()(p)^hash<bool>()(p.m_right_end);} };

            bool m_right_end;
        };

        map<Node, list<SpiderSegmentPDir>> end_segments;                                          //maps kmer to segment with free end
        for(auto spideri = spiders.begin(); spideri != spiders.end(); ++spideri) {
            for(auto segmi = spideri->begin(); segmi != spideri->end(); ++segmi) {
                if(segmi->m_right_connections.empty()) {
                    for(Node& node : segmi->m_seq.back().m_right_kmers)
                        end_segments[node.ReverseComplement()].emplace_back(spideri, segmi, true);
                }
                if(segmi->m_left_connections.empty()) {
                    for(Node& node : segmi->m_seq.front().m_left_kmers)
                        end_segments[node].emplace_back(spideri, segmi, false);
                }
            }
        }        

        unordered_set<SpiderSegmentP, SpiderSegmentP::Hash> erased_segments;                               //segments completely moved to contig ends
        map<string, list<SpiderSegmentPDir>> left_contig_links;                                            //contig ID, links
        map<string, list<SpiderSegmentPDir>> right_contig_links;                                           //contig ID, links
        unordered_map<SpiderSegmentPDir, list<tuple<string,bool>>, SpiderSegmentPDir::Hash> segment_links; //segment to contig ID, 'true' if right contig end

        for(auto& contig : contigs) {
            auto& acc = contig.first;
            auto& seq = contig.second;
            if(rkmers.count(acc)) {
                Node& rnode = rkmers[acc];
                auto rslt = end_segments.find(rnode);
                if(rslt != end_segments.end()) {
                    auto lst = rslt->second;
                    if(lst.size() > 1 || end_kmers_count[rnode] > 1) { // multiple connection - clip kmer from contig and keep segment
                        seq.erase(seq.end()-kmer_len, seq.end());
                    } else {                                           // one contig, one segment - include segment in contig and mark segment for deletion
                        SpiderSegmentPDir& segmp = lst.front();
                        auto& seg_seq = segmp.m_segmi->m_seq;
                        int seg_len = seg_seq.size();
                        if(seg_len <= kmer_len) {
                            seq.erase(seq.end()-kmer_len+seg_len, seq.end());
                        } else {                            
                            string s;
                            if(segmp.m_right_end) {                       // right contig to right segment - different strands
                                for(int i = 0; i < seg_len-kmer_len; ++i)
                                    s.push_back(seg_seq[i].m_nt);
                                ReverseComplementSeq(s.begin(), s.end());
                            } else {                                      // right contig to left segment - same strands
                                for(int i = kmer_len; i < seg_len; ++i)
                                    s.push_back(seg_seq[i].m_nt);
                            }
                            seq += s;
                        }
                        erased_segments.insert(segmp);
                        seg_seq.clear();
                    }
                    for(SpiderSegmentPDir& segmp : lst)
                        segment_links[segmp].emplace_back(acc, true);
                    auto& rcl = right_contig_links[acc];
                    rcl.splice(rcl.end(), lst);
                }
            }
            if(lkmers.count(acc)) {
                Node lnoder = lkmers[acc];
                auto rslt = end_segments.find(lnoder);
                if(rslt != end_segments.end()) {
                    auto lst = rslt->second;
                    if(lst.size() > 1 || end_kmers_count[lnoder] > 1) { // multiple connection - clip kmer from contig and keep segment
                        seq.erase(seq.begin(), seq.begin()+kmer_len);
                    } else {                                            // one contig, one segment - include segment in contig and mark segment for deletion
                        SpiderSegmentPDir& segmp = lst.front();
                        auto& seg_seq = segmp.m_segmi->m_seq;
                        int seg_len = seg_seq.size();
                        if(seg_len <= kmer_len) {
                            seq.erase(seq.begin(), seq.begin()+kmer_len-seg_len);
                        } else {                            
                            string s;
                            if(segmp.m_right_end) {                        // left contig to right segment - same strands
                                for(int i = 0; i < seg_len-kmer_len; ++i)
                                    s.push_back(seg_seq[i].m_nt);
                            } else {                                       // left contig to left segment - different strands                 
                                for(int i = kmer_len; i < seg_len; ++i)
                                    s.push_back(seg_seq[i].m_nt);
                                ReverseComplementSeq(s.begin(), s.end());
                            }
                            seq = s+seq;
                        }
                        erased_segments.insert(segmp);
                        seg_seq.clear();
                    }
                    for(SpiderSegmentPDir& segmp : lst)
                        segment_links[segmp].emplace_back(acc, false);
                    auto& lcl = left_contig_links[acc];
                    lcl.splice(lcl.end(), lst);
                }
            }
        }
                      
        ostream& out = gfa_out.is_open() ? gfa_out : cout;
        for(auto& contig : contigs) {
            auto& acc = contig.first;
            auto& seq = contig.second;
            CReadHolder rh(false);
            rh.PushBack(seq);
            size_t count = 0;
            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len); ik != rh.kend(); ++ik)
                count += graphp->Abundance(graphp->GetNode(*ik));
            out << "S\t" << acc << "\t" << seq <<  "\tKC:i:" << count << "\n";

            if(right_contig_links.count(acc)) {
                for(SpiderSegmentPDir& segmp : right_contig_links[acc]) {
                    list<GFAIterator> linked_segments;
                    list<GFAIterator> dead_segments;
                    if(erased_segments.count(segmp))  { 
                        dead_segments.push_back(segmp.m_segmi);
                        auto& connections = segmp.m_right_end ? segmp.m_segmi->m_left_connections : segmp.m_segmi->m_right_connections;
                        for(auto cnt : connections) {
                            if(erased_segments.count(SpiderSegmentP(segmp.m_spideri, cnt))) 
                                dead_segments.push_back(cnt);
                            else
                                linked_segments.push_back(cnt);
                        }
                    } else {
                        linked_segments.push_back(segmp.m_segmi);
                    }
                    list<tuple<string,bool>> linked_contigs;
                    for(auto ds : dead_segments) {
                        auto rslt = segment_links.find(SpiderSegmentPDir(segmp.m_spideri, ds, !segmp.m_right_end));
                        if(rslt != segment_links.end())
                            linked_contigs.insert(linked_contigs.end(), rslt->second.begin(), rslt->second.end());
                    }

                    for(GFAIterator i : linked_segments) { //connections to spiders
                        char dir = segmp.m_right_end ? '-' : '+';
                        out << "L\t" << acc << "\t+\t" << segmp.m_spideri->SegId(*i) << "\t" << dir << "\t0M\n";
                    }
                    for(auto& lc : linked_contigs) {       //direct connections to other contigs
                        char dir = get<1>(lc) ? '-' : '+';
                        out << "L\t" << acc << "\t+\t" << get<0>(lc) << "\t" << dir << "\t0M\n";
                    }
                }
            }

            if(left_contig_links.count(acc)) {
                for(SpiderSegmentPDir& segmp : left_contig_links[acc]) {
                    list<GFAIterator> linked_segments;
                    list<GFAIterator> dead_segments;
                    if(erased_segments.count(segmp))  { 
                        dead_segments.push_back(segmp.m_segmi);
                        auto& connections = segmp.m_right_end ? segmp.m_segmi->m_left_connections : segmp.m_segmi->m_right_connections;
                        for(auto cnt : connections) {
                            if(erased_segments.count(SpiderSegmentP(segmp.m_spideri, cnt))) 
                                dead_segments.push_back(cnt);
                            else
                                linked_segments.push_back(cnt);
                        }
                    } else {
                        linked_segments.push_back(segmp.m_segmi);
                    }
                    list<tuple<string,bool>> linked_contigs;
                    for(auto ds : dead_segments) {
                        auto rslt = segment_links.find(SpiderSegmentPDir(segmp.m_spideri, ds, !segmp.m_right_end));
                        if(rslt != segment_links.end())
                            linked_contigs.insert(linked_contigs.end(), rslt->second.begin(), rslt->second.end());
                    }                    

                    for(GFAIterator i : linked_segments) { //connections to spiders
                        char dir = segmp.m_right_end ? '-' : '+';
                        out << "L\t" << acc << "\t-\t" << segmp.m_spideri->SegId(*i) << "\t" << dir << "\t0M\n";
                    }
                    for(auto& lc : linked_contigs) {       //direct connections to other contigs
                        char dir = get<1>(lc) ? '-' : '+';
                        out << "L\t" << acc << "\t-\t" << get<0>(lc) << "\t" << dir << "\t0M\n";
                    }
                }
            }
        }        

        for(auto& segmp : erased_segments) {
            segmp.m_spideri->RemoveSegment(segmp.m_segmi);
            if(segmp.m_spideri->empty()) {
                cerr << "ErasedC " << segmp.m_spideri->Target() << endl;
                spiders.erase(segmp.m_spideri);
            }
        }
        for(auto& spider : spiders)
            spider.PrintGFA(out);

        if(gfa_out.is_open()) {
            gfa_out.close();
            if(!gfa_out) {
                cerr << "Can't write to file " << argm["gfa_out"].as<string>() << endl;
                exit(1);
            }
        } else {
            cout.flush();
            if(!cout) {
                cerr << "Write failed " << endl;
                exit(1);
            }
        }

        if(csv_out.is_open()) {
            csv_out << "Name,Color\n";
            for(auto& contig : contigs)
                csv_out << contig.first << "," <<  argm["contig_color"].as<string>() << "\n";
            for(auto& spider : spiders) {
                for(auto& segm : spider)
                    csv_out << spider.SegId(segm)  << "," <<  argm["connector_color"].as<string>() << "\n";
            }
            csv_out.close();
            if(!csv_out) {
                cerr << "Can't write to file " << argm["csv_out"].as<string>() << endl;
                exit(1);
            }
        }

        cerr << "Combined GFA in " << timer.Elapsed();
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        cerr << all << "\n";
        return 1;
    }
    return 0;
}

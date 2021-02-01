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

#include "readsgetter.hpp"
#include "assembler.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

template <class DBGraph>
//void PrintRslt(CDBGAssembler<DBGraph>& assembler, ofstream& contigs_out, ofstream& all_out, ofstream& hist_out, ofstream& connected_reads_out, ofstream& dbg_out, int mincontig) {
void PrintRslt(CDBGAssembler<DBGraph>& assembler, variables_map& argm) {

    ofstream contigs_out;
    ofstream all_out;
    ofstream hist_out;
    ofstream connected_reads_out;
    ofstream dbg_out;

    if(argm.count("contigs_out")) {
        contigs_out.open(argm["contigs_out"].as<string>());
        if(!contigs_out.is_open()) {
            cerr << "Can't open file " << argm["contigs_out"].as<string>() << endl;
            exit(1);
        }
    }
    
    if(argm.count("all")) {
        all_out.open(argm["all"].as<string>());
        if(!all_out.is_open()) {
            cerr << "Can't open file " << argm["all"].as<string>() << endl;
            exit(1);
        }
    }
    
    if(argm.count("hist")) {
        hist_out.open(argm["hist"].as<string>());
        if(!hist_out.is_open()) {
            cerr << "Can't open file " << argm["hist"].as<string>() << endl;
            exit(1);
        }
    }

    if(argm.count("connected_reads")) {
        connected_reads_out.open(argm["connected_reads"].as<string>());
        if(!connected_reads_out.is_open()) {
            cerr << "Can't open file " << argm["connected_reads"].as<string>() << endl;
            exit(1);
        }
    }

    if(argm.count("dbg_out")) {
        dbg_out.open(argm["dbg_out"].as<string>(), ios::binary | ios::out);
        if(!dbg_out.is_open()) {
            cerr << "Can't open file " << argm["dbg_out"].as<string>() << endl;
            exit(1);
        }
    }

    int mincontig = argm["min_contig"].as<int>();
    if(mincontig <= 0) {
        cerr << "Value of --min_contig must be > 0" << endl;
        exit(1);
    }
    
    DBGraph& first_graph = *assembler.Graphs().begin()->second;
    int first_kmer_len = first_graph.KmerLen();
    int num = 0; 
    ostream& out = contigs_out.is_open() ? contigs_out : cout;
    auto contigs = assembler.Contigs();

    contigs.sort();
    for(auto& contig : contigs) {
        if((int)contig.LenMin() >= mincontig) {
            deque<list<pair<double, string>>> scored_contig;
            for(unsigned chunk = 0; chunk < contig.size(); ++chunk) {
                scored_contig.emplace_back();
                if(contig.VariableChunk(chunk)) {
                    double total_abundance = 0.;
                    for(auto& variant : contig[chunk]) {
                        TVariation seq = variant;
                        if(chunk < contig.size()-1) {
                            auto a = contig[chunk+1].front().begin();
                            auto b = contig[chunk+1].front().end();
                            if((int)contig.ChunkLenMax(chunk+1) > first_kmer_len-1)
                                b = a+first_kmer_len-1;
                            seq.insert(seq.end(), a, b);
                        }
                        if(chunk > 0) {
                            auto b = contig[chunk-1].front().end();
                            auto a = contig[chunk-1].front().begin();
                            if((int)contig.ChunkLenMax(chunk-1) > first_kmer_len-1)
                                a = b-first_kmer_len+1;
                            seq.insert(seq.begin(), a, b);
                        }
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        double abundance = 0;
                        for(CReadHolder::kmer_iterator itk = rh.kbegin(first_graph.KmerLen()); itk != rh.kend(); ++itk) {
                            typename DBGraph::Node node = first_graph.GetNode(*itk);
                            abundance += first_graph.Abundance(node);
                        }
                        total_abundance += abundance;
                        double score = abundance;
                        string var_seq(variant.begin(), variant.end());
                        scored_contig.back().emplace_back(score, var_seq);
                    }

                    for(auto& score_seq : scored_contig.back())
                        score_seq.first /= total_abundance;
                    scored_contig.back().sort();
                    scored_contig.back().reverse();
                } else {
                    double score = 1.;
                    string var_seq(contig[chunk].front().begin(), contig[chunk].front().end());
                    scored_contig.back().emplace_back(score, var_seq);
                }
            }

            string first_variant;
            for(auto& lst : scored_contig)
                first_variant += lst.front().second;

            CReadHolder rh(false);
            if(contig.m_circular)
                first_variant += first_variant.substr(0, first_graph.KmerLen()-1);
            rh.PushBack(first_variant);
            double abundance = 0; // average count of kmers in contig   
            for(CReadHolder::kmer_iterator itk = rh.kbegin(first_graph.KmerLen()); itk != rh.kend(); ++itk) {
                typename DBGraph::Node node = first_graph.GetNode(*itk);
                abundance += first_graph.Abundance(node);
            }
            abundance /= first_variant.size()-first_graph.KmerLen()+1;
            out << ">Contig_" << ++num << "_" << abundance;
            if(contig.m_circular) {
                out << "_Circ [topology=circular]";
                first_variant.erase(first_variant.size()-first_graph.KmerLen()+1, first_graph.KmerLen()-1); 
            }
            out << "\n" << first_variant << "\n";
            int pos = 0;
            for(unsigned chunk = 0; chunk < scored_contig.size(); ++chunk) { //output variants  
                int chunk_len = scored_contig[chunk].front().second.size();
                if(contig.VariableChunk(chunk)) {
                    int left = 0;
                    if(chunk > 0)
                        left = min(100,(int)scored_contig[chunk-1].front().second.size());
                    int right = 0;
                    if(chunk < scored_contig.size()-1)
                        right = min(100,(int)scored_contig[chunk+1].front().second.size());
                    int var = 0;
                    auto it = scored_contig[chunk].begin();
                    for(++it; it != scored_contig[chunk].end(); ++it) {
                        double score = it->first;
                        string& variant = it->second;
                        out << ">Variant_" << ++var << "_for_Contig_" << num << ":" << pos-left+1 << "_" << pos+chunk_len+right << ":" << score << "\n";
                        if(chunk > 0) {                                
                            for(int l = left ; l > 0; --l)
                                out << *(scored_contig[chunk-1].front().second.end()-l);
                        }
                        out << variant;
                        if(chunk < scored_contig.size()-1) {
                            for(int r = 0; r < right; ++r)
                                out << scored_contig[chunk+1].front().second[r];
                        }
                        out << "\n";
                    }
                }
                pos += chunk_len;
            }            
        } 
    } 

    if(contigs_out.is_open()) {
        contigs_out.close();
        if(!contigs_out) {
            cerr << "Can't write to file " << argm["contigs_out"].as<string>() << endl;
            exit(1);
        }
    } else {
        cout.flush();
        if(!cout) {
            cerr << "Write failed " << endl;
            exit(1);
        }
    }

    if(all_out.is_open()) {
        auto graphp = assembler.Graphs().begin();
        auto it = assembler.AllIterations().begin();
        if(argm.count("seeds")) {
            auto& contigs = *it;
            int nn = 0;
            for(auto& contig : contigs) {
                string first_variant;
                for(auto& lst : contig)
                    first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                all_out << ">Seed_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
            }
            ++it;
        }
        for( ; graphp != assembler.Graphs().end(); ++it, ++graphp) {
            auto& contigs = *it;
            int nn = 0;
            for(auto& contig : contigs) {
                string first_variant;
                for(auto& lst : contig)
                    first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                all_out << ">kmer" << graphp->first << "_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
            }
        }
        if(argm.count("allow_snps")) {
            auto graphpr = assembler.Graphs().rbegin();
            for( ; graphpr != assembler.Graphs().rend(); ++it, ++graphpr) {
                auto& contigs = *it;
                int nn = 0;
                for(auto& contig : contigs) {
                    string first_variant;
                    for(auto& lst : contig)
                        first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                    all_out << ">SNP_recovery_kmer" << graphpr->first << "_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
                }
            } 
        } 

        all_out.close();
        if(!all_out) {
            cerr << "Can't write to file " << argm["all"].as<string>() << endl;
            exit(1);
        }
    }

    if(hist_out.is_open()) {
        for(auto& gr : assembler.Graphs()) {
            const TBins& bins = gr.second->GetBins();
            for(auto& bin : bins)
                hist_out << gr.first << '\t' << bin.first << '\t' << bin.second << "\n";
        }
        hist_out.close();
        if(!hist_out) {
            cerr << "Can't write to file " << argm["hist"].as<string>() << endl;
            exit(1);
        }
    }

    if(connected_reads_out.is_open()) {
        CReadHolder connected_reads = assembler.ConnectedReads();
        int num = 0;
        for(CReadHolder::string_iterator is = connected_reads.sbegin(); is != connected_reads.send(); ++is) {
            string s = *is;
            connected_reads_out << ">ConnectedRead_" << ++num << "\n" << s << "\n";
        }
        connected_reads_out.close();
        if(!connected_reads_out) {
            cerr << "Can't write to file " << argm["connected_reads"].as<string>() << endl;
            exit(1);
        }
    }
    
    if(dbg_out.is_open()) {
        for(auto& gr : assembler.Graphs())
            gr.second->Save(dbg_out);
        dbg_out.close();
        if(!dbg_out) {
            cerr << "Can't write to file " << argm["dbg_out"].as<string>() << endl;
            exit(1);
        }
    }    
}


int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    int ncores;
    int steps;
    double fraction;
    double vector_percent;
    int jump;
    int min_count;
    int min_kmer;
    int max_kmer = 0;
    bool usepairedends;
    bool forcesinglereads;
    int maxkmercount;
    int max_kmer_paired = 0;
    vector<string> sra_list;
    vector<string> file_list;
    bool allow_snps;
    bool estimate_min_count = true;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]")
        ("memory", value<int>()->default_value(32), "Memory available (GB, only for sorted counter) [integer]")
        ("hash_count", "Use hash counter [flag]")
        ("estimated_kmers", value<int>()->default_value(100), "Estimated number of distinct kmers for bloom filter (millions, only for hash counter) [integer]")
        ("skip_bloom_filter", "Don't do bloom filter; use --estimated_kmers as the hash table size (only for hash counter) [flag]");

    options_description input("Input/output options : at least one input providing reads for assembly must be specified");
    input.add_options()
        ("reads", value<vector<string>>(), "Input fasta/fastq file(s) for reads (could be used multiple times for different runs, could be gzipped) [string]")
        ("use_paired_ends", "Indicates that single (not comma separated) fasta/fastq files contain paired reads [flag]")
#ifndef NO_NGS
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
#endif
        ("contigs_out", value<string>(), "Output file for contigs (stdout if not specified) [string]");

    options_description assembly("Assembly options");
    assembly.add_options()
        ("kmer", value<int>()->default_value(21), "Minimal kmer length for assembly [integer]")
        ("min_count", value<int>(), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("max_kmer", value<int>(), "Maximal kmer length for assembly [integer]")
        ("max_kmer_count", value<int>(), "Minimum acceptable average count for estimating the maximal kmer length in reads [integer]")
        ("vector_percent", value<double>()->default_value(0.05, "0.05"), "Percentage of reads containing 19-mer for the 19-mer to be considered a vector (1. disables) [float (0,1]]")
        ("insert_size", value<int>(), "Expected insert size for paired reads (if not provided, it will be estimated) [integer]")
        ("steps", value<int>()->default_value(11), "Number of assembly iterations from minimal to maximal kmer length in reads [integer]")
        ("fraction", value<double>()->default_value(0.1, "0.1"), "Maximum noise to signal ratio acceptable for extension [float [0,1)]")
        ("max_snp_len", value<int>()->default_value(150), "Maximal snp length [integer]")
        ("min_contig", value<int>()->default_value(200), "Minimal contig length reported in output [integer]")
        ("allow_snps", "Allow additional step for snp discovery [flag]");

    options_description debug("Debugging options");
    debug.add_options()
        ("force_single_ends", "Don't use paired-end information [flag]")
        ("seeds", value<string>(), "Input file with seeds [string]")
        ("all", value<string>(), "Output fasta for each iteration [string]")
        ("dbg_out", value<string>(), "Output kmer file [string]")
        ("hist", value<string>(), "File for histogram [string]")
        ("connected_reads", value<string>(), "File for connected paired reads [string]");

    options_description deprecated("");
    deprecated.add_options()
        ("fasta", value<vector<string>>(), "Input fasta file(s) (could be used multiple times for different runs) [string]")
        ("fastq", value<vector<string>>(), "Input fastq file(s) (could be used multiple times for different runs) [string]")
        ("gz", "Input fasta/fastq files are gzipped [flag]");

    options_description all("");
    all.add(general).add(input).add(assembly).add(debug).add(deprecated); 

    options_description visible("");
    visible.add(general).add(input).add(assembly).add(debug);

    /*
    // gcc memory menagement
    if (mallopt(M_MMAP_THRESHOLD, 128*1024) != 1)
        cerr << "M_MMAP_THRESHOLD failed" << endl;
    if (mallopt(M_ARENA_MAX, 100) != 1)
        cerr << "M_ARENA_MAX failed" << endl;
    */

    try {
        variables_map argm;                                // boost arguments
        store(parse_command_line(argc, argv, all), argm);
        notify(argm);    

        if(argm.count("gz"))
            cerr << "WARNING: option --gz is deprecated - gzipped files are now recognized automatically" << endl;       

        if(argc == 1 || argm.count("help")) {
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cout << visible << "\n";
            return 0;
        }

        if(argm.count("version")) {
            cout << "SKESA 2.5.0" << endl;
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl;
#endif
            return 0;
        }

        if(!argm.count("reads") && !argm.count("fasta") && !argm.count("fastq") 
#ifndef NO_NGS
           && !argm.count("sra_run")
#endif
                                                     ) {
            cerr << "Provide some input reads" << endl;
            cerr << visible << "\n";
            return 1;
        }

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
            file_list = argm["reads"].as<vector<string>>(); 
        }       
        if(argm.count("fasta")) {
            auto& lst = argm["fasta"].as<vector<string>>();
            file_list.insert(file_list.end(), lst.begin(), lst.end());
        }
        if(argm.count("fastq")) {
            auto& lst = argm["fastq"].as<vector<string>>();
            file_list.insert(file_list.end(), lst.begin(), lst.end());
        }

        if(!file_list.empty()){
            unsigned num = file_list.size();
            sort(file_list.begin(), file_list.end());
            file_list.erase(unique(file_list.begin(),file_list.end()), file_list.end());
            if(file_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from file list" << endl; 
        }

        allow_snps = argm.count("allow_snps");
   
        ncores = thread::hardware_concurrency();
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

        steps = argm["steps"].as<int>();
        if(steps <= 0) {
            cerr << "Value of --steps must be > 0" << endl;
            exit(1);
        }
        fraction = argm["fraction"].as<double>();
        if(fraction >= 1.) {
            cerr << "Value of --fraction must be < 1 (more than 0.25 is not recommended)" << endl;
            exit(1);
        }
        if(fraction < 0.) {
            cerr << "Value of --fraction must be >= 0" << endl;
            exit(1);
        }
        jump = argm["max_snp_len"].as<int>();
        if(jump < 0) {
            cerr << "Value of --max_snp_len must be >= 0" << endl;
            exit(1);
        }
        if(argm.count("insert_size"))
            max_kmer_paired = argm["insert_size"].as<int>();

        min_count = 2;
        if(argm.count("min_count")) {
            min_count = argm["min_count"].as<int>();
            estimate_min_count = false;
        }
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }
        maxkmercount = 10;
        if(argm.count("max_kmer_count")) {
            maxkmercount = argm["max_kmer_count"].as<int>();
            estimate_min_count = false;
        }
        if(maxkmercount <= 0) {
            cerr << "Value of --max_kmer_count must be > 0" << endl;
            exit(1);
        }

        if(argm.count("max_kmer")) {
            max_kmer = argm["max_kmer"].as<int>();
            if(max_kmer <= 0) {
                cerr << "Value of --max_kmer must be > 0" << endl;
                exit(1);
            }            
        }

        if(max_kmer_paired < 0) {
            cerr << "Value of --insert_size must be >= 0" << endl;
            exit(1);
        }
        min_kmer = argm["kmer"].as<int>();
        if(min_kmer < 21 || min_kmer%2 ==0) {
            cerr << "Kmer must be an odd number >= 21" << endl;
            return 1;
        }
        vector_percent = argm["vector_percent"].as<double>();
        if(vector_percent > 1.) {
            cerr << "Value of --vector_percent  must be <= 1" << endl;
            exit(1);
        }
        if(vector_percent <= 0.) {
            cerr << "Value of --vector_percent  must be > 0" << endl;
            exit(1);
        }

        usepairedends = argm.count("use_paired_ends");
        forcesinglereads = argm.count("force_single_ends");

        TStrList seeds;
        if(argm.count("seeds")) {
            ifstream seeds_in;
            seeds_in.open(argm["seeds"].as<string>());
            if(!seeds_in.is_open()) {
                cerr << "Can't open file " << argm["seeds"].as<string>() << endl;
                exit(1);
            }
            char c;
            if(!(seeds_in >> c)) {
                cerr << "Empty fasta file for seeds" << endl;
            } else if(c != '>') {
                cerr << "Invalid fasta file format in " << argm["seeds"].as<string>() << endl;
                exit(1);
            }

            string record;
            while(getline(seeds_in, record, '>')) {                       
                size_t first_ret = min(record.size(),record.find('\n'));
                if(first_ret == string::npos) {
                    cerr << "Invalid fasta file format in " << argm["seeds"].as<string>() << endl;
                    exit(1);
                }
                string sequence = record.substr(first_ret+1);
                sequence.erase(remove(sequence.begin(),sequence.end(),'\n'), sequence.end()); 
                if(sequence.find_first_not_of("ACGTYRWSKMDVHBN") != string::npos) {
                    cerr << "Invalid fasta file format in " << argm["seeds"].as<string>() << endl;
                    exit(1);
                } 
                seeds.push_back(sequence);
            }
        }

        int low_count = max(min_count, 2); 
        CReadsGetter readsgetter(sra_list, file_list, ncores, usepairedends);

        if(argm.count("hash_count")) {
            int estimated_kmer_num =  argm["estimated_kmers"].as<int>();
            if(estimated_kmer_num <= 0) {
                cerr << "Value of --estimated_kmers must be > 0" << endl;
                exit(1);
            }
            bool skip_bloom_filter = argm.count("skip_bloom_filter");
            if(vector_percent < 1.) {
                readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmer_num, skip_bloom_filter);
                readsgetter.PrintAdapters();
            } else {
                cerr << "Adapters clip is disabled" << endl;
            }
            CDBGAssembler<CDBHashGraph> assembler(fraction, jump, low_count, steps, min_count, min_kmer, max_kmer, forcesinglereads, max_kmer_paired, 
                                                  maxkmercount, ncores, readsgetter.Reads(), seeds, allow_snps, estimate_min_count, 
                                                  estimated_kmer_num, skip_bloom_filter); 
            PrintRslt(assembler, argm);
        } else {
            int memory = argm["memory"].as<int>();
            if(memory <= 0) {
                cerr << "Value of --memory must be > 0" << endl;
                exit(1);
            }
            if(vector_percent < 1.) {
                readsgetter.ClipAdaptersFromReads_SortedCounter(vector_percent, memory);
                readsgetter.PrintAdapters();
            } else {
                cerr << "Adapters clip is disabled" << endl;
            }
            CDBGAssembler<CDBGraph> assembler(fraction, jump, low_count, steps, min_count, min_kmer, max_kmer, forcesinglereads, max_kmer_paired, 
                                              maxkmercount, ncores, readsgetter.Reads(), seeds, allow_snps, estimate_min_count, memory); 
            PrintRslt(assembler, argm);
        }

        cerr << "DONE" << endl;
        exit(0);
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }


    return 0;
}

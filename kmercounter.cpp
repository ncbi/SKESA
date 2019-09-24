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
#include "concurrenthash.hpp"
#include "graphdigger.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    int kmer;
    int ncores;
    double vector_percent;
    int min_count;
    vector<string> sra_list;
    vector<string> file_list;
    size_t estimated_kmer_num;

    options_description all("Options");
    all.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
#ifndef NO_NGS
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
#endif
        ("reads", value<vector<string>>(), "Input fasta/fastq file(s) (could be used multiple times for different runs) [string]")
        ("kmer", value<int>()->default_value(21), "Kmer length [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("vector_percent", value<double>()->default_value(0.05, "0.05"), "Count for  vectors as a fraction of the read number (1. disables) [float (0,1]]")

        ("estimated_kmers", value<int>()->default_value(100), "Estimated number of unique kmers for bloom filter (millions) for hash count [integer]")
        ("skip_bloom_filter", "Don't do bloom filter; use --estimated_kmers as the hash table size for hash count [flag]")

        ("dbg_out", value<string>(), "De Bruijn graph output")
        ("text_out", value<string>(), "Text kmer output")
        ("hist", value<string>(), "File for histogram [string]")

        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]");

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
            cout << "kmercounter v.2.1.0";
#ifdef SVN_REV
            cout << "-SVN_" << SVN_REV;
#endif
            cout << endl;
            return 0;
        }

        if(!argm.count("reads")
#ifndef NO_NGS
           && !argm.count("sra_run")
#endif
                                                     ) {
            cerr << "Provide some input reads" << endl;
            cerr << all << "\n";
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
            unsigned num = file_list.size();
            sort(file_list.begin(), file_list.end());
            file_list.erase(unique(file_list.begin(),file_list.end()), file_list.end());
            if(file_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from file list" << endl; 
        }

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

        min_count = argm["min_count"].as<int>();
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }

        estimated_kmer_num =  argm["estimated_kmers"].as<int>();
        kmer = argm["kmer"].as<int>();

        CReadsGetter readsgetter(sra_list, file_list, vector<string>(), ncores, false);
        vector_percent = argm["vector_percent"].as<double>();
        if(vector_percent > 1.) {
            cerr << "Value of --vector_percent  must be <= 1" << endl;
            exit(1);
        }
        if(vector_percent <= 0.) {
            cerr << "Value of --vector_percent  must be > 0" << endl;
            exit(1);
        }
 
        if(vector_percent < 1.) {
            readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmer_num, argm.count("skip_bloom_filter"));
            readsgetter.PrintAdapters();
        } else {
            cerr << "Adapters clip is disabled" << endl;
        }

        size_t MB = 1000000;
        CKmerHashCounter counter(readsgetter.Reads(), kmer, min_count, estimated_kmer_num*MB, true, ncores, argm.count("skip_bloom_filter"));

        if(argm.count("text_out")) {
            ofstream out(argm["text_out"].as<string>());
            if(!out.is_open()) {
                cerr << "Can't open file " << argm["text_out"].as<string>() << endl;
                return 1;             
            }
            
            CKmerHashCount& hash = counter.Kmers();
            for(auto index = hash.Begin(); index != hash.End(); ++index) {
                auto rslt = index.GetElement();
                out << rslt.first.toString(kmer) << "\t" << rslt.second->Count() << "\t" << (rslt.second->m_data.Load() >> 32) << endl;
            }
            out.close();
            if(!out) {
                cerr << "Can't write to file " << argm["text_out"].as<string>() << endl;
                return 1;             
            }
        }

        if(argm.count("hist")) {
            ofstream out(argm["hist"].as<string>());
            if(!out.is_open()) {
                cerr << "Can't open file " << argm["hist"].as<string>() << endl;
                return 1;             
            }
            TBins bins = counter.Kmers().GetBins();
            for(auto& bin : bins)
                out << bin.first << '\t' << bin.second << endl;
            out.close();
            if(!out) {
                cerr << "Can't write to file " << argm["hist"].as<string>() << endl;
                return 1;             
            }
        }

        if(argm.count("dbg_out")) {
            counter.GetBranches();
            CDBHashGraph graph(move(counter.Kmers()), true);
            ofstream dbg_out(argm["dbg_out"].as<string>(), ios::binary | ios::out);
            if(!dbg_out.is_open()) {
                cerr << "Can't open file " << argm["dbg_out"].as<string>() << endl;
                exit(1);
            }
            graph.Save(dbg_out);
            dbg_out.close();
            if(!dbg_out) {
                cerr << "Can't write to file " << argm["dbg_out"].as<string>() << endl;
                return 1;             
            }
        }                
        
        cerr << "DONE" << endl;
        exit(0);
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }

}


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

#include <unordered_map>
#include <math.h>
#include "readsgetter.hpp"
#include "counter.hpp"
#include "gfa.hpp"
#include "guidedpath_naa.hpp"
#include "guidedgraph.hpp"
#include "guidedassembler.hpp"

using namespace boost::program_options;
using namespace std;
using namespace DeBruijn;

int main(int argc, const char* argv[]) {

    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]")
        ("estimated_kmers", value<int>()->default_value(1000), "Estimated number of distinct kmers for bloom filter (millions) [integer]");

    options_description input("Input/output options : target file, at least one input for reads and output for gfa must be specified");
    input.add_options()
        ("gfa", value<string>(), "Output file for GFA graph [string]")
        ("all_variants", value<string>(), "Output file for sequences with all variant combinations [string]")
        ("max_variants", value<int>()->default_value(1000), "Restricts the number of variants outputted in --all_variants [integer]")
        ("selected_variants", value<string>(), "Output file for selected sequences representing all graph elements [string]")
        ("targets", value<string>(), "Input file with reference sequences [string]")
#ifndef NO_NGS
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
#endif
        ("reads", value<vector<string>>(), "Input fasta/fastq file(s) for reads (could be used multiple times for different runs, could be gzipped) [string]")
        ("use_paired_ends", "Indicates that single (not comma separated) fasta/fastq files contain paired reads [flag]");

    options_description assembly("Assembly options");
    assembly.add_options()
        ("extend_ends", "Unambiguously extend graph ends using de-novo assembly [flag]")
        ("kmer", value<int>(), "Kmer length for assembly (default automatic) [integer]")
        ("secondary_kmer", value<int>(), "Shorter kmer length for low coverage spots (default automatic) [integer]")
        ("secondary_kmer_threshold", value<int>()->default_value(1), "Coverage threshold for using shorter kmer [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained in graph [integer]")
        ("vector_percent", value<double>()->default_value(0.05, "0.05"), "Percentage of reads containing 19-mer for the 19-mer to be considered a vector (1. disables) [float (0,1]]")
        ("fraction", value<double>()->default_value(0.05, "0.05"), "Maximum noise to signal ratio acceptable for extension [float [0,1)]")
        ("word", value<int>()->default_value(8), "Word size for seeds [integer <= 16]")
        ("protect_reference_ends", "Near the reference ends, don't check if reads support some minimal extension of the fork's branches [flag]")        
        ("keep_subgraphs", "Don't remove redundant subgraphs [flag]")        
        ("match", value<int>()->default_value(2), "Bonus for match [integer]")
        ("mismatch", value<int>()->default_value(3), "Penalty for mismatch [integer]")
        ("gap_open", value<int>()->default_value(5), "Penalty for gap opening [integer]")
        ("gap_extend", value<int>()->default_value(3), "Penalty for gap extension [integer]")
        ("drop_off", value<int>()->default_value(300), "Maximal decrease of score before abandoning assembly path [integer]")
        ("kmer_complexity", value<int>()->default_value(2000), "Hard mask reference areas with hight kmer complexity (0 disables masking) [integer]")
        ("max_fork_density", value<double>()->default_value(0.1, "0.1"), "Maximal fork density averaged over --buf_length before abandoning assembling (0 disables) [float]")
        ("buf_length", value<int>()->default_value(200), "Buffer length for fork density [integer]");

    options_description filter("Graph cleaning options");
    filter.add_options()
        ("not_aligned_len", value<int>()->default_value(10), "Not aligned read length for break count [integer]")
        ("not_aligned_count", value<int>()->default_value(3), "Number of not aligned reads to make a break [integer]")
        ("aligned_count", value<int>()->default_value(2), "Number of aligned reads to confirm a connection [integer]")
        ("target_coverage", value<double>()->default_value(0.5, "0.5"), "Result is kept if has alignment to the reference which is at least target_coverage*(reference length) long [float (0,1]]")
        ("min_hit_len", value<int>(), "If a path is shorter than target_coverage*(reference length), use this length threshold for keeping paths (optional) [integer]")
        ("max_path", value<int>()->default_value(1000), "Maximal number of path extensions allowed for a single filtering check [integer]")
        ("no_filter_by_reads", "Don't use full length reads for variants filtering [flag]")
        ("no_filter_by_pairs", "Don't use mate pairs for variants filtering [flag]")
        ("remove_homopolymer_indels", "Remove homopolymer indels [flag]")
        ("homopolymer_len", value<int>()->default_value(4), "Minimal length of homopolymer [integer]")
        ("homopolymer_ratio", value<double>()->default_value(0.33, "0.33"), "Coverage ratio threshold for removing homopolymer indels [float]")
        ("use_ambiguous_na", "Use ambiguous nucleotide codes for SNPs [flag]");

    options_description all("");
    all.add(general).add(input).add(assembly).add(filter); 

    CStopWatch timer;
    timer.Stop();
    try {
        variables_map argm;                                // boost arguments
        store(parse_command_line(argc, argv, all), argm);
        notify(argm);    

        if(argc == 1 || argm.count("help")) {
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cout << all << "\n";
            return 0;
        }

        if(argm.count("version")) {
            cout << "saute 1.3.0" << endl;
#ifdef SVN_REV
            cout << "SVN revision:" << SVN_REV << endl;
#endif
            return 0;
        }

        ofstream gfa_out;
        if(!argm.count("gfa")) {
            cerr << "Provide file for GFA graph" << endl;
            cerr << all << "\n";
            return 1;
        } else {
            gfa_out.open(argm["gfa"].as<string>());
            if(!gfa_out.is_open()) {
                cerr << "Can't open file " << argm["gfa"].as<string>() << endl;
                exit(1);
            }
        }

        ofstream all_variants;
        if(argm.count("all_variants")) {
            all_variants.open(argm["all_variants"].as<string>());
            if(!all_variants.is_open()) {
                cerr << "Can't open file " << argm["all_variants"].as<string>() << endl;
                exit(1);
            }
        }
        int max_variants = argm["max_variants"].as<int>();

        ofstream selected_variants;
        if(argm.count("selected_variants")) {
            selected_variants.open(argm["selected_variants"].as<string>());
            if(!selected_variants.is_open()) {
                cerr << "Can't open file " << argm["selected_variants"].as<string>() << endl;
                exit(1);
            }
        }

        ifstream targets_in;
        if(!argm.count("targets")) {
            cerr << "Provide target sequences" << endl;
            cerr << all << "\n";
            return 1;
        } else {
            targets_in.open(argm["targets"].as<string>());
            if(!targets_in.is_open()) {
                cerr << "Can't open file " << argm["targets"].as<string>() << endl;
                exit(1);
            }
        }

        bool no_reads = argm.count("no_filter_by_reads");
        bool no_pairs = argm.count("no_filter_by_pairs");
        bool use_ambiguous = argm.count("use_ambiguous_na");
        
        if(!argm.count("reads") 
#ifndef NO_NGS
           && !argm.count("sra_run")
#endif
                                                                  ) {
            cerr << "Provide some input reads" << endl;
            cerr << all << "\n";
            return 1;
        }

        vector<string> sra_list;
        vector<string> file_list;

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

        int match = argm["match"].as<int>();
        int mismatch = argm["mismatch"].as<int>();
        int gap_open = argm["gap_open"].as<int>();
        int gap_extend = argm["gap_extend"].as<int>();
        int drop_off = argm["drop_off"].as<int>();
        int word_size = min(16, argm["word"].as<int>());

        int not_aligned_len = argm["not_aligned_len"].as<int>(); 
        int not_aligned_count = argm["not_aligned_count"].as<int>(); 
        int aligned_count = argm["aligned_count"].as<int>(); 
        int maxp = argm["max_path"].as<int>();
        double target_coverage = argm["target_coverage"].as<double>();
        if(target_coverage > 1.) {
            cerr << "Value of --target_coverage  must be <= 1" << endl;
            exit(1);
        }
        if(target_coverage <= 0.) {
            cerr << "Value of --target_coverage  must be > 0" << endl;
            exit(1);
        }
        int min_hit_len = numeric_limits<int>::max();
        if(argm.count("min_hit_len"))
            min_hit_len = argm["min_hit_len"].as<int>();
        if(min_hit_len <= 0) {
            cerr << "Value of --min_hit_len must be > 0" << endl;
            exit(1);
        } 

        int estimated_kmer_num =  argm["estimated_kmers"].as<int>();
        if(estimated_kmer_num <= 0) {
            cerr << "Value of --estimated_kmers must be > 0" << endl;
            exit(1);
        }

        if(argm.count("kmer") != argm.count("secondary_kmer")) {
            cerr << "Either specify both --kmer and --secondary_kmer or use automatic for both" << endl;
            return 1;
        }
        int kmer = 0;
        int secondary_kmer = 0;
        if(argm.count("kmer")) {
            kmer = argm["kmer"].as<int>();
            if(kmer%2 ==0) {
                cerr << "Kmer must be an odd number" << endl;
                return 1;
            }
            secondary_kmer = argm["secondary_kmer"].as<int>();
            if(secondary_kmer > kmer || secondary_kmer%2 == 0) {
                cerr << "Secondary kmer must be an odd number smaller or equal kmer" << endl;
                return 1;
            }
        }
        int seed_prec = 4;

        list<array<CReadHolder,2>> reads;
        CReadsGetter readsgetter(sra_list, file_list, ncores, usepairedends);
        {
            double length = 0;
            size_t reads_num = 0;
            for(auto& r : readsgetter.Reads()) {
                length += r[0].TotalSeq()+r[1].TotalSeq();
                reads_num += r[0].ReadNum()+r[1].ReadNum();
            }
            int read_len = length/reads_num+0.5;

            seed_prec = read_len/10;
            if(kmer == 0) {
                kmer = read_len/2;
                if(kmer%2 == 0)
                    --kmer;
                secondary_kmer = read_len/5;
                if(secondary_kmer%2 == 0)
                    --secondary_kmer;
                secondary_kmer = max(21, secondary_kmer);
                if(secondary_kmer > kmer) {
                    cerr << "Automatic kmer selection failed. Use --kmer --secondary_kmer" << endl;
                    return 1;
                }
            }
            cerr << "Read length: " << read_len << " Kmer: " << kmer << " Secondary kmer: " << secondary_kmer << " Seed precision: " << seed_prec << endl;
        }
        bool skip_bloom_filter = false;
        if(vector_percent < 1.) {
            readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmer_num, skip_bloom_filter);
            readsgetter.PrintAdapters();
        } else {
            cerr << "Adapters clip is disabled" << endl;
        }
        reads.splice(reads.end(), readsgetter.Reads());        
        
        bool extend_ends = argm.count("extend_ends");
        bool protect_reference_ends = argm.count("protect_reference_ends");
        bool keep_subgraphs = argm.count("keep_subgraphs");
        int kmer_complexity = argm["kmer_complexity"].as<int>();
        double max_fork_density = argm["max_fork_density"].as<double>();
        int buf_length = argm["buf_length"].as<int>();
        int secondary_kmer_threshold = argm["secondary_kmer_threshold"].as<int>();
        bool remove_homopolymer_indels = argm.count("remove_homopolymer_indels");
        int homopolymer_len = argm["homopolymer_len"].as<int>();
        double homopolymer_ratio = argm["homopolymer_ratio"].as<double>();
        CGuidedAssemblerNA gassembler(kmer, secondary_kmer, extend_ends, protect_reference_ends, keep_subgraphs, min_count, fraction, seed_prec, word_size, match, mismatch, gap_open, gap_extend, drop_off, 
                                      kmer_complexity, max_fork_density, buf_length, ncores, reads, targets_in, estimated_kmer_num, skip_bloom_filter, not_aligned_len, 
                                      not_aligned_count, aligned_count, maxp, target_coverage, min_hit_len, no_reads, no_pairs, secondary_kmer_threshold, remove_homopolymer_indels, homopolymer_len, homopolymer_ratio);

        gassembler.PrintRslt(gfa_out, all_variants, max_variants, selected_variants, use_ambiguous);
  
        timer.Resume();
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }

    cerr << "Cleanup and exit in " << timer.Elapsed();
    cerr << "DONE" << endl;

    return 0;
}

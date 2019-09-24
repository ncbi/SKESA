
For questions, please contact
    Alexandre Souvorov (souvorov@ncbi.nlm.nih.gov)
    Richa Agarwala (agarwala@ncbi.nlm.nih.gov)

## Citation

Alexandre Souvorov, Richa Agarwala and David J. Lipman. 
**SKESA: strategic k-mer extension for scrupulous assemblies.**
*Genome Biology* 2018 **19**:153.
[doi.org/10.1186/s13059-018-1540-z](https://doi.org/10.1186/s13059-018-1540-z)

## Compilation

Download current source code

       $ git clone https://github.com/ncbi/SKESA

Alternatively, download last stable release from https://github.com/ncbi/SKESA/releases Releases also include test data and precompiled binary. Test data is available in example subdirectory that has the command in file run.test for generating the SKESA assembly using the test data.

Do following:

       $ cd SKESA


If you would like to build NGS library for accessing reads from SRA, then do

       $ make

Otherwise, if reading inputs only from files, do

       $ make -f Makefile.nongs

BOOST install is expected by makefiles in the SKESA release. If you do not have BOOST on the system path, please specify BOOST_PATH using a command like

       setenv BOOST_PATH /netopt/ncbi_tools64/boost-1.62.0-ncbi1

before running make.

These make files have been tested with BOOST v 1.62.0 and gcc v 4.9.

# SKESA - Strategic K-mer Extension for Scrupulous Assemblies
Version 2.4.0

## Synopsis
    Running
       skesa 
    or 
       skesa -h
    or 
       skesa --help
    gives information about options and produces the following:

    General options:
      -h [ --help ]                 Produce help message
      -v [ --version ]              Print version
      --cores arg (=0)              Number of cores to use (default all) [integer]
      --memory arg (=32)            Memory available (GB, only for sorted counter) 
                                    [integer]
      --hash_count                  Use hash counter [flag]
      --estimated_kmers arg (=100)  Estimated number of unique kmers for bloom 
                                    filter (M, only for hash counter) [integer]
      --skip_bloom_filter           Don't do bloom filter; use --estimated_kmers as
                                    the hash table size (only for hash counter) 
                                    [flag]
    
    Input/output options: at least one input providing reads for assembly must be specified:
      --reads arg                   Input fasta/fastq file(s) for reads (could be 
                                    used multiple times for different runs, could 
                                    be gzipped) [string]
      --use_paired_ends             Indicates that a single (not comma separated) 
                                    fasta/fastq file contains paired reads [flag]
      --sra_run arg                 Input sra run accession (could be used multiple
                                    times for different runs) [string]
      --contigs_out arg             Output file for contigs (stdout if not 
                                    specified) [string]
    
    Assembly options:
      --kmer arg (=21)              Minimal kmer length for assembly [integer]
      --min_count arg               Minimal count for kmers retained for comparing 
                                    alternate choices [integer]
      --max_kmer_count arg          Minimum acceptable average count for estimating
                                    the maximal kmer length in reads [integer]
      --vector_percent arg (=0.05)  Count for  vectors as a fraction of the read 
                                    number (1. disables) [float (0,1]]
      --insert_size arg             Expected insert size for paired reads (if not 
                                    provided, it will be estimated) [integer]
      --steps arg (=11)             Number of assembly iterations from minimal to 
                                    maximal kmer length in reads [integer]
      --fraction arg (=0.1)         Maximum noise to signal ratio acceptable for 
                                    extension [float [0,1)]
      --max_snp_len arg (=150)      Maximal snp length [integer]
      --min_contig arg (=200)       Minimal contig length reported in output 
                                    [integer]
      --allow_snps                  Allow additional step for snp discovery [flag]
    
    Debugging options:
      --force_single_ends           Don't use paired-end information [flag]
      --seeds arg                   Input file with seeds [string]
      --all arg                     Output fasta for each iteration [string]
      --dbg_out arg                 Output kmer file [string]
      --hist arg                    File for histogram [string]
      --connected_reads arg         File for connected paired reads [string]

    Note that --sra_run option is not available if SKESA is built using Makefile.nongs
    

## Short description
SKESA is a de-novo sequence read assembler for microbial genomes based on DeBruijn graphs. It uses conservative heuristics and is designed to create breaks at repeat regions in the genome. This leads to excellent sequence quality. Using k-mers longer than mate length and up to insert size also allows SKESA to attain good contiguity as determined by the N50 statistic. It is a multi-threaded application that scales well with the number of processors. For different runs with the same inputs, including the order of reads, the order and orientation of contigs in the output is deterministic. 

SKESA can process read information by accessing reads from SRA (option `--sra_run`) or from files (option `--reads`) format. Any combination of input streams is allowed. Files could be if fasta or fastq format and gzipped, which is recognized automatically. 
    
When accessing reads from SRA SKESA automatically determines if the read set consists of paired-end or single-end reads. For fasta/fastq input of paired reads with separate files for each mate, filenames separated by a comma for first mate followed by the second mate are listed and in this case, the order of reads is expected to be same in files for both mates. Alternatively, a single file with both mates could be specified. In this case the reads are expected to be interleaved with first mate followed by the second, and the option `--use_paired_ends` must be used.
    
A limitation of the current release is that in case multiple streams of paired reads are provided, it is assumed that all streams have the same insert size. User can explicitly specify expected insert size for the reads (option --insert_size). Otherwise, a sample of input reads is used to estimate the expected insert size. This sampling may lead to very small differences in assembly of the same read set if the order of reads is different and selected sample gives a difference in expected insert size.
    
Two additional options users may wish to specify depending on the resources available to them are as follows:
1. the number of cores (option `--cores`) and
2. total amount of memory in Gb (option `--memory`)

Remaining options are for debugging or modifying algorithm parameters.

Output of assembly is contigs in fasta format. The definition line for contig has format Contig_<N>_<cnt> where <N> is consecutive integers starting from one for numbering the contigs and <cnt> is the average count of k-mers in the contig using minimal k-mer length used in the assembly. Contigs are ordered lexicographically.

Limitations:
1. SKESA is designed for haploid genomes. If it is used for diploid genomes or RNAseq reads, it should create breaks at all heterozygous sites in the genome and sites for alternative splicing, respectively. The allow_snps option can be used to make some joins at well separated heterozygous sites.
2. SKESA is designed for ILLUMINA reads that do not have systematic homopolymer errors. The assembly for reads that do not have properties similar to ILLUMINA is likely to be quite fragmented.
3. Forward-reverse orientation for paired reads is assumed at this time. If this is not true, steps using paired reads are unlikely to change/improve the assembly.
4. Requesting expected insert size to be estimated using a sample is guaranteed to give the same result, including the order of contigs, for the same order of reads but may give very small differences if read order is changed and insert size estimate is different.


## Usage examples

     In all the examples below, we are providing 4 cores and have 48 Gb of memory.

     Example of an assembly that directly accesses SRA for an unpaired read set SRR867211 is:

       $ skesa --sra_run SRR867211 --cores 4 --memory 48 > SRR867211.skesa.fa

     Example of an assembly that directly accesses SRA for a paired read set SRR1960353 is:

       $ skesa --sra_run SRR1960353 --cores 4 --memory 48 > SRR1960353.skesa.fa

     Example of an assembly that uses separate fastq files for each mate of SRR1703350 is:

       $ skesa --reads SRR1703350_1.fq,SRR1703350_2.fq --cores 4 --memory 48 > SRR1703350.skesa.fa

     Example of an assembly that uses interleaved mates for SRR1703350 as fastq input is:

       $ skesa --reads SRR1703350.fq --use_paired_ends --cores 4 --memory 48 > SRR1703350.skesa.fa

     Example of an assembly that uses reads from SRA for SRR1695624 and gzipped fasta for SRR1745628 is:

       $ skesa --sra_run SRR1695624 --reads SRR1745628.fa.gz --use_paired_ends --cores 4 --memory 48 > SAMN03218571.skesa.fa

     Example of the same assembly as above done with both runs accessed from SRA is:

       $ skesa --sra_run SRR1695624 --sra_run SRR1745628 --cores 4 --memory 48 > SAMN03218571.skesa.fa

# GFA_connector - create GFA graph strating from a set of contigs
Version 1.1.0
## Synopsis
    Running
       gfa_connector 
    or 
       gfa_connector -h
    or 
       gfa_connector --help
    gives information about options and produces the following:

    General options:
      -h [ --help ]                        Produce help message
      -v [ --version ]                     Print version
      --cores arg (=0)                     Number of cores to use (default all) 
                                           [integer]
      --estimated_kmers arg (=100)         Estimated number of unique kmers for 
                                           bloom filter (millions) [integer]
      --skip_bloom_filter                  Don't do bloom filter; use 
                                           --estimated_kmers as the hash table size
                                           [flag]

    Input/output options:
      --sra_run arg                        Input sra run accession (could be used 
                                           multiple times for different runs) 
                                           [string]
      --reads arg                          Input fasta/fastq file(s) (could be used
                                           multiple times for different runs, could
                                           be gzipped) [string]
      --use_paired_ends                    Indicates that a single (not comma 
                                           separated) fasta/fastq file contains 
                                           paired reads [flag]
      --contigs arg                        Input file with contigs [string]
      --gfa arg                            GFA graph output (stdout if not 
                                           specified) [string]
      --dbg arg                            Input de Bruijn graph (optional) 
                                           [string]
      --csv arg                            CSV file ouput (optional) [string]
      --contig_color arg (=Purple)         Color for contigs [string]
      --connector_color arg (=Lime Green)  Color for connectors [string]

    Assembly options:
      --kmer arg (=41)                     Kmer length for assembly [integer]
      --min_count arg (=2)                 Minimal count for kmers retained for 
                                           comparing alternate choices [integer]
      --vector_percent arg (=0.05)         Count for vectors as a fraction of the 
                                           read number (1. disables) [float (0,1]]
      --fraction arg (=0.1)                Threshold for extension
      --entropy arg (=0.51)                Minimal entropy for a seed kmer [float]
      --ext_len arg (=2000)                Maximal length for extension

    Graph cleaning options:
      --not_aligned_len arg (=10)          Not aligned read length for break count
      --not_aligned_count arg (=3)         Number of not aligned reads to make a 
                                           break
      --aligned_count arg (=2)             Number of aligned reads to confirm a 
                                           connection
      --max_path arg (=1000)               Maximal number of paths allowed in 1 
                                           step of filtering
      --no_filter_by_reads                 Don't use full length reads for variants
                                           filtering [flag]
      --no_filter_by_pairs                 Don't use mate pairs for variants 
                                           filtering [flag]

    Note that --sra_run option is not available if gfa_connector is built using Makefile.nongs

## Short description
GFA_connector is a SKESA companion which can input the SKESA contigs (or any other contigs for that matter), find supported by reads connections between them and produce a [GFA graph]( http://gfa-spec.github.io/GFA-spec/) which can be viewed and analyzed in [Bandage](https://rrwick.github.io/Bandage/).

GFA_connector constructs a de Bruijn graph for a single k-mer (option `--kmer`) and starts assembling from each end of the input contigs. It continues each possible path up to the maximal length for extension (option `--ext_len`) until it is connected to another contig or starts following another path. In the last case two paths are connected via a branch. At the end of assembling all loose paths not connected to contigs are discarded. This process creates a number of graphs each connecting two or more contigs. 

Before the input contigs and connectors are combined into a single GFA graph for output, the connector graphs are analyzed using the reads. First the input reads are trimmed to the longest stretches compatible with the de Bruijn graph, and then they are aligned to the connectors. Chimeric connections which do not agree with read sequences or paired ends are removed from the graphs. This step could be partially or completely disabled using options `--no_filter_by_reads` and `--no_filter_by_pairs`.

To run GFA_connector one has to provide reads using options `--sra_run` or `--reads` (see SKESA description) and contigs using the option `--contigs`.

## Usage example
The following commands will produce a GFA graph representing a Klebsiella pneumoniae KPR0928 genome (NZ_CP008831.1), two plasmids (NZ_CP008832.1 and NZ_CP008833.1) and a few contigs which the program was not able to connect
```
skesa --sra_run SRR1510963 --hash_count --contigs_out SRR1510963_skesa.fa >& SRR1510963_skesa.log
gfa_connector --sra_run SRR1510963 --contigs SRR1510963_skesa.fa --gfa SRR1510963.gfa --csv SRR1510963.csv --kmer 75 >& SRR1510963_gfa.log
```
# Kmercounter - effectively count k-mers in a read set
Version 2.1.0
## Synopsis
    Running
       kmercounter 
    or 
       kmercounter -h
    or 
       kmercounter --help
    gives information about options and produces the following:
```
Options:
  -h [ --help ]                Produce help message
  -v [ --version ]             Print version
  --sra_run arg                Input sra run accession (could be used multiple 
                               times for different runs) [string]
  --reads arg                  Input fasta/fastq file(s) (could be used 
                               multiple times for different runs) [string]
  --kmer arg (=21)             Kmer length [integer]
  --min_count arg (=2)         Minimal count for kmers retained for comparing 
                               alternate choices [integer]
  --vector_percent arg (=0.05) Count for vectors as a fraction of the read 
                               number (1. disables) [float (0,1]]
  --estimated_kmers arg (=100) Estimated number of unique kmers for bloom 
                               filter (millions) for hash count [integer]
  --skip_bloom_filter          Don't do bloom filter; use --estimated_kmers as 
                               the hash table size for hash count [flag]
  --dbg_out arg                De Bruijn graph output
  --text_out arg               Text kmer output
  --hist arg                   File for histogram [string]
  --cores arg (=0)             Number of cores to use (default all) [integer]
  
     Note that --sra_run option is not available if kmercounter is built using Makefile.nongs
```
## Short description
Kmercounter creates and outputs a de Bruijn graph (option `--dbg_out`) for a specified kmer size (option `--kmer`). This file could be used as an optional input for GFA_connector in which case the kmer counting step will be skipped. 

All accepted kmers and their counts could be output in a text form (option `--text_out`). The tab delimited output file contains:
- Kmer sequence
- Number of reads in which this kmer was seen
- Number of reads in which this kmer was seen in positive direction

The output order is random and could be different for different runs.

With option `--hist` one can get a histogram for the k-mer counts.

## Usage example
```
kmercounter --sra_run SRR1510963 --kmer 75 --dbg_out SRR1510963_k75.hdbg --text_out SRR1510963_k75.kmers --hist SRR1510963_k75.hist >& SRR1510963_counter.log
```

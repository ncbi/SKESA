
For questions, please contact
    Alexandre Souvorov (souvorov@ncbi.nlm.nih.gov)
    Richa Agarwala (agarwala@ncbi.nlm.nih.gov)

## Citation

Alexandre Souvorov, Richa Agarwala and David J. Lipman. 
**SKESA: strategic k-mer extension for scrupulous assemblies.**
*Genome Biology* 2018 **19**:153.
[doi.org/10.1186/s13059-018-1540-z](https://doi.org/10.1186/s13059-018-1540-z)

Alexandre Souvorov, Richa Agarwala.
**SAUTE: sequence assembly using target enrichment.**
BMC Bioinformatics **22**, 375 (2021). 
[doi.org/10.1186/s12859-021-04174-9](https://doi.org/10.1186/s12859-021-04174-9)

## Introduction
This repositiry is the home for two different short read assemblers and helper programs. These programs are based on de Bruijn graphs and share k-mer related code.

[SKESA](#skesa---strategic-k-mer-extension-for-scrupulous-assemblies) is a de-novo sequence read assembler for microbial genomes. It uses conservative heuristics and is designed to create breaks at repeat regions in the genome. This leads to excellent sequence quality without significantly compromising contiguity. If desired, SKESA contigs could be connected into a GFA graph using [GFA connector](#gfa-connector---create-gfa-graph-starting-from-a-set-of-contigs).

[SAUTE](#saute---sequence-assembly-using-target-enrichment) is a de Bruijn graph based target enriched de-novo assembler designed for assembling genomic and RNA-seq reads sequenced using Illumina. The result is reported as a [GFA graph]( http://gfa-spec.github.io/GFA-spec/) and two nucleotide fasta sequence files for assemblies in the graph.

## Contents
- [Compilation](#compilation)
- [SKESA - Strategic K-mer Extension for Scrupulous Assemblies](#skesa---strategic-k-mer-extension-for-scrupulous-assemblies)
- [SAUTE - Sequence assembly using target enrichment](#saute---sequence-assembly-using-target-enrichment)
- [GFA connector - create GFA graph starting from a set of contigs](#gfa-connector---create-gfa-graph-starting-from-a-set-of-contigs)
- [Kmercounter - effectively count k-mers in a read set](#kmercounter---effectively-count-k-mers-in-a-read-set)

## Compilation

Download current source code

       $ git clone https://github.com/ncbi/SKESA

Alternatively, download last stable release from https://github.com/ncbi/SKESA/releases. Releases also include test data and precompiled binaries. Test data is available in example subdirectory that has commands in file run.test for generating assemblies using the test data.

To build all executables, do following:

       $ cd SKESA


If you would like to build NGS library for accessing reads from SRA, then do

       $ make

Otherwise, if reading inputs only from files, do

       $ make -f Makefile.nongs

BOOST install is expected by makefiles in the SKESA release. If you do not have BOOST on the system path, please specify BOOST_PATH using a command like

       setenv BOOST_PATH /netopt/ncbi_tools64/boost-1.62.0-ncbi1

before running make.

These make files have been tested with BOOST v 1.72.0 and gcc v 7.3.

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
      --estimated_kmers arg (=100)  Estimated number of distinct kmers for bloom 
                                    filter (M, only for hash counter) [integer]
      --skip_bloom_filter           Don't do bloom filter; use --estimated_kmers as
                                    the hash table size (only for hash counter) 
                                    [flag]
    
    Input/output options: at least one input providing reads for assembly must be specified:
      --reads arg                   Input fasta/fastq file(s) for reads (could be 
                                    used multiple times for different runs, could 
                                    be gzipped) [string]
      --use_paired_ends             Indicates that single (not comma separated) 
                                    fasta/fastq files contain paired reads [flag]
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
      --vector_percent arg (=0.05)  Percentage of reads containing 19-mer for the
                                    19-mer to be considered a vector (1. disables) [float (0,1]]
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
SKESA is a de-novo sequence read assembler for microbial genomes based on DeBruijn graphs. It uses conservative heuristics and is designed to create breaks at repeat 
regions in the genome. This leads to excellent sequence quality. Using k-mers longer than mate length and up to insert size also allows SKESA to attain good contiguity 
as determined by the N50 statistic. It is a multi-threaded application that scales well with the number of processors. For different runs with the same inputs, 
including the order of reads, the order and orientation of contigs in the output is deterministic. 

SKESA can process read information by accessing reads from SRA (option `--sra_run`) or from files (option `--reads`) format. Any combination of input streams is 
allowed. Files could be in fasta or fastq format and gzipped, which is recognized automatically. 
    
When accessing reads from SRA, SKESA automatically determines if the read set consists of paired-end or single-end reads. For fasta/fastq input of paired reads 
with separate files for each mate, filenames separated by a comma for first mate followed by the second mate are listed and in this case, the order of reads 
is expected to be same in files for both mates. Alternatively, a single file with both mates could be specified. In this case the reads are expected to be 
interleaved with first mate followed by the second, and the option `--use_paired_ends` must be used.
    
A limitation of the current release is that in case multiple streams of paired reads are provided, it is assumed that all streams have the same insert size. 
User can explicitly specify expected insert size for the reads (option --insert_size). Otherwise, a sample of input reads is used to estimate the expected 
insert size. This sampling may lead to very small differences in assembly of the same read set if the order of reads is different and selected sample gives 
a difference in expected insert size.
    
Two additional options users may wish to specify depending on the resources available to them are as follows:
1. the number of cores (option `--cores`) and
2. total amount of memory in Gb (option `--memory`)

Remaining options are for debugging or modifying algorithm parameters.

Output of assembly is contigs in fasta format. The definition line for contig has format Contig_<N>_<cnt> where <N> is consecutive integers starting from one 
for numbering the contigs and <cnt> is the average count of k-mers in the contig using minimal k-mer length used in the assembly. Contigs are ordered lexicographically.

Limitations:
1. SKESA is designed for haploid genomes. If it is used for diploid genomes or RNAseq reads, it should create breaks at all heterozygous sites in the genome 
and sites for alternative splicing, respectively. The allow_snps option can be used to make some joins at well separated heterozygous sites.
2. SKESA is designed for ILLUMINA reads that do not have systematic homopolymer errors. The assembly for reads that do not have properties similar to ILLUMINA 
is likely to be quite fragmented.
3. Forward-reverse orientation for paired reads is assumed at this time. If this is not true, steps using paired reads are unlikely to change/improve the assembly.
4. Requesting expected insert size to be estimated using a sample is guaranteed to give the same result, including the order of contigs, for the same order 
of reads but may give very small differences if read order is changed and insert size estimate is different.


## SKESA usage examples

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
       
# SAUTE - Sequence Assembly Using Target Enrichment
Version 1.3.0
## Short description
SAUTE is a de Bruijn graph based target enriched de-novo assembler designed for assembling genomic and
RNA-seq reads sequenced using Illumina. It aligns user-provided reference sequences in the target to the de Bruijn graph
for reads and finds all paths with good similarity to each reference. Assembled paths are filtered by SAUTE using full length reads and read pairs to remove possibly chimeric connections. The filtering steps could be partially or completely disabled using options `--no_filter_by_reads` and `--no_filter_by_pairs`.

Similar to [SKESA](#skesa-usage-examples), the user can provide reads using options `--sra_run` or `--reads`. The fasta file for reference sequences are provided using the option `--targets`. The results are reported as a [GFA graph]( http://gfa-spec.github.io/GFA-spec/) (optin `--gfa`) and fasta files (options `--all_variants` and `--selected_variants`). A GFA graph is a convenient way for viewing and analyzing all variants using [Bandage](https://rrwick.github.io/Bandage/).

SAUTE has two flavors - one uses nucleotide reference sequences and the other protein sequences. In both cases, the output is nucleotide assemblies. Depending on the type of reference sequences, do one of the following:
- Nucleotide references
  saute --reads mates1.fa,mates2.fa --targets nuc_references.fa --gfa graph.gfa --all_variants assembly.fa
- Protein references
  saute_prot --reads mates1.fa,mates2.fa --targets prot_references.fa --gfa graph.gfa --all_variants assembly.fa --genetic_code 1
  The genetic code is one of the [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The protein version translates on the fly the segments of de Bruijn graph to the protein space. Any indel changing the frame will corrupt the translated protein and will stop the assembly. There is an option `--allow_frameshifts` which will resume the correct translation after such an indel. It could be useful, for example, for genes with ribosomal slippage.  

## Command line options
### Options shared by saute and saute_prot
    General options:
    -h [ --help ]                    Produce help message
    -v [ --version ]                 Print version
    --cores (=0)                     Number of cores to use (default all) [integer]
    --estimated_kmers (=1000)        Estimated number of distinct kmers for bloom filter (millions) [integer]
                                     To avoid expensive rehashing of the de Bruijn graph, SAUTE uses a light-weight counting Bloom filter to find 
                                     the number of different kmers satisfying the --min_count criterion. Parameter --estimated_kmers provides the 
                                     number (in millions) of different kmers with any count present in the reads. Overestimating will result in a 
                                     small waste of memory. Underestimating will trigger iterative Bloom filter recalculation until a good estimate
                                     is found.

    Input/output options : target file, at least one input for reads and output for gfa must be specified:
    --gfa                            Output file for GFA graph [string]
    --all_variants                   Output file for sequences with all variant combinations [string]
    --max_variants (=1000)           Restricts the number of variants reported in --all_variants per graph [integer]
    --selected_variants              Output file for selected sequences representing all graph elements [string]
    --targets                        Input file with reference sequences [string]
    --sra_run                        Input sra run accession (could be used multiple times for different runs) [string]
    --reads                          Input fasta/fastq file(s) for reads (could be used multiple times for different runs, could be gzipped) [string]
    --use_paired_ends                Indicates that single (not comma separated) fasta/fastq files contain paired reads [flag]

    Assembly options:
    --vector_percent arg (=0.05)  Percentage of reads containing 19-mer for the 19-mer to be considered a vector (1. disables) [float (0,1]]
                                     Before doing any assembly SAUTE removes vectors from the reads. It finds and clips read after all 19-mers found in a relatively high number of reads.

    --min_count (=2)                 Minimal count for kmers retained in graph [integer]
                                     This option keeps noise kmers out of the de Bruijn graphs used for assembling. By default, a kmer has to be seen in 
                                     at least 2 reads to be included. It is allowed to use --min_count 1, and it could help in assembling low 
                                     coverage spots. The drawback of this setting is that all read errors satisfying --fraction criterion will be 
                                     recognized as valid SNPs. To alleviate this, the program will remove all sections supported by only 1 read 
                                     if they were not necessary for a connection. 
    --fraction (=0.05)               Maximum noise to signal ratio acceptable for extension [float [0,1)]
                                     All fork branches with relative abundance less than this are ignored.
    --kmer                           Primary kmer length for assembly (default automatic) [integer]
    --secondary_kmer                 Shorter kmer length for low coverage spots (default automatic) [integer]
                                     For assembling SAUTE uses two de Bruijn graphs with different kmer lengths. It uses the longer kmer most of the 
                                     time and switches to the shorter kmer to assemble the low coverage spots. Both kmers must be odd numbers. In case 
                                     of protein references, both kmers must be divisible by 3. By default, the program will use the closest valid values
                                     not exceeding half and a fifth of the read length.
    --secondary_kmer_threshold (=1)  Coverage threshold for using shorter kmer [integer]
                                     This parameter defines what is considered to be low coverage for the secondary kmer use. With default, these are spots with no 
                                     possible extension or extension supported by only 1 read (could happen only with --min-count 1). Specifying 
                                     a higher number may help with detecting some low coverage forks. With more work shifted toward the shorter 
                                     kmer, there is a higher chance that the program will enter a highly repetitive area and will create a very 
                                     complex graph.

    --word                           Word size for seeds [integer <= 16]
                                     To start assembling, SAUTE scans the main de Bruijn graph and finds seed kmers with high similarity to the reference.
                                     These kmers are used as starting points for the assembly. To trigger a kmer comparison to the reference, one of its ends
                                     has to have an exact match of --word length. For a protein reference, the word must be divisible by 3. Defaults are 8 for 
                                     nucleotide references and 12 for protein references. Latter translates to 4 aa match.

    --kmer_complexity (=2000)        Hard mask reference areas with high number of variant seed possibilities (0 disables masking) [integer]
                                     Repetitive regions could result in very complex output graph and excessive calculation time. All reference areas
                                     for which SAUTE found in each position more different seeds than this parameter are internally hard masked and
                                     will not be assembled. 
    --max_fork_density (=0.1)        Maximal fork density averaged over --buf_length before abandoning assembling (0 disables) [float]
    --buf_length (=200)              Buffer length for fork density [integer]
                                     If assembling enters an unmasked repetitive area, the number of forks and paths to follow may become very high. 
                                     SAUTE calculates the average fork density for the last --buf_length long portion of the assembly. If the fork 
                                     density is above the --max_fork_density threshold, assembling for the current seed is stopped and the program 
                                     moves to the next seed.

    --target_coverage (=0.5)         Keep a path in output if it has alignment to the reference that is at least target_coverage*(reference length) long [float (0,1]]
    --min_hit_len                    If a path is shorter than target_coverage*(reference length), use this length threshold for keeping paths (optional) [integer]

    --extend_ends                    Unambiguously extend graph ends using de-novo assembly [flag]
    --protect_reference_ends         Near the reference ends, don't check if reads support some minimal extension of the fork's branches [flag]
                                     Each fork is analysed using the algorithm described in the SKESA paper, and less supported branches are deleted. 
                                     There is one check which could be detrimental to assembling RNAseq -- each valid branch has to be extendable by ~100bp.
                                     If this is close to a transcript 5' or 3' end, extension is expected to fail. This option disables the extension
                                     check at the reference ends.
    --keep_subgraphs                 Don't remove redundant subgraphs [flag]
                                     If the reference file contains similar references, SAUTE may create multiple graphs which are either identical to
                                     or are sub-graphs of some other graphs. By default, the program will remove redundant
                                     graphs and sub-graphs. This option disables sub-graph removal.
    --use_ambiguous_na               Use ambiguous nucleotide codes for SNPs [flag]
                                     If the option is used, SAUTE will replace SNPs with ambiguous nucleotide codes. This will simplify the graph
                                     and will reduce the number of sequence variants. The counts for individual SNPs will be lost.

    --gap_open                       Penalty for gap opening [integer] (default of 5 for nucleotide references and 11 for protein references)
    --gap_extend                     Penalty for gap extension [integer] (default of 3 for nucleotide references and 1 for protein references)
    --drop_off (=300)                Maximal decrease of score before abandoning assembly path [integer]
                                     SAUTE uses a Smith-Waterman type algorithm for aligning references to the de Bruijn graph. For nucleotide references,
                                     the usual match bonus and mismatch penalty are used. They could be changed from the input (see below Options specific 
                                     for nucleotide targets). For protein references BLOSUM62 substitution matrix is used. In both cases, penalties for gaps 
                                     are specified using --gap_open and --gap_extend parameters. The --drop_off parameter determines how far from the 
                                     reference the assembly may go before the extension is stopped, and the assembly is trimmed to the aligned part.
                                     Larger gap open penalty is used for opening frameshifts (see below).

    Graph cleaning options:
    --no_filter_by_reads             Don't use full length reads for variants filtering [flag]
    --no_filter_by_pairs             Don't use mate pairs for variants filtering [flag]
                                     After assembling, SAUTE will use full length reads and read pairs to remove chimeric connections from the graph. 
                                     Above options will partially or completely disable these steps.
    --max_path (=1000)               Maximal number of path extensions allowed for a single filtering check [integer]
                                     For each fork in graph, SAUTE expands all sequences starting at this fork using either the read length or the insert size.
                                     If the number of expanded sequences exceed this parameter value, the extension length is reduced.
    --not_aligned_len (=10)          Not aligned read length for break count [integer]
    --not_aligned_count (=3)         Number of not aligned reads to make a break [integer]
    --aligned_count (=2)             Number of aligned reads to confirm a connection [integer]
                                     For each expanded sequence, SAUTE counts how many reads or pairs confirm this path and how many contradict it. Based on
                                     these counts, the path is either kept or discarded.
 
    --remove_homopolymer_indels     Remove homopolymer indels [flag]
    --homopolymer_len (=4)          Minimal length of homopolymer [integer]
    --homopolymer_ratio (=0.33)     Coverage ratio threshold for removing homopolymer indels [float]


### Options specific for nucleotide targets    
    --match (=2)                    Bonus for match [integer]
    --mismatch (=3)                 Penalty for mismatch [integer]

### Options specific for protein targets                                         
    --genetic_code                   Genetic code [integer]
                                     Accepted genetic codes:
                                     1 Standard Code
                                     2 Vertebrate Mitochondrial Code
                                     3 Yeast Mitochondrial Code
                                     4 Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
                                     5 Invertebrate Mitochondrial Code
                                     6 Ciliate, Dasycladacean and Hexamita Nuclear Code
                                     9 Echinoderm and Flatworm Mitochondrial Code
                                    10 Euplotid Nuclear Code
                                    11 Bacterial, Archaeal and Plant Plastid Code
                                    12 Alternative Yeast Nuclear Code
                                    13 Ascidian Mitochondrial Code
                                    14 Alternative Flatworm Mitochondrial Code
                                    16 Chlorophycean Mitochondrial Code
                                    21 Trematode Mitochondrial Code
                                    22 Scenedesmus obliquus Mitochondrial Code
                                    23 Thraustochytrium Mitochondrial Code
                                    24 Pterobranchia Mitochondrial Code
                                    25 Candidate Division SR1 and Gracilibacteria Code
                                    26 Pachysolen tannophilus Nuclear Code
                                    27 Karyorelict Nuclear Code
                                    28 Condylostoma Nuclear Code
                                    29 Mesodinium Nuclear Code
                                    30 Peritrich Nuclear Code
                                    31 Blastocrithidia Nuclear Code
                                    33 Cephalodiscidae Mitochondrial UAA-Tyr Code

    --allow_frameshifts              Allow frameshifts [flag]
    --fs_open (=25)                  Penalty for opening a gap that results in a frameshift [integer]
    --na_targets                     References are in-frame nucleotide sequences that will be translated before use [flag]


# GFA connector - create GFA graph starting from a set of contigs
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
      --estimated_kmers arg (=100)         Estimated number of distinct kmers for 
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
      --use_paired_ends                    Indicates that single (not comma 
                                           separated) fasta/fastq files contain 
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
      --vector_percent arg (=0.05)         Percentage of reads containing 19-mer for the
                                           19-mer to be considered a vector (1. disables) [float (0,1]]
      --fraction arg (=0.1)                Threshold for extension
      --entropy arg (=0.51)                Minimal entropy for a seed kmer [float]
      --ext_len arg (=2000)                Maximal length for extension

    Graph cleaning options:
      --not_aligned_len arg (=10)          Not aligned read length for break count
      --not_aligned_count arg (=3)         Number of not aligned reads to make a 
                                           break
      --aligned_count arg (=2)             Number of aligned reads to confirm a 
                                           connection
      --max_path (=1000)                   Maximal number of path extensions allowed for a single filtering check [integer]
      --no_filter_by_reads                 Don't use full length reads for variants
                                           filtering [flag]
      --no_filter_by_pairs                 Don't use mate pairs for variants 
                                           filtering [flag]

    Note that --sra_run option is not available if gfa_connector is built using Makefile.nongs

## Short description
GFA_connector is a SKESA companion program that takes contigs as input, finds pairs of contigs connected by reads or read pairs,
and produces a [GFA graph]( http://gfa-spec.github.io/GFA-spec/). The graph can be viewed and analyzed in [Bandage](https://rrwick.github.io/Bandage/).

GFA_connector constructs a de Bruijn graph (option `--kmer`) and starts extending from each end of the input contigs. It continues each 
possible path up to the maximal length for extension (option `--ext_len`) or until it is connected to another contig. In the 
latter case, two contigs are connected by an edge. At the end of process, all paths not connected to contigs are discarded. This process creates a 
number of graphs each connecting two or more contigs. 

Before the input contigs and connectors are combined into a single GFA graph for output, the connector graphs are analyzed using the reads. First the input 
reads are trimmed to the longest stretches compatible with the de Bruijn graph, and then they are aligned to the connectors. Chimeric connections which do 
not agree with read sequences or paired ends are removed from the graphs. This step could be partially or completely disabled using 
options `--no_filter_by_reads` and `--no_filter_by_pairs`.

To run GFA_connector, one has to provide reads using options `--sra_run` or `--reads` (see SKESA description) and contigs using the option `--contigs`.

## Usage example
The following commands will produce a GFA graph representing a Klebsiella pneumoniae KPR0928 genome (NZ_CP008831.1), two plasmids (NZ_CP008832.1 
and NZ_CP008833.1) and a few contigs which the program was not able to connect
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
  --min_count arg (=2)         Minimal count for kmers retained [integer]
  --vector_percent arg (=0.05) Percentage of reads containing 19-mer for the
                               19-mer to be considered a vector (1. disables) [float (0,1]]
  --estimated_kmers arg (=100) Estimated number of distinct kmers for bloom 
                               filter (millions) for hash count [integer]
  --skip_bloom_filter          Don't do bloom filter; use --estimated_kmers as 
                               the hash table size for hash count [flag]
  --dbg_out arg                de Bruijn graph output
  --text_out arg               Text kmer output
  --hist arg                   File for histogram [string]
  --cores arg (=0)             Number of cores to use (default all) [integer]
  
     Note that --sra_run option is not available if kmercounter is built using Makefile.nongs
```
## Short description
Kmercounter creates and outputs a de Bruijn graph (option `--dbg_out`) for a specified kmer size (option `--kmer`). This file could be used as an optional 
input for GFA_connector in which case the kmer counting step will be skipped. 

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

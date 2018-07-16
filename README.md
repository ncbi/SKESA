# SKESA - Strategic Kmer Extension for Scrupulous Assemblies
Version 2.2

For questions regarding SKESA, please contact
    Alexandre Souvorov (souvorov@ncbi.nlm.nih.gov)
    Richa Agarwala (agarwala@ncbi.nlm.nih.gov)

======================================================================

I.  COMPILATION

    Download current source code for SKESA
       $ git clone https://github.com/ncbi/SKESA
    Alternatively, download last stable release from https://github.com/ncbi/SKESA/releases
    Releases also include test data and precompiled binary.

    Do following:
       $ cd SKESA

    If you would like to build NGS library for accessing reads from SRA,
    then do
       $ make
    Otherwise, if reading inputs only from files, do
       $ make -f Makefile.nongs

    BOOST install is expected by makefiles in the SKESA release. If you
    do not have BOOST on the system path, please specify BOOST_PATH using a command like
       setenv BOOST_PATH /netopt/ncbi_tools64/boost-1.62.0-ncbi1
    before running make.

    These make files have been tested with BOOST v 1.62.0 and gcc v 4.9.

======================================================================

II. SYNOPSIS

    Running
       skesa 
    or 
       skesa -h
    or 
       skesa --help
    gives information about options and produces the following:

    --------------------------------------------------------------
    
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
    
    Input/output options : at least one input providing reads for assembly must be specified:
      --fasta arg                   Input fasta file(s) (could be used multiple 
                                    times for different runs) [string]
      --fastq arg                   Input fastq file(s) (could be used multiple 
                                    times for different runs) [string]
      --gz                          Input fasta/fastq files are gzipped [flag]
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
      --use_paired_ends             Use pairing information from paired reads in 
                                    input [flag]
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
      --seeds arg                   Input file with seeds [string]
      --all arg                     Output fasta for each iteration [string]
      --dbg_out arg                 Output kmer file [string]
      --hist arg                    File for histogram [string]
      --connected_reads arg         File for connected paired reads [string]
    
    --------------------------------------------------------------

    Note that --sra_run option is not available if SKESA is built using Makefile.nongs
    
======================================================================

II. SHORT DESCRIPTION
    
    SKESA is a de-novo sequence read assembler for microbial genomes
    based on DeBruijn graphs. It uses conservative heuristics and is designed to
    create breaks at repeat regions in the genome. This leads to excellent sequence
    quality. Using kmers longer than mate length and up to insert size also allows
    SKESA to attain good contiguity as determined by the N50 statistic. It is a multi-threaded
    application that scales well with the number of processors. For different runs
    with the same inputs, including the order of reads, the order and orientation
    of contigs in the output is deterministic. 

    SKESA can process read information by accessing reads from SRA (option --sra_run)
    or from files in fasta (option --fasta) or fastq (option --fastq) format. Any
    combination of input streams is allowed. For paired reads, if a single file is
    specified, reads are expected to be interleaved with first mate followed by the
    second. To specify a separate file for each mate, filenames separated by a comma
    for first mate followed by the second mate are listed and in this case, the order
    of reads is expected to be same in files for both mates. A limitation of the current
    release is that in case multiple streams of paired reads are provided, it is assumed
    that all streams have the same insert size. User can explicitly specify expected
    insert size for the reads (option --insert_size). Otherwise, a sample of input
    reads is used to estimate the expected insert size. This sampling may lead to very
    small differences in assembly of the same read set if the order of reads is different
    and selected sample gives a difference in expected insert size.

    Two additional options users may wish to specify depending on the resources
    available to them are as follows:
        1. the number of cores (option --cores) and
        2. total amount of memory in Gb (option --memory)

    Remaining options are for debugging or modifying algorithm parameters.

    Output of assembly is contigs in fasta format. The definition line for contig has
    format Contig_<N>_<cnt> where <N> is consecutive integers starting from one for
    numbering the contigs and <cnt> is the average count of kmers in the contig using
    minimal kmer length used in the assembly. Contigs are ordered lexicographically.

    Limitations:

     1. SKESA is designed for haploid genomes. If it is used for diploid genomes or
        RNAseq reads, it should create breaks at all heterozygous sites in the genome
        and sites for alternative splicing, respectively. The allow_snps option can
        be used to make some joins at well separated heterozygous sites.

     2. SKESA is designed for ILLUMINA reads that do not have systematic homopolymer
        errors. The assembly for reads that do not have properties similar to ILLUMINA
        is likely to be quite fragmented.

     3. Forward-reverse orientation for paired reads is assumed at this time. If this
        is not true, steps using paired reads are unlikely to change/improve the assembly.

     4. Requesting expected insert size to be estimated using a sample is guaranteed
        to give the same result, including the order of contigs, for the same order of
        reads but may give very small differences if read order is changed and insert size
        estimate is different.

======================================================================

III. USAGE EXAMPLES

     In all the examples below, we are providing 4 cores and have 48 Gb of memory.

     Example of an assembly that directly accesses SRA for an unpaired read set SRR867211 is:

       $ skesa --sra_run SRR867211 --cores 4 --memory 48 > SRR867211.skesa.fa

     Example of an assembly that directly accesses SRA for a paired read set SRR1960353 is:

       $ skesa --sra_run SRR1960353 --cores 4 --memory 48 --use_paired_ends > SRR1960353.skesa.fa

     Example of an assembly that uses separate fastq files for each mate of SRR1703350 is:

       $ skesa --fastq SRR1703350_1.fq,SRR1703350_2.fq --cores 4 --memory 48 --use_paired_ends > SRR1703350.skesa.fa

     Example of an assembly that uses interleaved mates for SRR1703350 as fastq input is:

       $ skesa --fastq SRR1703350.fq --cores 4 --memory 48 --use_paired_ends > SRR1703350.skesa.fa

     Example of an assembly that uses reads from SRA for SRR1695624 and gzipped fasta for SRR1745628 is:

       $ skesa --sra_run SRR1695624 --fasta SRR1745628.fa.gz --gz --cores 4 --memory 48 --use_paired_ends > SAMN03218571.skesa.fa

     Example of the same assembly as above done with both runs accessed from SRA is:

       $ skesa --sra_run SRR1695624 --sra_run SRR1745628 --cores 4 --memory 48 --use_paired_ends > SAMN03218571.skesa.fa

======================================================================

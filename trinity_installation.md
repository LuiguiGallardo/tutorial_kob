### Installation of Miniconda

To check if you have onda installed, try one of the following your terminal: 

```bash
 which conda
 # or
 conda --help
```

If you get an error, like "command not found", you must to install Miniconda (light version) or Anaconda (complete version). I recommend to you install Miniconda. You could find the instructions in this link: https://docs.conda.io/en/latest/miniconda.html.

### Installation and usage of Trinity

This is the software to assembly *de novo* transcriptomes from RNA-seq data. It depends on multiple programs, that could have been change over time, so the best way to install it is with a new enviroment of conda. To accomplish this you must run something like this:

```bash
conda create -n trinity -c bioconda trinity
```

Where **create** is the option to generate a new environment; **-n** is the name that you want to put to you environment; **-c** is the channel to search for packages; lastly, we call to the package we want to install.

After, we must to call to this new environment:

```bash
conda activate trinity
Trinity
```

You should see something like this:

```bash
###############################################################################
#

     ______  ____   ____  ____   ____  ______  __ __
    |      ||    \ |    ||    \ |    ||      ||  |  |
    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
      |__|  |__|\_||____||__|__||____|  |__|  |____/

    Trinity-v2.8.5


#
#
# Required:
#
#  --seqType <string>      :type of reads: ('fa' or 'fq')
#
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provided in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
#
#  Or,
#      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.
#
####################################
##  Misc:  #########################
#
#  --include_supertranscripts      :yield supertranscripts fasta and gtf files as outputs.
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#
#  --CPU <int>                     :number of CPUs to use, default: 2
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=200)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#                                   (** note: experimental parameter **, this functionality continues to be under development)
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format
#                                   for reads).
#                                   (note: jaccard_clip is an expensive
#                                   operation, so avoid using it unless
#                                   necessary due to finding excessive fusion
#                                   transcripts w/o it.)
#
#  --trimmomatic                   :run Trimmomatic to quality trim reads
#                                        see '--quality_trimming_params' under full usage info for tailored settings.
#                                  
#
#  --no_normalize_reads            :Do *not* run in silico normalization of reads. Defaults to max. read coverage of 200.
#                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
#                                       (note, as of Sept 21, 2016, normalization is on by default)
#     
#  --no_distributed_trinity_exec   :do not run Trinity phase 2 (assembly of partitioned reads), and stop after generating command list.
#
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/home/luigui/trinity_out_dir" 
#                                    note: must include 'trinity' in the name as a safety precaution! )
#             
#  --workdir <string>              :where Trinity phase-2 assembly computation takes place (defaults to --output setting).
#                                  (can set this to a node-local drive or RAM disk)     
#  
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (Trinity-v2.8.5) and exits.
#
#  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
#
#
###############################################################################
#
#  *Note, a typical Trinity command might be:
#
#        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
#
#            (if you have multiple samples, use --samples_file ... see above for details)
#
#    and for Genome-guided Trinity, provide a coordinate-sorted bam:
#
#        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
#
#     see: /home/luigui/miniconda3/envs/trinity/opt/trinity-2.8.5/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#
#     For more details, visit: http://trinityrnaseq.github.io
#
###############################################################################

```

If you have this output, the installation of Trinity in a conda environment is complete! 

To exit from the environment, you must run this line:

```bash
conda deactivate
```

### Optional: create shortcuts to the environment and utilities

The next lines are not necessary, but could help to make things easier in the analysis. To create a shortcut to the conda environment:

```bash
echo "alias trinity='conda activate trinity'" >> ~/.bash_aliases
echo "alias deactivate-conda'conda deactivate'" >> ~/.bash_aliases
```

Now you could get into and exit from the environment simply tipping:

```bash
# To get into:
trinity
# To exit:
deactivate-conda
```

To create a shotcut to the Trinity utilities, you mus to know where Trinity was installed, for example in my computer the path is:

```bash
/home/luigui/miniconda3/envs/trinity/opt/trinity-2.8.5
```

You could find this information when you run **Trinity** (see the ouput above, almost at the end). After, run this line:

```bash
echo "export TRINITY_HOME=/home/luigui/miniconda3/envs/trinity/opt/trinity-2.8.5"
```

And that's it!


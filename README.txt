This is NOT a complete (or finished) workflow intended to provide gene functional annotations for lncRNAs.

This is the directory structure and basic info about the setup.


/xdisk/fionamcc/amcooksey/lnc_nf_proj
(puma) -bash-4.2$ ls -l
total 88
-rw-r--r--.  1 amcooksey shaneburgess  2905 May  4  2023 lnc.config
-rw-r--r--.  1 amcooksey shaneburgess 10274 May  4  2023 lnc_pipeline.nf
-rwx--x--x.  1 amcooksey shaneburgess 14624 May  4  2023 nextflow
-rw-r--r--.  1 amcooksey shaneburgess  4721 May  4  2023 nextflow.config
drwxr-sr-x.  4 amcooksey shaneburgess  5632 May  4  2023 project
drwxr-sr-x. 12 amcooksey shaneburgess  5120 May  4  2023 publish
-rw-r--r--.  1 amcooksey shaneburgess   174 Apr 23 13:41 README.txt
drwxr-sr-x.  2 amcooksey shaneburgess  4096 May  4  2023 ref
drwxr-sr-x.  2 amcooksey shaneburgess  3072 May  4  2023 scripts
-rw-r--r--.  1 amcooksey shaneburgess   258 May  4  2023 slurm_lnc_nf
-rw-r--r--.  1 amcooksey shaneburgess   298 May  4  2023 stats_test.sh
drwxr-sr-x. 82 amcooksey shaneburgess 40960 May  4  2023 work


--lnc.config
This file specifies the containers used for each step of the analyses process along with the slurm params necessary for each

--lnc_pipeline.nf
This fil define the Nextflow process (what comes out of and into each process) and define some basic variables (eg. ref genome).

--nextflow
This is the installed Nextflow executable (could probably use a module going forward??).

--project
This directory contains the input files (FASTQ and subsets made for testing). The 'tools' subdir contains the singularity images used in
the pipeline.

ls -l project/
-rw-r--r-- 1 amcooksey shaneburgess 3012219907 May  4  2023 M20_od1_0.1_1.fq
-rw-r--r-- 1 amcooksey shaneburgess 3012219907 May  4  2023 M20_od1_0.1_2.fq
-rw-r--r-- 1 amcooksey shaneburgess 2078599785 May  4  2023 M20_od2_0.1_1.fq
-rw-r--r-- 1 amcooksey shaneburgess 2078599785 May  4  2023 M20_od2_0.1_2.fq
-rw-r--r-- 1 amcooksey shaneburgess 1639681308 May  4  2023 M6_od1_0.1_1.fq
-rw-r--r-- 1 amcooksey shaneburgess 1639681308 May  4  2023 M6_od1_0.1_2.fq
-rw-r--r-- 1 amcooksey shaneburgess 1580774483 May  4  2023 M6_od2_0.1_1.fq
-rw-r--r-- 1 amcooksey shaneburgess 1580774483 May  4  2023 M6_od2_0.1_2.fq
drwxr-sr-x 2 amcooksey shaneburgess      17408 May  4  2023 not_in_use
-rw-r--r-- 1 amcooksey shaneburgess        656 May  4  2023 slurm_subset_fq
drwxr-sr-x 2 amcooksey shaneburgess       3072 May  4  2023 tools

ls -l tools/
-rwxr-xr-x 1 amcooksey shaneburgess  408743936 May  4  2023 feelnc_0.1.1--pl526_5.sif
-rwxr-xr-x 1 amcooksey shaneburgess  185671680 May  4  2023 filterlncpullgo_1.0.sif
-rwxr-xr-x 1 amcooksey shaneburgess   16990208 May  4  2023 gffread_0.11.7--h8b12597_0.sif
-rwxr-xr-x 1 amcooksey shaneburgess   82722816 May  4  2023 gmap_2019-03-15.sif
-rwxr-xr-x 1 amcooksey shaneburgess    7135232 May  4  2023 star_2.7.3a--0.sif
-rwxr-xr-x 1 amcooksey shaneburgess 1235496960 May  4  2023 trinity_rsem_quantexp.sif

--publish
This is the directory for outputs from each step. The publication dir is specified in lnc_pipeline.nf.  There is one subdir for each process.

--README.txt
This file.

--ref
This directory contains the reference genome, gtf, ref GAF file and ref Brenda tissue ontology files.

ls -l ref
-rw-r--r-- 1 amcooksey shaneburgess         52 May  4  2023 bto.txt
-rw-r--r-- 1 amcooksey shaneburgess 1083128675 May  4  2023 galgal6_nowhitespace.fa
-rw-r--r-- 1 amcooksey shaneburgess      16945 May  4  2023 galgal6_nowhitespace.fa.fai
-rw-r--r-- 1 amcooksey shaneburgess  304159205 May  4  2023 Gallus_gallus.GRCg6a.104.chr.gtf
-rw-r--r-- 1 amcooksey shaneburgess 1083155078 May  4  2023 Gallus_gallus.GRCg6a.dna.toplevel.fa
-rw-r----- 1 amcooksey shaneburgess    1806458 May  4  2023 gga_ens_uni_mapping_15DEC21.csv
-rw-r--r-- 1 amcooksey shaneburgess   35349202 May  4  2023 goa_chicken.gaf

--Scripts
These are some of the original scripts used for executing this type of anaylsis (pre-Nextflow). The Nextflow process is neccesaily 
redundant to a lot of these scripts but these were used for personal reference when setting up the Nextflow pipeline.

ls -l scripts
-rw-r--r-- 1 amcooksey shaneburgess 1688 May  4  2023 slurm_filter_feelnc_topX
-rw-r--r-- 1 amcooksey shaneburgess 2116 May  4  2023 slurm_gtf2fa
-rw-r--r-- 1 amcooksey shaneburgess 1116 May  4  2023 slurm_gtfmerge
-rw-r--r-- 1 amcooksey shaneburgess 1517 May  4  2023 slurm_makeUNILNCmapping
-rw-r--r-- 1 amcooksey shaneburgess  962 May  4  2023 slurm_makeUNILNCmapping_quant
-rw-r--r-- 1 amcooksey shaneburgess 5744 May  4  2023 slurm_pull_go_topX_gafout_feelnc

--slurm_lnc_nf
This is the slurm script to run the Nextflow pipeline on the UA HPC.

--stats_test.sh
This is a script I used to pull the GTF stats for each sample for larger log files.

--work
This is the work dir for Nextflow.

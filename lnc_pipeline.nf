/*Also note that the pipeline work directory is intended to be used as a temporary scratch area. 
*The final workflow outputs are expected to be stored in a different location specified using the publishDir directive.
*/

/*
 *  Parse the input parameters
 */

params.reads = '/xdisk/shaneburgess/amcooksey/lnc_nf_proj/project/*0.25_{1,2}.fq'
params.annot = "$projectDir/ref/Gallus_gallus.GRCg6a.104.chr.gtf"
params.genome = "$projectDir/ref/Gallus_gallus.GRCg6a.dna.toplevel.fa"
params.genome_nowht = "$projectDir/ref/galgal6_nowhitespace.fa"
params.uniprot = "$projectDir/ref/gga_ens_uni_mapping_15DEC21.csv"
params.goa = "$projectDir/ref/goa_chicken.gaf"

genome_file     =  file(params.genome)
genome_nowht	=  file(params.genome_nowht)
annot_file 	=  file(params.annot) 
reads_ch        =  Channel.fromFilePairs(params.reads) 
uniprot_file 	=  file(params.uniprot)
goa_file	=  file(params.goa)

/*
 * Process 1C: Create the genome index file for STAR
 */

process 'STARIndex' {
  input:
      path(genome) from genome_file 

  output:
      path(genome_dir) into genome_dir_ch 

  script:
  """
  mkdir genome_dir 

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}

process 'STAR' {
  input:
      path genome from genome_file 
      path genome_dir from genome_dir_ch 
      path annot from annot_file
      tuple val(replicateId), path(reads) from reads_ch 

  output:
      tuple val(replicateId), path('Aligned.sortedByCoord.out.bam') into aligned_bam_ch 

  script:
  """
  STAR --genomeDir genome_dir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --outSAMtype BAM SortedByCoordinate \
       --sjdbGTFfile $annot \
       --outFilterScoreMinOverLread 0 \
       --outFilterMatchNminOverLread 0 \
       --outFilterMatchNmin 0 \
       --outReadsUnmapped Fastx 

  """
}

process 'TrinityGG' {
  input:
      tuple val(replicateId), path('Aligned.sortedByCoord.out.bam') from aligned_bam_ch
  
  output:
      path('trinity/Trinity-GG.fasta') into trinity_ch

  script:
  """
  Trinity --genome_guided_bam Aligned.sortedByCoord.out.bam \
        --genome_guided_max_intron 10000 \
        --max_memory 3000G \
        --min_kmer_cov 3 \
	--bflyCalculateCPU \
        --output trinity
  """
}

process 'GMAPIndex' {
  input:
      path(genome_nowht) from genome_file

  output:
      path(gmap_index) into gmap_idx_ch

  script:
  """       
  gmap_build \
	-d gmap_indexes \
	-D gmap_index \
	${genome_nowht}
  """
}


process 'GMAP' {
   input:
      path('trinity/Trinity-GG.fasta') from trinity_ch
      path gmap_index from gmap_idx_ch

  output:
      file '*.gff3' into gmap_ch

  script:
  """
  gmap \
	-D gmap_index \
	-d gmap_indexes \
        -f 2 \
	-t 26 \
        trinity/Trinity-GG.fasta > gmap.gff3
  """
}

process 'mapped transcripts GFF to GTF' {

  input:
      path('gmap.gff3') from gmap_ch

  output:
      path('gmap.gtf') into gtf_ch

  script:
  """
      gffread \
      -F \
      -T \
      gmap.gff3 \
      -o gmap.gtf
  """
}

process 'FEELnc_filter' {

  input:
      path('gmap.gtf') from gtf_ch
      path(annot) from annot_file

  output:
      path('filter.log') into filterlog_ch
      path('candidate_lncRNA.gtf') into candidate_ch

  script:
  """
        FEELnc_filter.pl \
	-i gmap.gtf \
	-a $annot \
	-o filter.log \
	> candidate_lncRNA.gtf
  """
}

/*
*process 'make coding annot' {

*   input:
*       path(annot) from annot_file

*   output:
*       path('Gallus_gallus.GRCg6a.104.protein_coding.gtf') into coding_annot_ch

*   script:
*   """
*       grep "protein_coding" $annot \
*	> Gallus_gallus.GRCg6a.104.protein_coding.gtf
*   """
*}

*process 'FEELnc_codpot' {

*  input:
*      file('candidate_lncRNA.gtf') from candidate_ch
*      file(genome) from genome_file
*      file('Gallus_gallus.GRCg6a.104.protein_coding.gtf') from coding_annot_ch
*  output:
*       file('candidate_lncRNA.gtf.lncRNA.gtf') into lncrna_ch
*       file('candidate_lncRNA.gtf.lncRNA.gtf') into lnc_gtf_ch
*       file('candidate_lncRNA.gtf.mRNA.gtf') into mrna_ch

*  script:
*  """
*	export FEELNCPATH=/usr/local/
*        FEELnc_codpot.pl \
*        -i 'candidate_lncRNA.gtf' \
*        -a 'Gallus_gallus.GRCg6a.104.protein_coding.gtf' \
*        -b transcript_biotype=protein_coding \
*        -g $genome \
*	--mode=shuffle
*  """
*}


*process 'FEELnc_classifer' {

*  input:
*      path('*candidate_lncRNA.gtf.lncRNA.gtf') from candidate_ch
       path(proteincoding_"$annot") from coding_annot_ch

*  output:
*      path(class) into feelncclass_ch

*  script:
*  """
* FEELnc_classifier.pl \
*        -i '*candidate_lncRNA.gtf.lncRNA.gtf' \
*        -a path(proteincoding_"$annot") from coding_annot_ch \
*        -b
*  """
*}

*THIS PROCESS NEEDS TO USE AN INPUT FOR A SECOND TIME--HOW?
*If you need to connect a process output channel to more than one process
*or operator use the into operator to create two (or more) copies of the same channel and use each of them to connect a separate process.

*process 'gffcompare' {

*  input:
*      each('*candidate_lncRNA.gtf.lncRNA.gtf') from lnc_gtf_ch

*  output:
*      path(unified_lncrna_ids) into mergedgtf_ch

*  script:
*  """
*    gffcompare \
*        -p UNILNC \
*        -o unified_lncrna_ids \
*        --strict-match \
*        --debug \
*        '*candidate_lncRNA.gtf.lncRNA.gtf'	
*  """
*}

*ADD GTF 2 FASTA FOR CANDIDATE FILES--FOR WHAT--nothing, just in case
*ADD GFFCOMPARE, UNILNC MAPPING
*ADD FILTER FEELNC TARGETS
*ADD PULLGO GAFOUT

*/


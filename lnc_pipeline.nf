/*
 *  Parse the input parameters
 */

params.reads = '/xdisk/shaneburgess/amcooksey/lnc_nf_proj/project/*0.1_{1,2}.fq'
params.annot = "$projectDir/ref/Gallus_gallus.GRCg6a.104.chr.gtf"
params.genome = "$projectDir/ref/Gallus_gallus.GRCg6a.dna.toplevel.fa"

genome_file     =  file(params.genome)
annot_file 	=  file(params.annot) 
reads_ch        =  Channel.fromFilePairs(params.reads) 

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
      path(trinity) into trinity_ch

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
      path(genome) from genome_file

  output:
      path(gmap_dir) into gmap_idx_ch

  script:
  """
  mkdir gmap_dir 
        
  gmap_build \
	-d gmap_index \
	-D gmap_dir \
	${genome}
  """
}

/*
*process 'GMAP' {
*  input:
*      tuple val(replicateId), path('Trinity-GG.fasta') from trinity_ch, path(gmap_dir) from gmap_idx_ch

*  output:
*      tuple val(replicateId), path(gmap_out) into gmap_ch

*  script:
*  """
*   gmap \
*	-d gmap_index \
*        -D gmap_dir \
*        -f 2 \
*	-t 26 \
*        Trinity-GG.fasta
*  """
*}

*process 'mapped transcripts GFF to GTF' {

*  input:
*      tuple val(replicateId), path(gmap_out) from gmap_ch

*  output:
*      tuple val(replicateId), path('gmap.gtf') into gtf_ch

*  script:
*  """
*   gffread \
*      -F \
*      -T \
*      $gmap_out\
*      -o gmap.gtf
*  """
*}

*process 'FEELncfilter' {

*  input:
*      tuple val(replicateId), path('gmap.gtf') from gtf_ch, path(annot) from annot_file

*  output:
*      tuple val(replicateId), path(log) into feelnclog_ch, path('candidate_lncRNA.gtf') into candidate_ch

*  script:
*  """
*        FEELnc_filter.pl \
*	-i gmap.gtf \
*	-a $annot \
*	-o $log \
*	> candidate_lncRNA.gtf
*  """
*}

*process 'FEELnccodpot' {

*  input:
*      tuple val(replicateId), path('candidate_lncRNA.gtf') from candidate_ch, path(genome) from genome_file, path(annot) from annot_file

*  output:
*      tuple val(replicateId), path(log) into feelnclog_ch, path('candidate_lncRNA.gtf.lncRNA.gtf') into candidate_ch

*  script:
*  """
*        FEELnc_codpot.pl \
*        -i 'candidate_lncRNA.gtf' \
*        -a $annot \
*        -b transcript_biotype=protein_coding \
*        -g $genome \
*	--mode=shuffle
*  """
*}

*process 'make coding annot' {

*   input:
*       path(annot) from annot_file

*   output:
*       path(codingannot) into coding_annot_ch

*   script:
*   """
*       grep "protein_coding" $annot
*   """
*}

*process 'FEELncclassifer' {

*  input:
*      tuple val(replicateId), path('candidate_lncRNA.gtf.lncRNA.gtf') from candidate_ch, path(codingannot) from coding_annot_ch

*  output:
*      tuple val(replicateId), path(class) into feelncclass_ch

*  script:
*  """
* FEELnc_classifier.pl \
*        -i candidate_lncRNA.gtf.lncRNA.gtf \
*        -a path(codingannot) from coding_annot_ch \
*        -b
*  """
*}
*/

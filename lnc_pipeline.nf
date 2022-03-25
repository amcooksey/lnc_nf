/*
 *  Parse the input parameters
 */

params.reads = '/xdisk/shaneburgess/amcooksey/lnc_nf_proj/project/*_{1,2}.fq'
params.annot = "$projectDir/ref/Gallus_gallus.GRCg6a.104.chr.gtf"
params.genome = "$projectDir/ref/Gallus_gallus.GRCg6a.dna.toplevel.fa"

genome_file     =  file(params.genome)
annot_file 	= file(params.annot) 
reads_ch        =  Channel.fromFilePairs(params.reads) 

/*
 * Process 1C: Create the genome index file for STAR
 */

process 'STARIndex' {
//  errorStrategy 'ignore'
  input:
      file(genome) from genome_file 

  output:
      file(genome_dir) into genome_dir_ch 

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
      file genome from genome_file 
      file genome_dir from genome_dir_ch 
      file annot from annot_file
      set replicateId, file(reads) from reads_ch 

  output:
      set replicateId, file('Aligned.sortedByCoord.out.bam') into aligned_bam_ch 

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
       --outReadsUnmapped Fastx \
       --outSAMunmapped Within \
       --outWigType wiggle \
       --outWigStrand Stranded \
       --outWigNorm RPM

  """
}



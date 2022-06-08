/*Also note that the pipeline work directory is intended to be used as a temporary scratch area. 
*The final workflow outputs are expected to be stored in a different location specified using the publishDir directive.
*/

/*
 *  Parse the input parameters
 */

params.reads = '/xdisk/shaneburgess/amcooksey/lnc_nf_proj/project/*0.01_{1,2}.fq'
params.annot = "$projectDir/ref/Gallus_gallus.GRCg6a.104.chr.gtf"
params.genome = "$projectDir/ref/Gallus_gallus.GRCg6a.dna.toplevel.fa"
params.genome_nowht = "$projectDir/ref/galgal6_nowhitespace.fa"
params.uniprot = "$projectDir/ref/gga_ens_uni_mapping_15DEC21.csv"
params.goa = "$projectDir/ref/goa_chicken.gaf"
params.bto = "$projectDir/ref/bto.txt"

genome_file     =  file(params.genome)
genome_nowht	=  file(params.genome_nowht)
annot_file 	=  file(params.annot) 
reads_ch        =  Channel.fromFilePairs(params.reads) 
uniprot_file 	=  file(params.uniprot)
goa_file	=  file(params.goa)
bto_file 	=  file(params.bto)


process 'STARIndex' {
  publishDir "$projectDir/publish/starindex", overwrite: true

  errorStrategy 'finish'

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
  publishDir "$projectDir/publish/star/$replicateId", overwrite: true

  errorStrategy 'finish'

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
  publishDir "$projectDir/publish/trinity/${replicateId}/", overwrite: true

  errorStrategy 'retry'
  maxRetries 3
  input:
      tuple val(replicateId), path('Aligned.sortedByCoord.out.bam') from aligned_bam_ch

  output:
      tuple val(replicateId), path('trinity/Trinity-GG.fasta') into trinity_ch

  script:
  """
  Trinity --genome_guided_bam Aligned.sortedByCoord.out.bam \
        --genome_guided_max_intron 10000 \
        --max_memory 3000G \
        --min_kmer_cov 3 \
	--output trinity
  """
}


process 'GMAPIndex' {
  publishDir "$projectDir/publish/gmapindex", overwrite: true

  errorStrategy 'finish'

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
  publishDir "$projectDir/publish/gmap/${replicateId}", overwrite: true

  errorStrategy 'finish'

   input:
      tuple val(replicateId), path('trinity/Trinity-GG.fasta') from trinity_ch
      path gmap_index from gmap_idx_ch

  output:
      tuple val(replicateId), path("${replicateId}_gmap.gff3") into gmap_ch

  script:
  """
  gmap \
	-D gmap_index \
	-d gmap_indexes \
        -f 2 \
	-t 26 \
        trinity/Trinity-GG.fasta > "${replicateId}_gmap.gff3"
  """
}


process 'mapped transcripts GFF to GTF' {
  publishDir "$projectDir/publish/gff2gtf/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("${replicateId}_gmap.gff3") from gmap_ch

  output:
      tuple val(replicateId), path("${replicateId}_gmap.gtf") into gtf_ch

  script:
  """
      gffread \
      -F \
      -T \
      "${replicateId}_gmap.gff3" \
      -o "${replicateId}_gmap.gtf"
  """
}


process 'FEELnc_filter' {
  publishDir "$projectDir/publish/FEELnc/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("${replicateId}_gmap.gtf") from gtf_ch
      path(annot) from annot_file

  output:
      tuple val(replicateId), path("${replicateId}_filter.log") into filterlog_ch
      tuple val(replicateId), path("${replicateId}_candidate_lncRNA.gtf") into candidate_ch

  script:
  """
        FEELnc_filter.pl \
	-i "${replicateId}_gmap.gtf" \
	-a $annot \
	-o "${replicateId}_filter.log" \
	> "${replicateId}_candidate_lncRNA.gtf"
  """
}

process 'make coding annot' {
  publishDir "$projectDir/publish/coding_annot", overwrite: true

  errorStrategy 'finish'

   input:
       path(annot) from annot_file

   output:
       path('Gallus_gallus.GRCg6a.104.protein_coding.gtf') into coding_annot_ch

   script:
   """
       grep "protein_coding" $annot \
	> Gallus_gallus.GRCg6a.104.protein_coding.gtf
   """
}

process 'FEELnc_codpot' {
  publishDir "$projectDir/publish/FEELnc/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("${replicateId}_candidate_lncRNA.gtf") from candidate_ch
      path(genome) from genome_file
      path(annot) from annot_file

  output:
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") into lncrna_ch
       path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") into lnc_gtf_ch
       path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") into mapping_ch
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") into lnc_fa_ch
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.mRNA.gtf") into mrna_ch

  script:
  """
	export FEELNCPATH=/usr/local/
        FEELnc_codpot.pl \
        -i "${replicateId}_candidate_lncRNA.gtf" \
        -a $annot \
        -b transcript_biotype=protein_coding \
        -g $genome \
	--mode=shuffle
  """
}

process 'FEELnc_classifier' {
  publishDir "$projectDir/publish/FEELnc/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") from lncrna_ch
      path('Gallus_gallus.GRCg6a.104.protein_coding.gtf') from coding_annot_ch

  output:
      tuple val(replicateId), path("${replicateId}_lncRNA_classes.txt") into feelncclass_ch
      tuple val(replicateId), path("${replicateId}_candidate_lncRNA.gtf.lncRNA.feelncclassifier.log") into classlog_ch

  script:
  """
 	FEELnc_classifier.pl \
        -i "feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf" \
        -a 'Gallus_gallus.GRCg6a.104.protein_coding.gtf' \
        -b \
	> "${replicateId}_lncRNA_classes.txt"
  """
}


process 'gtf2fa_lnc' {
  publishDir "$projectDir/publish/gffcompare/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input: 
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf") from lnc_fa_ch

  output: 
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.fa") into lnc_fasta_ch

  script:
  """
        gffread \
        -w "feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.fa" \
        -g $genome_nowht \
	"feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.lncRNA.gtf"
  """
}

process 'gtf2fa_mrna' {
  publishDir "$projectDir/publish/gffcompare/${replicateId}", overwrite: true 

  errorStrategy 'finish'

  input: 
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.mRNA.gtf") from mrna_ch

  output: 
       tuple val(replicateId), path("feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.mRNA.fa") into mrna_fasta_ch

  script:
  """
	gffread \
	-w "feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.mRNA.fa" \
        -g $genome_nowht \
	"feelnc_codpot_out/${replicateId}_candidate_lncRNA.gtf.mRNA.gtf"

  """
}


Channel
     .from lnc_gtf_ch
     .toSortedList()
     .set {merge_ch}

process 'gtfmerge' {
  publishDir "$projectDir/publish/gtfmerge/", overwrite: true

  errorStrategy 'finish'

  input:
      path(gtf) from merge_ch

  output:
      path('unified_lncrna_ids.combined.gtf') into mergedgtf_ch
      path('unified_lncrna_ids.loci') into mergeloci_ch
      path('unified_lncrna_ids.tracking') into mergetrack_ch
      path('unified_lncrna_ids.stats') into mergestats_ch

  script:
  """
    gffcompare \
        -p UNILNC \
        -o unified_lncrna_ids \
        --strict-match \
        --debug \
        *gtf	
  """
}


process 'UNILNCmapping' {
  publishDir "$projectDir/publish/gtfmerge", overwrite: true

  errorStrategy 'finish'

  input:
       path('unified_lncrna_ids.stats') from mergestats_ch
       path('unified_lncrna_ids.tracking') from mergetrack_ch

  output:
       path('*_lncRNA_mapping.txt') into mapped_ch
       path('*_lncRNA_mapping.txt') into maptogo_ch

  script:
  """
    makeUNILNCmapping.sh
  """
}

process 'top targets' {
  publishDir "$projectDir/publish/FEELnc/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("${replicateId}_lncRNA_classes.txt") from feelncclass_ch
      path('*_lncRNA_mapping.txt') from mapped_ch

  output:

      tuple val(replicateId), path("lncRNA_classes_UNILNC_top2.txt") into toptarg

  script:
  """
	filter_feelnc_topX.sh
  """
}

process 'pull GO each' {
  publishDir "$projectDir/publish/GO/${replicateId}", overwrite: true

  errorStrategy 'finish'

  input:
      tuple val(replicateId), path("lncRNA_classes_UNILNC_top2.txt") from toptarg
      path(uniprot) from uniprot_file
      path(goa) from goa_file
      path(bto) from bto_file
      path('*_lncRNA_mapping.txt') from maptogo_ch

  output:
      path("${replicateId}_feelnc_GO_2.gaf") into go_ch

  script:
  """
	pull_go_topX_gafout_feelnc_each.sh \
	-a AgBase \
	-t 9031 \
	-o lncRNA \
	-d Agbase \
	-T "${replicateId}"

  """
}


Channel
     .from go_ch
     .toSortedList()
     .set {mergego_ch}


process 'pull GO all' {
  publishDir "$projectDir/publish/GO/", overwrite: true

  errorStrategy 'finish'

  input:
      path('*_feelnc_GO_2.gaf') from mergego_ch

  output:
      path('feelnc_GO_2.gaf') into tot_GO_ch

  script:
  """
        pull_go_topX_gafout_feelnc_all.sh 

  """
}

/*


**ADD PULLGO GAFOUT


*THIS PROCESS NEEDS TO USE AN INPUT FOR A SECOND TIME--HOW?
*If you need to connect a process output channel to more than one process
*or operator use the into operator to create two (or more) copies of the same channel and use each of them to connect a separate process.


*/


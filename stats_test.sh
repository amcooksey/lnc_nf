#!/bin/bash

cd /xdisk/shaneburgess/amcooksey/lnc_nf_proj/publish/gtfmerge

grep 'Summary' unified_lncrna_ids.stats > stattmp.txt
sed -i 's/\#\= Summary for dataset\: //' stattmp.txt
sed -i 's/_candidate_lncRNA.gtf.lncRNA.gtf//' stattmp.txt

readarray -t tissuearray < stattmp.txt

rm stattmp.txt


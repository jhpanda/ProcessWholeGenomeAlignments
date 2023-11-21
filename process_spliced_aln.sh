

suffix="v6"
fgroup="last_orf_unannotated_nohomologs_${suffix}.txt"
seqdir="raw_hit_300_new_${suffix}"

## first extract sequence and perform splice alignments
echo "Started spliced alignment"
sbatch --wait --wrap="python extract_gff_hit.py ${fgroup} ${seqdir}"
echo "Finished spliced alignment"

## then translate spliced DNA sequence
ls ${seqdir}/*genewise623.fa > groups_new.txt
ls ${seqdir}/*spaln.fa >> groups_new.txt
python translate2pro.py
echo "Finished translation"

## align them to reference gene
echo "Started MAFFT-linsi"
sbatch --wait mafft.sh ${seqdir}
echo "Finished MAFFT-linsi"

## extract alignment score and summarize
python extract_alnscore_genewise.py $fgroup $seqdir $suffix
echo "All Done"

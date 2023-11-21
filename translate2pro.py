
from fasta2seq import *
from translate_cds import *

#groupfile = "groups_blastp_E0.05.txt"
groupfile = "groups_new.txt"
#groupfile = "test_new.txt"
#groupfile = "test.txt"

lines = open(groupfile,"r")
for line in lines:
    prefix = line[:-4]
    fdna   = line[:-1]
    fpro   = "%s_pro.fa"%prefix
    print(fpro)

    dnaSeq = fasta2seq(fdna)
    proSeq = dict() 

    for key in dnaSeq:
        proSeq[key] = translate(dnaSeq[key],name=prefix)
    seq2fasta(proSeq,fpro)

lines.close()

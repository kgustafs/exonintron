## The purpose of this python program is to extract information from the Ensembl annotation for the mouse genome mm9.
## Our desired information consists of the boundaries of all introns and exons in the genome.

## INPUT: ens.txt file with columns in the following fields:
## #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames

## OUTPUT: two files, one with intron boundaries matched to ENS gene id and one with exon boundaries also matched
## these files will be similar, except for the distinction between introns and exons.
## output fields: chrom exonStart exonEnd name1 exon# exonLength strand

import sys
import string
import sets

def exonintronlist():
    ens_file = open('../auxfiles/ens.txt','r').readlines()

    exonslist = open('exonslist.txt','w')
    intronslist = open('intronslist.txt','w')

    ## read in the geneselection file from Felix Naef
    geneset = set(line.strip() for line in open('../auxfiles/geneselection.txt'))

    ## read in the set of microarray ids for CycliX
    marray = [line.strip() for line in open('../auxfiles/microarrayid.txt')]

    ## read in the conversion from microarray id to ensembl id
    ens2marray = dict(line.strip().split() for line in open('../auxfiles/conv.txt'))
    ## invert the dictionary
    marray2ens = dict((v,k) for (k,v) in ens2marray.items())

    ## this is meant to be a list of all ensembl values in the conv.txt dictionary
    marrayens = []

    ## locate all microarray entries with an ensembl value and make a list of the ensembl values
    for i in marray:
        if i in marray2ens:
            marrayens.append(marray2ens[i])

    ## place the list of ensembl values in a set
    marrayens_set = set(marrayens)

    ## starting at 2nd element to avoid header line
    ## use string split to break lines into entries separated by spaces - entries counted from zero
    for line in ens_file[1:]:
        words = string.split(line)
        ## check whether the transcript is in the geneselection.txt and microarray.txt sets
        if words[1] in geneset and words[1] in marrayens_set:
            ## pull out exon start and end positions
            ex1 = str(words[9]).rstrip(',')
            ex2 = str(words[10]).rstrip(',')
            exon1 = ex1.split(',')
            exon2 = ex2.split(',')
            num_exon = len(exon1)
            for x in range(0,num_exon):
                ## compute exon length
                len_exon = int(exon2[x])-int(exon1[x])
                exonslist.write(str(words[2])+'\t'+exon1[x]+'\t'+exon2[x]+'\t'+
                                str(words[1])+'\t'+str(x+1)+'\t'+str(len_exon)+'\t'+words[3]+'\n')
                    ## this is necessary since the transcript ends with an exon rather than an intron
                if x < num_exon-1:
                    ## compute intron length
                    len_intron = int(exon1[x+1])-int(exon2[x])
                    intronslist.write(str(words[2])+'\t'+exon2[x]+'\t'+exon1[x+1]+'\t'+
                                      str(words[1])+'\t'+str(x+1)+'\t'+str(len_intron)+'\t'+words[3]+'\n')

exonintronlist()

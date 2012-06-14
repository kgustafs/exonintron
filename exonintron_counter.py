## DOUBLE ## indicates a durable comment, other # are temporary debugs

## This python program uses mapped ChIPseq data converted (elsewhere)
## from fastq to BAM to text format to give the occupancy of
## polII RPB2 units on the genome of mouse, sorted on the introns and exons
## for all genes. It is possible to select various subgroups of exons and introns to
## examine. Extensions of this program include separation between
## constitutive and alternative exons (separate module),
## and rescaling (separate module) of exons or introns to permit
## integration of their signals across the entire genome.

## Designed to take input from exonintron_finder.py, and coordinate with shell
## scripts exon_chr_mm9_merge and intron_chr_mm9_merge to parallelize across chromosomes.

## 'merge' refers to the HiSeq data with duplicates included:
## /archive/systemsx/cyclix/analysis/Mouse_CCNRC/POLR2BHS/bam/

import sys
import string
import numpy as np

def exoncount(chrid,delta,infile,outfile):

    ## arguments: chrid is the chromosome id; delta is the distance from the start and
    ## end of the exon/intron to include in window of reads; infile is the input exon/intron file with start and end positions;
    ## outfile is the ragged-end output file with all the read positions on the exon/intron;

    ## convert delta to integer
    delta = int(delta)

    ## create a file for storing the reads on the exons/introns (name specified by user)
    reads = open(outfile,'w')

    ## read file with a single column of read positions
    ## one for each chromosome allows embarrassing parallelization 
    chr_read_file = chrid+'/'+chrid+'_merge_col4.txt'
    ## strip end-of-lines and split into single read entries
    chr_read = [line.strip().split() for line in open(chr_read_file,'r')]

    ## convert all read strings to integers
    chr_int = [map(int,x) for x in chr_read]

    ## convert to a numpy array for using closestread
    chr_array = np.array(chr_int)

    ## this file must have the exon/intron start in the second position of each row,
    ## followed by the exon/intron end in the third position of each row
    startend = [line.strip() for line in open(infile,'r')]

    ## most important loop - line-by-line through sequencing data to match with exons/introns
    for line in startend:
        ## break the row up into pieces separated by spaces
        elements = string.split(line)
        ## only look at the user-specified chromosome
        if elements[0] == chrid:
            ## find the index of the closest read to delta from the start of the feature
            startidx = closestread(chr_array,np.int(elements[1])-delta)
            ## find the index of the closest read to delta from the end of the feature
            endidx   = closestread(chr_array,np.int(elements[2])+delta)
            ## record the reads within delta of the start and end of the feature, inclusive
            for i in range(startidx,endidx):
                ## get the read from the array at that index
                i_pos = chr_array[i,0]
                ## write the read to tab-spaced file
                reads.write(str(i_pos)+'\t')
                ## make new line for new exon
            reads.write('\n')

def closestread(array,value):
    ## function to locate the nearest value in an array to a set of values
    ## takes a numpy array and a numeric value for arguments
    idx = (np.abs(array-value)).argmin()
    ## returns the index of the entry of the array nearest to value
    return idx

## function chooser for taking arguments
func_arg = {"-a": exoncount}

## specify arguments for that function chosen above
if __name__ == "__main__":
    func_arg[sys.argv[1]](sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) 

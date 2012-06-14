
# script for running exonintron_counter.py in parallel across chromosomes for mouse

# put the chromosome name in the species context into a new file
#awk '{print $2}' ~/archive/chrmapmm9.b2w > file1

# initialize index
index=0

# set the (delta) distance from the end and start of exon/intron to count
delta=200

# put the file of chromosome names into an array
while read line; do
    prefixarray[$index]="$line"
    index=$(($index+1))
done < "./file1"

# create directory and col4 file for each chromosome if necessary, then run the python exon counter with bsub
for (( row=0;row<index;row++)); do
    echo ${prefixarray[row]}
    if [ ! -d "${prefixarray[row]}" ]; then
	mkdir "${prefixarray[row]}"
    fi
    if [ ! -f ${prefixarray[row]}/${prefixarray[row]}_merge_col4.txt ]; then
	# this is valid only on vital-it; NOTE that ZT02 is specified
	awk '{print $4}' ~/archive/rservscr/HiSeq-vader/RPB2_ZT02_R1/merge/txtfiles/${prefixarray[row]}_merge.txt > ${prefixarray[row]}/${prefixarray[row]}_merge_col4.txt
    fi
    bsub -M 16000000 -R "rusage[mem=16000]" "python exonintron_counter.py -a ${prefixarray[row]} $delta exonsbed.txt ${prefixarray[row]}/${prefixarray[row]}_exonreads.txt"
done

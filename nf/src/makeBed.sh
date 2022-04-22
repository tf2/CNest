#!/usr/bin/env bash
chromSize=../data/hg38.chrom.sizes
winSize=1000
blackList=../data/hg38.ENCODE.blacklist.ENCFF356LFX.bed
gapList=../data/hg38.gaps
out=../data/hg38.1kb.baits.bed 

# make genome wide baits
bedtools makewindows -g ${chromSize} -w ${winSize} > tmp1.bed

# remove black list and gap regions
# If the bait has >0 bp in blacklist or gap regions, the whole bait is removed
bedtools subtract -a tmp1.bed -b ${blackList} -A > tmp2.bed
bedtools subtract -a tmp2.bed -b ${gapList} -A > ${out} 

# remove chr
# sed 's/^chr//' hg38.1kb.bins.bed > tmp

# remove tmp files
rm tmp1.bed tmp2.bed
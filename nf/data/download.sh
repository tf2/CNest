# hg38 chrom sizes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg38.chromInfo" | head -25 | tail -n +2 | sort -k1,1 > hg38.chrom.sizes
# hg38 gaps
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, chromStart, chromEnd, type from hg38.gap" | tail -n +2 > hg38.gaps
# hg38 ENCODE blacklist regions
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz \
    -O hg38.ENCODE.blacklist.ENCFF356LFX.bed.gz
gzip -d hg38.ENCODE.blacklist.ENCFF356LFX.bed.gz
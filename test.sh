
scTE_build -g mm10 -te Data/TE.bed -gene Data/Gene.gtf -o Data/test

scTE -i Data/test.bam -p 12 --min_genes 1 -o out --genome mm10 -m exclusive -x  Data/test.exclusive.idx


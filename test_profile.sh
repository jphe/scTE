
python3 -m cProfile -s cumtime bin/scTE_build -g mm10 -te Data/TE.bed -gene Data/Gene.gtf -o Data/test >profile_scTE_build.txt

python3 -m cProfile -s cumtime bin/scTE -i Data/test.bam -p 1 --min_genes 1 -o out --genome mm10 -x Data/test.exclusive.idx >profile_scTE.txt


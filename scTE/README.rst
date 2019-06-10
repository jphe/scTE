
INSTALL:
python setup.py install

Usage:
scTE -i INPUT.bam -p 8 --min_genes 200 --min_counts 1000 --genome mm -m exclusive -te default -gene default  -o out 

Note:
The bam file can output by CellRanger or STARsolo, but it must be sorted by chromosome position 
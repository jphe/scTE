[![PyPI Version](https://img.shields.io/pypi/v/pyGenomeTracks.svg?style=plastic)](https://pypi.org/project/pyGenomeTracks/) [![bioconda-badge](https://img.shields.io/conda/vn/bioconda/pyGenomeTracks.svg?style=plastic)](https://anaconda.org/bioconda/pygenometracks) [![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=plastic)](http://bioconda.github.io)

scTE
==============

Quantifying transposable element expression from single-cell sequencing data
----------------------------------------------------------------------

scTE takes as input:

 * aligned sequence reads (bam/sam format)
 * the genomic location of TEs (gtf format)
 * the genomic location of genes (gtf format)


![scTE workflow](./docs/content/images/hic_example_nat_comm_small.png)


Installation
------------
scTE works with python >=3.6.

```bash
$ git clone https://github.com/jphe/scTE.git
$ cd scTE
$ python setup.py install
```

Usage
-----
scTE builds genome indices for the fast alignment of reads to genes and TEs. These indices can be automatically generated using the commands:

```bash
$ scTE_build -g mm10 # mouse genome
$ scTE_build -g hg38 # human genome
```

These two scripts will automatically download the genome annotations, for mouse:

```bash
$ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.g tf.gz
$ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
```

Or for human:

```bash
$ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf .gz
$ http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
```

`mm10, hg38` is the genome assembly version, scTE currently support mm10 and hg38. 
If you want to use your customs reference, you can use the ` -gene -te` options:

```
scTE_build -te TEs.bed -gene Genes.gtf -o custome.idx

-te Six columns bed file for transposable elements annotation. Need the -gene option. For bed file format, see from [UCSC](https://genome.ucsc.edu/FAQ/FAQformat)
-gene Gtf file for genes annotation. Need the -te option. For gtf file format, see from [UCSC](https://genome.ucsc.edu/FAQ/FAQformat)
```

These annotations are then processed and converted into genome indices. The scTE algorithm will allocate reads first to gene exons, and then to TEs by default. Hence TEs inside exon/UTR regions of genes annotated in GENCODE will only contribute to the gene, and not to the TE score. This feature can be changed by setting ‘–mode/-m exclusive’ in scTE, which will instruct scTE to assign the reads to both TEs and genes if a read comes from a TE inside exon/UTR regions of genes.


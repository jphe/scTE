import os,sys,gzip,time
import numpy as np
from miniglbase import genelist, glload, location
#from glbase3 import genelist, glload, location

chr_list = list(range(1,50)) + ['X','Y'] # delete M, MT

def cleanexon(exons):
    #oh=gzip.open('_tmp/%s.bed.gz'%(genefilename),'wt')
    tmp = []
    for k in sorted(exons):
        E=[]
        for it in exons[k]:
            E+=list(range(it[1],it[2]))
        E=sorted(set(E))

        s=0
        #tmp=[]
        for id in range(0,len(E)-1):
            if E[id+1]-E[id] >1:
                en=id
                tmp.append({'loc': location(chr=it[0], left=E[s], right=E[en]), 'annot': k}) #[E[s],E[en]])
                s=en+1
        tmp.append({'loc': location(chr=it[0], left=E[s], right=E[id+1]), 'annot': k}) # [E[s],E[id+1]])

        #for item in tmp:
        #    oh.write('%s\t%s\t%s\t%s\n'%(it[0],item[0],item[1],k))
    #oh.close()
    return tmp

def genomeIndex(genome, mode, geneurls, teurls):
    print(mode)
    if not os.path.exists('_tmp'):
        os.system('mkdir -p _tmp')

    #os.system('wget -c %s'%geneurls)
    #os.system('wget -c %s'%teurls)

    geneform ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}
    teform ={'force_tsv': True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])', 'annot': 10}

    genefilename = geneurls.split('/')[-1:][0]
    tefilename = teurls.split('/')[-1:][0]
    print(genefilename, tefilename)

    raw = {}
    clean = {}
    if '.gz' in genefilename:
        o = gzip.open(genefilename,'rb')
    else:
        o = open(genefilename,'r')

    for idx, l in enumerate(o):
        if '.gz' in genefilename:
            l=l.decode('ascii')
        if l.startswith('#'):
            continue
        t=l.strip().split('\t')
        if t[2]=='exon' or t[2]=='UTR':
            chr = t[0]
            left = int(t[3])
            riht =  int(t[4])
            name=t[8].split('gene_name "')[1].split('";')[0]

            if name not in raw:
                raw[name] = []
            raw[name].append([chr,left,riht])

            if 'protein_coding' not in l and 'lincRNA' not in l:
                continue
            if name not in clean:
                clean[name] = []
            clean[name].append([chr,left,riht])

        if idx > 1000:
            break
    o.close()

    raw = cleanexon(raw) # '%s.raw' % genefilename, raw)
    clean = cleanexon(clean) # '%s.clean' % genefilename, clean)

    if mode == 'exclusive':
        gene = {}
        #cleangl # o = gzip.open('_tmp/%s.clean.bed.gz'%(genefilename),'rb')
        for l in clean:
            #t = l.decode('ascii').strip().split('\t')
            chr = l['loc'].loc['chr'] # t[0]
            if chr not in chr_list:
                continue
            left = l['loc'].loc['left']
            rite = l['loc'].loc['right']

            left_buck = ((left-1)//10000) * 10000
            right_buck = (rite//10000) * 10000
            buckets_reqd = range(left_buck, right_buck+10000, 10000)

            if chr not in gene:
                gene[chr] = {}

            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck not in gene[chr]:
                        gene[chr][buck] = []
                    gene[chr][buck].append([left, rite])
        #o.close()

        # Process the TEs:
        noverlap = []
        print(tefilename)
        if '.gz' in tefilename:
            o = gzip.open(tefilename,'rt')
        else:
            o = open(tefilename,'rU')

        for n, l in enumerate(o):
            t = l.strip().split('\t')
            chr = t[5].replace('chr', '')

            if chr not in chr_list:
                continue

            left = int(t[6])
            rite = int(t[7])

            left_buck = ((left-1)//10000) * 10000
            right_buck = (rite//10000) * 10000
            buckets_reqd = range(left_buck, right_buck+10000, 10000)

            if buckets_reqd:
                i = 1
                for buck in buckets_reqd:
                    if chr not in gene: # Should be very rare
                        continue
                    if buck not in gene[chr]:
                        continue

                    for k in gene[chr][buck]:
                        if left < k[1] and rite > k[0]:
                            i = 0
                            break
                    if i == 0: # already found an overlap, so quit out;
                        break
                if i == 1:
                    noverlap.append({'loc': location(chr=chr, left=left, right=rite), 'annot': t[10]})

            if n > 1000:
                break

        TEs = genelist()
        TEs.load_list(noverlap)

        #oh = gzip.open('_tmp/%s.exclusive.gz'%(tefilename),'wt')
        #for k in noverlap:
        #    oh.write(k)
        #oh.close()

        genes = genelist() # genelist('_tmp/%s.raw.bed.gz'%(genefilename), format=geneform, gzip=True)
        genes.load_list(raw)
        #TEs = genelist('_tmp/%s.exclusive.gz'%(tefilename), format=geneform, gzip=True)

        all_annot = genes + TEs
        if genome == 'mm':
            all_annot.save('mm10.exclusive.glb')
        elif genome == 'hs':
            all_annot.save('hg38.exclusive.glb')

    elif mode == 'inclusive':
        genes = genelist() # genelist('_tmp/%s.raw.bed.gz'%(genefilename), format=geneform, gzip=True)
        genes.load_list(raw)

        if tefilename.endswith('.gz'):
            TEs = genelist(tefilename, format=teform, gzip=True)
        else:
            TEs = genelist(tefilename, format=teform)
        all_annot = genes + TEs

        if genome == 'mm':
            all_annot.save('mm10.inclusive.glb')
        elif genome == 'hs':
            all_annot.save('hg38.inclusive.glb')

#     os.system('rm -r _tmp')

genomeIndex('mm','exclusive',
            'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz',
            'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz')

genomeIndex('mm','inclusive',
            'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz',
            'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz')

genomeIndex('hs','exclusive',
            'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz',
            'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz')

genomeIndex('hs','inclusive',
            'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz',
            'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz')



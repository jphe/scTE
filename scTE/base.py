import multiprocessing
import argparse
from functools import partial
import logging
import os, sys, glob, datetime, time, gzip
import collections
from collections import defaultdict
from math import log
from scTE.miniglbase import genelist, glload, location
from scTE.annotation import annoGtf

def read_opts(parser):
    args = parser.parse_args()
    if args.format == "BAM" :
        args.parser = "BAM"
    elif args.format == "SAM" :
        args.parser = "SAM"
    else :
        logging.error("The input file must be SAM/BAM format: %s !\n" % (args.format))
        sys.exit(1)

#     if not args.annoglb:
#         if not args.genefile:
#             logging.error("-gene option are needed without index refernece annotation file")
#             sys.exit(1)
#         else:
#             genefile = args.genefile[0]
#
#         if not args.tefile:
#             logging.error("-te option are needed without index refernece annotation file")
#             sys.exit(1)
#
#     if args.annoglb:
#         if args.genefile or args.tefile:
#             logging.error("-x option are exclusive to -te/-gene options")
#             sys.exit(1)
#
#     if args.mode not in ['inclusive', 'exclusive'] :
#         logging.error("Counting mode %s not supported\n" % (args.mode))
#         parser.print_help()
#         sys.exit(1)

    args.error = logging.critical
    args.warn = logging.warning
    args.debug = logging.debug
    args.info = logging.info

    args.argtxt ="\n".join(("Parameter list:", \
                "Sample = %s" % (args.out), \
                "Genome = %s" % (args.genome), \
#                 "TE file = %s" % (args.tefile), \
#                 "Gene file = %s" % (args.genefile), \
                "Reference annotation index = %s" %(args.annoglb), \
                "Minimum number of genes required = %s" % (args.genenumber), \
                "Minimum number of counts required = %s"% (args.countnumber),\
#                 "Mode = %s " % (args.mode), \
                "Number of threads = %s " % (args.thread),\
    ))
    return args

def getanno(filename, genefile, tefile, genome, mode):
    form ={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])', 'annot': 3}

    if genefile == 'default' and tefile == 'default':
        if genome == 'mm10':
            chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]
            if mode == 'exclusive':
                if not os.path.exists('mm10.exclusive.glb'):
                    logging.error("Did not find the annotation index mm10.exclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
                    sys.exit(1)
                all_annot = 'mm10.exclusive.glb'
                allelement = set(glload(all_annot)['annot'])

            elif mode == 'inclusive':
                if not os.path.exists('mm10.inclusive.glb'):
                    logging.error("Did not find the annotation index mm10.inclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
                    sys.exit(1)
                all_annot = 'mm10.inclusive.glb'
                allelement = set(glload(all_annot)['annot'])

        elif genome == 'hg38':
            chr_list = ['chr'+ str(i) for i in range(1,23) ] + [ 'chrX','chrY', 'chrM' ]
            if mode == 'exclusive':
                if not os.path.exists('hg38.exclusive.glb'):
                    logging.error("Did not find the annotation index hg38.exclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n" )
                    sys.exit(1)
                all_annot = 'hg38.exclusive.glb'
                allelement = set(glload(all_annot)['annot'])

            elif mode == 'inclusive':
                if not os.path.exists('hg38.inclusive.glb'):
                    logging.error("Did not find the annotation index hg38.inclusive.glb, you can download it from scTE github (www....) or either give the annotation with -te and -gene option \n")
                    sys.exit(1)
                all_annot = 'hg38.inclusive.glb'
                allelement = set(glload(all_annot)['annot'])
    else:
        if genome in ['hg38']:
            chr_list = ['chr'+ str(i) for i in range(1,23) ] + [ 'chrX','chrY', 'chrM' ]

        elif genome in ['mm10']:
            chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]

        if not os.path.isfile(tefile) :
            logging.error("No such file: %s !\n" %(tefile))
            sys.exit(1)

        if not os.path.isfile(genefile) :
            logging.error("No such file: %s !\n" % (genefile))
            sys.exit(1)

        all_annot = annoGtf(filename, genefile=genefile, tefile=tefile, mode=mode)
        allelement = set(glload(all_annot)['annot'])

    return(allelement,chr_list,all_annot)

def Readanno(filename, annoglb, genome):
    all_annot = glload(annoglb)
    allelement = set(all_annot['annot'])
    if genome in ['mm10']:
        chr_list = ['chr'+ str(i) for i in range(1,20) ] + [ 'chrX','chrY', 'chrM' ]
    elif genome in ['hg38']:
        chr_list = ['chr'+ str(i) for i in range(1,22) ] + [ 'chrX','chrY', 'chrM' ]
    return(allelement, chr_list, all_annot)

def Bam2bed(filename,out):
    if not os.path.exists('%s_scTEtmp/o1'%out):
        os.system('mkdir -p %s_scTEtmp/o1'%out)

    # Test to see which columns the barcode is in:
    o = open('%s.test.sh'%out,'w')
    st = 'samtools view -@ 2 %s | head | awk \'{OFS="\t"}{for(i=1;i<=NF;i++)if($i~/CR:Z:/)n=i}{for(i=1;i<=NF;i++)if($i~/UR:Z:/)m=i}{print n,m}\' > %s.test'%(filename,out)
    o.write(st)
    o.close()

    os.system('sh %s.test.sh'%out)

    o = open('%s.test'%out,'rU')
    for l in o:
        t = l.strip().split('\t')
        n=int(t[0])
        m=int(t[1])
    o.close()

    os.system('rm %s.test*'%out)

    if sys.platform == 'darwin': # Mac OSX has BSD sed
        os.system('samtools view -@ 2 %s | awk \'{OFS="\t"}{print $3,$4,$4+100,$%s,$%s}\' | sed -E \'s/CR:Z://g\' | sed -E \'s/UR:Z://g\'| awk \'!x[$0]++\' | gzip -c > %s_scTEtmp/o1/%s.bed.gz'%(filename,n,m,out,out))
    else:
        os.system('samtools view -@ 2 %s | awk \'{OFS="\t"}{print $3,$4,$4+100,$%s,$%s}\' | sed -r \'s/CR:Z://g\' | sed -r \'s/UR:Z://g\'| awk \'!x[$0]++\' | gzip > %s_scTEtmp/o1/%s.bed.gz'%(filename,n,m,out,out))
#     os.system('samtools view -@ 2 %s | awk \'{OFS="\t"}{for(i=1;i<=NF;i++)if($i~/CR:Z:/)n=i}{for(i=1;i<=NF;i++)if($i~/UR:Z:/)m=i}{print $3,$4,$4+100,$n,$m}\' | sed -r \'s/CR:Z://g\' | sed -r \'s/UR:Z://g\'| gzip > %s_scTEtmp/o1/%s.bed.gz'%(filename,out,out)) # need ~triple time


def splitAllChrs(chromosome_list, filename, genenumber, countnumber):
    '''
    **Purpose**
        Split the data into separate beds, and count up all the times each barcode appears

        This variant uses more memory, but does it all at the same time and gets the filtered whitelist for free

    **Returns**
        The barcode whitelist
    '''
    if not os.path.exists('%s_scTEtmp/o2' % filename):
        os.system('mkdir -p %s_scTEtmp/o2'%filename)

    file_handle_in = gzip.open('%s_scTEtmp/o1/%s.bed.gz' % (filename,filename), 'rt')
    file_handles_out= {chr: gzip.open('%s_scTEtmp/o2/%s.%s.bed.gz' % (filename,filename,chr), 'wt') for chr in chromosome_list}
    uniques = {chrom: set([]) for chrom in chromosome_list}

    CRs = defaultdict(int)

    # Make a unique BED
    for line in file_handle_in:
        t = line.strip().split('\t')
        chrom = t[0]
        CR = t[3]

        #if line in uniques[chrom]:
        #    continue
        #uniques[chrom].add(line)
        file_handles_out[chrom].write(line)

        CRs[CR] += 1
    [file_handles_out[k].close() for k in file_handles_out]
    file_handle_in.close()
    del uniques

    # Because this does it all in one go, you can just filter the whitelist here now, and don't need the .count. file;
    sortcb = sorted(CRs.items(), key=lambda item:item[1], reverse=True) # Sorts by the count;

    if not countnumber:
        mincounts = 2 * genenumber
    else:
        mincounts = countnumber

    whitelist = []
    for n, k in enumerate(sortcb):
        if k[1] < mincounts:
            break
        whitelist.append(k[0])

    return set(whitelist)
'''
def splitChr(chr,filename):
    if not os.path.exists('%s_scTEtmp/o2'%filename):
        os.system('mkdir -p %s_scTEtmp/o2'%filename)

    if chr == 'chr1':
        os.system('gunzip -c -f %s_scTEtmp/o1/%s.bed.gz | grep -v chr1\'[0-9]\' | grep %s | awk \'!x[$0]++\' | gzip -c > %s_scTEtmp/o2/%s.%s.bed.gz'%(filename,filename,chr,filename,filename,chr))
    elif chr == 'chr2':
        os.system('gunzip -c -f %s_scTEtmp/o1/%s.bed.gz | grep -v chr2\'[0-9]\' | grep %s | awk \'!x[$0]++\' | gzip -c > %s_scTEtmp/o2/%s.%s.bed.gz'%(filename,filename,chr,filename,filename,chr))
    else:
        os.system('gunzip -c -f %s_scTEtmp/o1/%s.bed.gz | grep %s | awk \'!x[$0]++\' | gzip -c > %s_scTEtmp/o2/%s.%s.bed.gz'%(filename,filename,chr,filename,filename,chr))

    CRs = {}
    o = gzip.open('%s_scTEtmp/o2/%s.%s.bed.gz'%(filename,filename,chr),'rb')
    for l in o:
        t = l.decode('ascii').strip().split('\t')
        if t[3] not in CRs:
            CRs[t[3]] = 0
        CRs[t[3]] += 1
    o.close()

    o = gzip.open('%s_scTEtmp/o2/%s.%s.count.gz'%(filename,filename,chr),'wt')
    for k in CRs:
        o.write('%s\t%s\n'%(k,CRs[k]))
    o.close()
'''
def align(chr, filename, all_annot, whitelist):
    '''
    **Purpose**
        For each read, align it to the index and assign a TE, gene.

    This is the speed critical part.

    '''
    s1 = time.time()

    buckets = all_annot.buckets[chr.replace('chr', '')]
    all_annot = all_annot.linearData

    oh = gzip.open('%s_scTEtmp/o2/%s.%s.bed.gz' % (filename, filename, chr), 'rt')
    res = {}
    for line in oh:
        t = line.strip().split('\t')

        barcode = t[3]
        if barcode not in whitelist:
            continue
        if barcode not in res:
            res[barcode] = defaultdict(int)

        chrom = t[0].replace('chr', '')
        left = int(t[1])
        rite = int(t[2])

        #loc = location(chr=chrom, left=left, right=rite)
        left_buck = ((left-1)//10000) * 10000
        right_buck = ((rite)//10000) * 10000
        buckets_reqd = range(left_buck, right_buck+10000, 10000)

        if buckets_reqd:
            loc_ids = set()
            loc_ids_update = loc_ids.update

            # get the ids reqd.
            [loc_ids_update(buckets[buck]) for buck in buckets_reqd if buck in buckets]
            #for buck in buckets_reqd:
            #    if buck in all_annot.buckets[chrom]:
            #        loc_ids.update(all_annot.buckets[chrom][buck])

            result = [all_annot[index]['annot'] for index in loc_ids if (rite >= all_annot[index]['loc'].loc['left'] and left <= all_annot[index]['loc'].loc["right"])]
            #for index in loc_ids:
            #    if rite >= all_annot[index]["loc"].loc['left'] and left <= all_annot[index]["loc"].loc["right"]:
            #        result.append(all_annot[index]['annot'])

            if result:
                for gene in result:
                    res[barcode][gene] += 1

    oh.close()

    if not os.path.exists('%s_scTEtmp/o3'%filename):
        os.system('mkdir -p %s_scTEtmp/o3'%filename)

    oh = gzip.open('%s_scTEtmp/o3/%s.%s.bed.gz' % (filename,filename,chr), 'wt')
    for bc in sorted(res):
        for gene in sorted(res[bc]):
            oh.write('%s\t%s\t%s\n' % (bc, gene, res[bc][gene]))
    oh.close()


def Countexpression(filename, allelement, genenumber, cellnumber):
    gene_seen = allelement

    whitelist={}
    o = gzip.open('%s_scTEtmp/o4/%s.bed.gz'%(filename, filename), 'rb')
    for n,l in enumerate(o):
        t = l.decode('ascii').strip().split('\t')
        if t[0] not in whitelist:
            whitelist[t[0]] = 0
        whitelist[t[0]] += 1
    o.close()

    CRlist = []
    sortcb = sorted(whitelist.items(), key=lambda item:item[1], reverse=True)
    for n,k in enumerate(sortcb):
        if k[1] < genenumber:
            break
        if n >= cellnumber:
            break
        CRlist.append(k[0])
    CRlist = set(CRlist)

    res = {}
    genes_oh = gzip.open('%s_scTEtmp/o4/%s.bed.gz' % (filename,filename), 'rb')
    for n, l in enumerate(genes_oh):
        t = l.decode('ascii').strip().split('\t')
        if t[0] not in CRlist:
            continue
        if t[0] not in res:
            res[t[0]] = {}
        if t[1] not in res[t[0]]:
            res[t[0]][t[1]] = 0
        res[t[0]][t[1]] += int(t[2])

    genes_oh.close()

    s=time.time()

    # Save out the final file
    res_oh = open('%s.csv'%filename, 'w')
    res_oh.write('barcodes,')
    res_oh.write('%s\n' % (','.join([str(i) for i in gene_seen])))

    gene_seen = list(gene_seen) # Do the sort once;
    gene_seen.sort()

    for k in sorted(res):
        l = ["0"] * len(gene_seen) # Avoid all the appends
        for idx, gene in enumerate(gene_seen):
            if gene in res[k]:
                l[idx] = str(res[k][gene])
        res_oh.write('%s,%s\n' % (k, ','.join(l)))
    res_oh.close()

    print('Detect %s cells expressed at least %s genes, results output to %s.csv'%(len(res),genenumber,filename))

'''
def filterCRs(filename,genenumber,countnumber):
    CRs = {}
    for f in glob.glob('%s_scTEtmp/o2/%s*.count.gz'%(filename,filename)):
        o = gzip.open(f,'rb')
        for l in o:
            t = l.decode('ascii').strip().split('\t')
            if t[0] not in CRs:
                CRs[t[0]] = 0
            CRs[t[0]] += int(t[1])
        o.close()

    sortcb=sorted(CRs.items(),key=lambda item:item[1],reverse=True)

    if not countnumber:
        mincounts = 2* genenumber
    else:
        mincounts = countnumber

    whitelist=[]
    for n,k in enumerate(sortcb):
        if k[1] < mincounts:
            break
        whitelist.append(k[0])

    return set(whitelist)
'''

def timediff(timestart, timestop):
        t  = (timestop-timestart)
        time_day = t.days
        s_time = t.seconds
        ms_time = t.microseconds / 1000000
        usedtime = int(s_time + ms_time)
        time_hour = int(usedtime / 60 / 60 )
        time_minute = int((usedtime - time_hour * 3600 ) / 60 )
        time_second =  int(usedtime - time_hour * 3600 - time_minute * 60 )
        retstr = "%dd %dh %dm %ds"  %(time_day, time_hour, time_minute, time_second,)
        return retstr

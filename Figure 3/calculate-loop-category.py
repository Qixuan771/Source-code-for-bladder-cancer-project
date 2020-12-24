import sys, bisect
import numpy as np
from collections import Counter

def check_in(p, List):

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        
        cache.add(tuple(q)) # return ATAC peak IDs
    
    return cache

def parse_bedpe1(fil, chip, gene_db, both=True):

    pool = []
    with open(fil, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            chrom = p[0]
            loci1 = [int(p[1])-10000, int(p[2])+10000]
            loci2 = [int(p[4])-10000, int(p[5])+10000]
            cache1 = check_in(loci1, chip[chrom])
            cache2 = check_in(loci2, chip[chrom])

            if both:
                if len(cache1) and len(cache2):
                    pool.append(p)
            else:
                if len(cache1) or len(cache2):
                    pool.append(p)
    
    return pool

def parse_bedpe(fil, chips1, chips2):

    classes = []
    total = 0
    with open(fil, 'r') as source:
        for line in source:
            total += 1
            p = line.rstrip().split()
            chrom = p[0]
            loci1 = [int(p[1])-25000, int(p[2])+25000]
            loci2 = [int(p[4])-25000, int(p[5])+25000]
            cache1_1 = check_in(loci1, chips1[chrom]) # promoter
            cache1_2 = check_in(loci2, chips1[chrom])
            cache2_1 = check_in(loci1, chips2[chrom]) # enhancer
            cache2_2 = check_in(loci2, chips2[chrom])
            label1 = []
            label2 = []
            if len(cache1_1):
                label1.append('promoter')
            if len(cache2_1):
                label1.append('enhancer')
            if len(cache1_2):
                label2.append('promoter')
            if len(cache2_2):
                label2.append('enhancer')
            if not label1:
                label1.append('None')
            if not label2:
                label2.append('None')
            keys = []
            for i in label1:
                for j in label2:
                    keys.append(tuple(sorted((i, j))))
            keys = set(keys)
            for k in keys:
                classes.append(k)
    
    counts = Counter(classes)
    
    return counts, total

def parse_bed(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            if not p[0] in D:
                D[p[0]] = []
            D[p[0]].append([int(p[1]), int(p[2])])

    for c in D:
        D[c].sort()

    return D

chips1 = parse_bed(sys.argv[2])
chips2 = parse_bed(sys.argv[3])

counts, total = parse_bedpe(sys.argv[1], chips1, chips2)
L = [round(counts[i]/total, 4) for i in [('enhancer', 'promoter'), ('enhancer', 'enhancer'), ('promoter', 'promoter'), ('None', 'enhancer'), ('None', 'promoter'), ('None', 'None')]]
print(L)

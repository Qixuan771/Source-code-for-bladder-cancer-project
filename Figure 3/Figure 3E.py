import bisect, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def get_loop_loci(loops):

    loci = {}
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        if not c1 in loci:
            loci[c1] = set()
        if not c2 in loci:
            loci[c2] = set()
        loci[c1].add((s1, e1))
        loci[c2].add((s2, e2))

    return loci

def parseLoops(loopfil):

    loops = []
    with open(loopfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            c1 = parse[0].lstrip('chr')
            c2 = parse[3].lstrip('chr')
            loops.append([c1, int(parse[1]), int(parse[2]), c2, int(parse[4]), int(parse[5])])
    
    return loops

def parseChIP(fil):
    
    data = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0].lstrip('chr')
            if '_' in chrom:
                continue
            if chrom in ['Y','M']:
                continue
            if chrom in data:
                data[chrom].append([int(parse[1]), int(parse[2])])
            else:
                data[chrom] = [[int(parse[1]), int(parse[2])]]
    for chrom in data:
        data[chrom].sort()
    
    return data

def centerLoci(loci, chip, halfr=50000, bins=50):

    arrays = []
    for c in loci:
        if not c in chip:
            continue
        bychrom = chip[c]
        for p in loci[c]:
            center = (p[0] + p[1]) // 2
            r = [c, center-halfr, center+halfr]
            data = np.zeros(r[2]-r[1])
            idx = max(0, bisect.bisect(bychrom, r[1:])-1)
            for q in bychrom[idx:]:
                if q[1] <= r[1]:
                    continue
                if q[0] >= r[2]:
                    break
                s = q[0]-r[1]
                if s < 0:
                    s = 0
                e = q[1]-r[1]
                data[s:e] += 1
            idx = np.linspace(0,r[2]-r[1],bins+1).astype(int)
            arr = []
            for i, j in zip(idx[:-1], idx[1:]):
                arr.append(data[i:j].mean())
            arr = np.r_[arr]
            arrays.append(arr)
    
    arr = np.r_[arrays].mean(axis=0)
    num = len(arrays)
    
    return arr, num

loopfil = sys.argv[1]
loops = parseLoops(loopfil)
loop_loci = get_loop_loci(loops)

chip1 = parseChIP('RT4_FOXA1_idr.txt')
chip2 = parseChIP('RT4_GATA3_1.txt')
halfr = 300000
arr1, num1 = centerLoci(loop_loci, chip1, halfr=halfr, bins=51)
arr2, num2 = centerLoci(loop_loci, chip2, halfr=halfr, bins=51)

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111)

l1, = ax.plot(arr1, color='#7FC97F', lw=3)
ax.set_xticks([0,25,50])
ax.set_xticklabels(['-300K','loop anchor','+300K'], fontsize=18)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel('Number of FOXA1 peaks per 10Kb', fontsize=20)

ax2 = ax.twinx()
l2, = ax2.plot(arr2, color='#BEAED4', lw=3)
#ax2.set_xticks([0,25,50])
#ax2.set_xticklabels(['-300K','loop anchor','+300K'], fontsize=18)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax2.set_ylabel('Number of GATA3 peaks per 10Kb', fontsize=20)

for spine in ['right', 'top']:
    ax.spines[spine].set_visible(False)

ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.yaxis.set_tick_params(width=2, labelsize=15)

for spine in ['left', 'top', 'bottom']:
    ax2.spines[spine].set_visible(False)

ax2.spines['right'].set_linewidth(2)
ax2.xaxis.set_tick_params(bottom=False, labelbottom=False)
ax2.yaxis.set_tick_params(width=2, labelsize=15)

ax.legend([l1, l2], ['FOXA1', 'GATA3'], frameon=False, fontsize=15)

#ax.xaxis.set_tick_params(width=2, labelsize=15)


#ax.set_title('YY1 binding at Peakachu loop loci', fontsize=20, pad=16)

plt.savefig('chip-enrich.pdf', dpi=300, bbox_inches='tight')
plt.close()

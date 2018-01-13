# -*- coding: utf-8 -*-
import sys, random
reload(sys)
sys.setdefaultencoding('utf-8')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def Load_Gene_Annotations(filename):
    Genes = []
    for line in open(filename, 'r').readlines():
        toks = line.split('\t')
        if len(toks) == 9 and toks[2] == "gene" and toks[0] != "Mito": # ignore mito genes
            Genes.append([toks[0], int(toks[3]), int(toks[4]), toks[8].split('"')[1], toks[6], [], toks[8].split('"')[3]])
        elif len(toks) == 9 and toks[2] == "CDS" and toks[0] != "Mito":
            gname = toks[8].split('"')[1]
            if gname == Genes[-1][3]:
                Genes[-1][5].append([int(toks[3]), int(toks[4])])
    return Genes

def Load_CDS_fa(filename):
    fasta = {}
    last = ''
    with open(filename, 'r') as fi:
        for i in fi:
            if i[0] == '>':
                last = i[1:-1].split(' ')[0]
                fasta[last] = ''
            else:
                fasta[last] += i[:-1]
    return fasta


def BEDGC(fib, genome):
    intervals = []
    total = 0.0; gc = 0.0
    with open(fib, 'r') as fin:
        for line in fin:
            cols = line[:-1].split('\t')
            if float(cols[4]) > 10.0:
                intervals.append([int(cols[1]), int(cols[2])])
    for i in intervals:
        for j in range(i[0], i[1]):
            total += 1
            if genome[j] == 'G' or genome[j] == 'C':
                gc += 1
    return gc / total

def PromotorGC(genes, genome):
    gc = 0.0; total = 0.0
    for g in genes:
        if g[4] == '+':
            for j in range(0, 200):
                total += 1
                if genome[g[0]][g[1]-j] == 'G' or genome[g[0]][g[1]-j] == 'C':
                    gc += 1
        else:
            for j in range(0, 200):
                total += 1
                if genome[g[0]][g[2]+j] == 'G' or genome[g[0]][g[2]+j] == 'C':
                    gc += 1
    return gc / total


genes = Load_Gene_Annotations('../../SacCer/Saccharomyces_cerevisiae.R64-1-1.86.gtf')
genes_cds_fa = Load_CDS_fa('../../SacCer/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa')
genome_fa = Load_CDS_fa('../../SacCer/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa')

gc_chr_nhp6 = []
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr1.bed', genome_fa['I']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr2.bed', genome_fa['II']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr3.bed', genome_fa['III']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr4.bed', genome_fa['IV']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr5.bed', genome_fa['V']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr6.bed', genome_fa['VI']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr7.bed', genome_fa['VII']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr8.bed', genome_fa['VIII']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr9.bed', genome_fa['IX']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr10.bed', genome_fa['X']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr11.bed', genome_fa['XI']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr12.bed', genome_fa['XII']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr13.bed', genome_fa['XIII']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr14.bed', genome_fa['XIV']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr15.bed', genome_fa['XV']))
gc_chr_nhp6.append(BEDGC('/Users/jalgard/Yandex.Disk.localized/Jalgard/Legacy Projects/Nastya/Supp_Nhp6A_BEDfiles/DowellSupp_Nhp6Chr16.bed', genome_fa['XVI']))

gc_nhp6 = np.median(gc_chr_nhp6)

global_gc = 0.0; global_len = 0.0
for i in genome_fa:
    for nuc in genome_fa[i]:
        global_len += 1
        if nuc == 'G' or nuc == 'C':
            global_gc += 1
global_gc = global_gc/global_len

promotGC = PromotorGC(genes, genome_fa)

gene_dict = defaultdict(lambda : [])
gene_per_chr = defaultdict(lambda: 0.0)
for gene in genes:
    gene_dict[gene[3]] = gene
    gene_per_chr[gene[0]] += 1


genes_nhp6_cds = [x[:-2] for x in open('../../SacCer/nhp6A_orfBinding_geneList.txt','r').readlines()[1:]]
genes_nhp6_prm = [x.split('\t')[4] for x in open('../../SacCer/nhp6A_bindingGrps_fig2kmeansClst.txt','r').readlines()[1:]]
genes_non_nhp6 = [x[:-2] for x in open('../../SacCer/nhp6A_unbound_geneList.txt','r').readlines()[1:]]

gene_len_bins = defaultdict(lambda: [])
for i in genes_cds_fa:
    l = len(genes_cds_fa[i])
    gene_len_bins[l/50].append(l)

gene_gc_bins = defaultdict(lambda: [])
for i in genes_cds_fa:
    l = len(genes_cds_fa[i])
    gc = 0.0
    for nuc in genes_cds_fa[i]:
        if nuc == 'G' or nuc == 'C':
            gc+=1
    gene_gc_bins[l/50].append(gc/len(genes_cds_fa[i]))

gene_name_bins = defaultdict(lambda: [])
for i in genes_cds_fa:
    l = len(genes_cds_fa[i])
    gene_name_bins[l/50].append(i)

gene_len_dist_nhp6 = defaultdict(lambda: 0)
for i in genes_cds_fa:
    l = len(genes_cds_fa[i])
    if i in genes_nhp6_cds:
        gene_len_dist_nhp6[l/50] += 1

def pick_short(nhp6, all):
    result = []
    for i in nhp6:
        result.extend(random.sample(all[i], nhp6[i]))
    return result


gc_content_cds = []
gc_content_all = []
gc_content_srt = []
for i in genes_cds_fa:
    if i in genes_nhp6_cds:
        gc = 0.0; at = 0.0; total = 0
        for nuc in genes_cds_fa[i]:
            total += 1
            if nuc == 'G' or nuc == 'C':
                gc+=1
            else:
                at+=1
        gc_content_cds.append(gc/total)

    else:
        gc = 0.0; at = 0.0; total = 0
        for nuc in genes_cds_fa[i]:
            total += 1
            if nuc == 'G' or nuc == 'C':
                gc+=1
            else:
                at+=1
        gc_content_all.append(gc/total)

# bootstrap short genes
short_gene_bootstrap = []
for i in range(100):
    short_gene_bootstrap.append(pick_short(gene_len_dist_nhp6, gene_gc_bins))

short_q10 = [np.percentile(np.array(x,dtype=np.float64), 10) for x in short_gene_bootstrap]
short_q90 = [np.percentile(np.array(x,dtype=np.float64), 90) for x in short_gene_bootstrap]
short_q25 = [np.percentile(np.array(x,dtype=np.float64), 25) for x in short_gene_bootstrap]
short_q50 = [np.percentile(np.array(x,dtype=np.float64), 50) for x in short_gene_bootstrap]
short_q75 = [np.percentile(np.array(x,dtype=np.float64), 75) for x in short_gene_bootstrap]

fig, ax = plt.subplots(2, 2)

plt.sca(ax[0, 1])
axes = plt.gca()
D = plt.boxplot([gc_content_cds, gc_content_all, short_gene_bootstrap[0]], patch_artist=True, showfliers=False)

D['medians'][2].set_ydata([np.median(short_q50), np.median(short_q50)])
path = D['boxes'][2].get_path()
path.vertices[0][1] = np.average(short_q25)
path.vertices[1][1] = np.average(short_q25)
path.vertices[2][1] = np.average(short_q75)
path.vertices[3][1] = np.average(short_q75)
path.vertices[4][1] = np.average(short_q25)
D['whiskers'][4].set_ydata(np.array([np.average(short_q10), np.average(short_q25)]))
D['whiskers'][5].set_ydata(np.array([np.average(short_q90), np.average(short_q75)]))
D['caps'][4].set_ydata(np.array([np.average(short_q10), np.average(short_q10)]))
D['caps'][5].set_ydata(np.array([np.average(short_q90), np.average(short_q90)]))


D['boxes'][0].set(facecolor="#FFFFFF", linewidth=2)
D['boxes'][1].set(facecolor="#444444", linewidth=1)
D['boxes'][2].set(facecolor="#888888", linewidth=1)

D['medians'][0].set(color="#000000")
D['medians'][1].set(color="#000000")
D['medians'][2].set(color="#000000")

plt.axhline(global_gc, ls='--', c='#000000')
plt.axhline(gc_nhp6, ls=':', c='#000000')
plt.axhline(promotGC, ls='-', c='#000000')


axes.set_title('GC состав')
axes.set_xlabel('Выборка генов')
axes.set_ylabel('% GC')




plt.sca(ax[0, 0])
axes = plt.gca()

gene_len_all = []
gene_len_nhp6 = []
gene_len_small = pick_short(gene_len_dist_nhp6, gene_len_bins)

for i in genes_cds_fa:
    if i in genes_nhp6_cds:
        gene_len_nhp6.append(len(genes_cds_fa[i]))

    else:
        gene_len_all.append(len(genes_cds_fa[i]))


DL = plt.boxplot([gene_len_nhp6, gene_len_all, gene_len_small], patch_artist=True, showfliers=False)

DL['boxes'][0].set(facecolor="#FFFFFF", linewidth=2)
DL['boxes'][1].set(facecolor="#444444", linewidth=1)
DL['boxes'][2].set(facecolor="#888888", linewidth=1)

DL['medians'][0].set(color="#000000")
DL['medians'][1].set(color="#000000")
DL['medians'][2].set(color="#000000")

axes.set_ylim([0,2300])
axes.set_title('Длина гена')
axes.set_xlabel('Выборка генов')
axes.set_ylabel('Длина, п.н.')



plt.sca(ax[1, 0])
axes = plt.gca()

short0 = pick_short(gene_len_dist_nhp6, gene_name_bins)

CUnph6 = defaultdict(lambda: 0)
CUshort = defaultdict(lambda: 0)
CUall = defaultdict(lambda: 0)

for i in genes_cds_fa:
    if i in genes_nhp6_cds:
        for j in range(0, len(genes_cds_fa[i]), 3):
            CUnph6[genes_cds_fa[i][j:j+3]] += 1
    elif i in short0:
        for j in range(0, len(genes_cds_fa[i]), 3):
            CUshort[genes_cds_fa[i][j:j+3]] += 1
    else:
        for j in range(0, len(genes_cds_fa[i]), 3):
            CUall[genes_cds_fa[i][j:j+3]] += 1

CUdata_nhp6 = [CUnph6[x] for x in sorted(CUnph6)]
CUdata_all = [CUall[x] for x in sorted(CUnph6)]
CUdata_short = [CUshort[x] for x in sorted(CUnph6)]

Trip1 = ['#FF0000' if 'CG' in x or 'GC' in x else '#000000' for x in sorted(CUnph6) ]
Trip2 = ['#FF0000' if 'G' in x else '#000000' for x in sorted(CUnph6) ]

CUdata_nhp6 = [float(x)/sum(CUdata_nhp6) for x in CUdata_nhp6]
CUdata_all = [float(x)/sum(CUdata_all) for x in CUdata_all]
CUdata_short = [float(x)/sum(CUdata_short) for x in CUdata_short]

plt.scatter(CUdata_all, CUdata_nhp6, c='#FFFFFF', edgecolor='#000000', marker='^', s=30)
plt.scatter(CUdata_all, CUdata_short, c='#000000', edgecolor='#000000', marker='s', s=10)

axes.set_ylim([0.0, 0.06])
axes.set_xlim([0.0, 0.06])
axes.set_title('Частота использования триплетов')
axes.set_xlabel('частота в (2)')
axes.set_ylabel('частота в (1) / (3)')

plt.show()

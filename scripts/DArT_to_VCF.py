#!/home/dariza/bin/python3

# -*- coding: utf-8 -*-

"""
Created on Mon Aug 10 12:25:00 2020

@author: darizasu
"""

import argparse
import sys


# %% Input arguments


prs = argparse.ArgumentParser(description="This script converts DArT 2-row format into a VCF file")

prs.add_argument("-dartG",
                 type = str,
                 help = "DArT 2-row file with genotypes",
                 required = True,
                 metavar = 'DArT')

prs.add_argument("-dartC",
                 type = str,
                 help = "DArT 2-row file with allele counts",
                 required = False,
                 metavar = 'DArT')

prs.add_argument("-vcf",
                 type = str,
                 help = "Name of the output VCF file",
                 required = True,
                 metavar = 'VCF')

prs.add_argument("-ref", type = str,
                 required = True,
                 help = 'A comma-separated file with SNP marker information. 1st col: CHROM. 2nd col: POS. 3rd col: ID. 4th col: REF. 5th col: ALT',
                 metavar = 'CSV')

prs.add_argument("-mn",
                 type = str,
                 help = "Name of the column with marker name",
                 required = False,
                 metavar = 'MarkerName',
                 default = 'MarkerName')

prs.add_argument("-at",
                 type = str,
                 help = "Name of the column with allele type",
                 required = False,
                 metavar = 'AlleleType',
                 default = 'AlleleType')

prs.add_argument("-aq",
                 type = str,
                 help = "Name of the column with allele sequence",
                 required = False,
                 metavar = 'AlleleSeq',
                 default = 'AlleleSequence')

# ar = prs.parse_args('-dartG Report_DCob20-5378_tgSNP_1.csv -dartC Report_DCob20-5378_tgRawCounts_2.csv -vcf Report_DCob20-5378_RELIABLE_SNPs.vcf -ref ref_markers.csv'.split())
# ar = prs.parse_args('-dartG Report_DCob20-5378_tgSNP_1.csv -vcf output.vcf -ref ref_markers.csv'.split())
ar = prs.parse_args()

# prs.print_help()
# print(ar)


# %% Functions


def ref (reference):

    with open(reference, 'r') as reference:

        r = dict()
        for m in reference:

            if m.startswith('#'):
                continue

            else:
                m = m.strip().split(',')
                r[m.pop(2)] = m

        print('Read ', len(r), ' markers from the reference file "',
              reference, '".', file = sys.stderr)

        return r , [ '##contig=<ID=' + i + '>\n' 
                   for i in sorted( {v[0] for k,v in r.items()} ) ]
        # {mrk: ['CHROM','POS','STRAND','REF','ALT']}


def findSNP (ref, snp):

    d = sum( x != y for x, y in zip(ref, snp) )

    if d > 1:
        ref = ref[::-1]
        ref = "".join([compSNP(i) for i in ref])
        d = sum( x != y for x, y in zip(ref, snp) )
        if d > 1:
            print('ERROR EXIT', ref, snp)

    elif d == 0:
        return None , None

    for i in range(len(ref)):
        if ref[i] != snp[i]:
            return ref[i] , snp[i]


def dartR (dart, mn, at, aq, first = False):

    h = False
    m = dict() if not bool(first) else first

    with open(dart, 'r') as file:

        for l in file:
            l = l.strip().split(',')

            if l[0] == '*':
                info = l.count('*')

            elif h:
                mt = [i if i != '-' else '9' for i in l[info:]]

                if not bool(first):

                    if l[mn] in m:
                        m[l[mn]] [l[at]] = [ l[aq] , mt ]
                        # print(l[mn])
                        r , s = findSNP(m[l[mn]] ['Ref'][0] , m[l[mn]] ['Snp'][0])

                        if not (r or s):
                            print('The marker', l[mn], 'does not show any variation in the reported sequence:', file = sys.stderr)
                            print(m[l[mn]] ['Ref'][0], m[l[mn]] ['Snp'][0], file = sys.stderr, sep = '\n')
                            del m[l[mn]]

                        else:
                            m[l[mn]] ['Ref'][0] , m[l[mn]] ['Snp'][0] = r , s

                    else:
                        m[l[mn]] = {l[at]: [ l[aq] , mt ]}

                else:
                    m[l[mn]] [l[at]].append(mt)

            else:
                h = l
                mn , at , aq = h.index(mn) , h.index(at) , h.index(aq)

        
        print('Read ', len(m), ' markers from the DArT file "',
              file, '".', file = sys.stderr)
        return m , h[info:]
        # {'mrk': {'ref': [[REF],[GT],[DP]],
        #          'snp': [[ALT],[GT],[DP]]}}


def compSNP(r, s = False):

    a = {'A':'T', 'T':'A', 'C':'G', 'G':'C',
         'R':'Y', 'Y':'R', 'K':'M', 'M':'K', 'S':'S', 'W':'W',
         'N':'N', 'B':'B', 'D':'D', 'H':'H', 'V':'V', 'X':'X', '-':'-'}
    if s:
        return a[r] , a[s]
    else:
        return a[r]


def errmsg(kerr, aerr, merr, awarn):

    if kerr:
        print('\nThe following markers were present in the DArT file but not in the reference markers file.\n',
              'They were not included in the final VCF:\n\n',
              '\n'.join(['\t' + i for i in kerr]), file = sys.stderr, end = '\n', sep = '')

    if aerr:
        print('\nThe following markers had inconsistent alleles between the DArT and the reference markers file.\n',
              'They were not included in the final VCF:\n\n',
              '\n'.join(['\t' + i for i in aerr]), file = sys.stderr, end = '\n', sep = '')

    if merr:
        print('\nThe following markers do not have genotype calls for any sample.\n',
              'They were not included in the final VCF:\n\n',
              '\n'.join(['\t' + i for i in merr]), file = sys.stderr, end = '\n', sep = '')

    if awarn:
        print('\nThe following markers have switched REF and ALT alleles between the DArT and the reference markers file.\n',
              'They were included in the VCF using the alleles in the reference markers file.\n\n',
              '\n'.join(['\t' + i for i in awarn]), file = sys.stderr, end = '\n', sep = '')

    return None


def infoW(f):

    NS = len(f) - sum('./.' in i for i in f)
    if NS == 0:
        return None 

    AC_Het = sum('0/1' in i for i in f)
    AC_Hom = sum('1/1' in i for i in f)
    MAF = ((min( sum('0/0' in i for i in f), AC_Hom ) * 2) + AC_Het) / (NS * 2)
    AC_Hom *= 2

    AN = 'AN=' + str(NS * 2)
    NS = 'NS=' + str(NS)
    MAF = 'MAF=' + str( round(MAF, 6) )
    AC = 'AC=' + str(AC_Hom + AC_Het)
    AC_Het = 'AC_Het=' + str(AC_Het)
    AC_Hom = 'AC_Hom=' + str(AC_Hom)

    return [';'.join([AN, NS, MAF, AC, AC_Het, AC_Hom])]


def lineW(mar, dartC):

    g = {'10': '0/0', '11': '0/1', '01': '1/1', '99': './.'}

    if dartC:

        f = list()
        ref , alt = mar['REF'][0] , mar['ALT'][0]

        for a , b , c , d in zip( mar['REF'][1], mar['ALT'][1],
                                  mar['REF'][2], mar['ALT'][2] ):
            bsdp = {'A':'0','C':'0','G':'0','T':'0'}
            bsdp[ref] , bsdp[alt] = c , d
            bsdp = ','.join([v for k,v in bsdp.items()])
            f.append(g[a + b] + ':' + str(int(c) + int(d)) + ':' + bsdp)

        return f

    else:
        return [g[a + b] for a, b in zip(mar['REF'][1], mar['ALT'][1])]



def vcfW(vcf, dartG, reference, mn, at, aq, dartC = False):

    m , h = dartR(dartG, mn, at, aq)

    header = ['##fileformat=VCFv4.2\n',
              '##FILTER=<ID=PASS,Description="All filters passed">\n',
              '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples genotyped">\n',
              '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">\n',
              '##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of alleles in called genotypes">\n',
              '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">\n',
              '##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">\n',
              '##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">\n',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(h) + '\n']

    r , c = ref(reference)
    header[9:9] = c

    line = list()
    FORMAT = 'GT'

    if dartC:
        m , h = dartR(dartC, mn, at, aq, first = m)
        header[9:9] = [ '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n',
                        '##FORMAT=<ID=BSDP,Number=4,Type=Integer,Description="Number of base calls (depth) for the 4 nucleotides in called SNVs sorted as A,C,G,T">\n' ]
        FORMAT += ':DP:BSDP'

    with open(vcf, 'w') as vcf:

        vcf.writelines(header)
        
        kerr = list()
        aerr = list()
        merr = list()
        awarn = list()

        for k in list(m): # Avoid Runtime error because of dict changing size

            if r.get(k):
                if r[k][2] == '-':
                    m[k]['Ref'][0] , m[k]['Snp'][0] = compSNP(m[k]['Ref'][0], m[k]['Snp'][0])

                if set(r[k][3:]) != {m[k]['Ref'][0], m[k]['Snp'][0]}:
                    aerr.append(k + '\t' + m[k]['Ref'][0] + ':' + m[k]['Snp'][0] + ' != ' + r[k][3] + ':' + r[k][4])
                    continue

                elif r[k][3] == m[k] ['Ref'][0]:
                    m[k] ['REF'] = m[k].pop('Ref')
                    m[k] ['ALT'] = m[k].pop('Snp')

                else:
                    awarn.append(k + '\t' + m[k]['Ref'][0] + ':' + m[k]['Snp'][0] + ' ~= ' + r[k][3] + ':' + r[k][4])
                    m[k] ['REF'] = m[k].pop('Snp')
                    m[k] ['ALT'] = m[k].pop('Ref')

                SAMPLE = lineW(m[k], dartC)
                INFO = infoW(SAMPLE)
                if not INFO:
                    merr.append(k)
                    continue

                CPIRA = r[k]
                CPIRA[2] = k

                line.append(CPIRA + ['255'] + ['PASS'] + INFO + [FORMAT] + SAMPLE)

            else:
                kerr.append(k)

        errmsg(kerr, aerr, merr, awarn)

        for marker in sorted(line, key = lambda x: (x[0], int(x[1]))):
            vcf.write('\t'.join(marker) + '\n')

        return None


# %% Run

if __name__ == '__main__':

    vcfW(vcf = ar.vcf,
         dartG = ar.dartG,
         dartC = ar.dartC,
         reference = ar.ref,
         mn = ar.mn,
         at = ar.at,
         aq = ar.aq)


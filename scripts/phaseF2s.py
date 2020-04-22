#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 20:21:48 2020

@author: darizasu
"""


import gzip
import argparse
import re


# %% Input arguments

prs = argparse.ArgumentParser(description="This script does...")

prs.add_argument("-vcf",
                 type = str,
                 help = "VCF file",
                 required = True,
                 metavar = 'VCF')

prs.add_argument("-out",
                 type = str,
                 help = "Name of the output VCF file",
                 required = True,
                 metavar = 'VCF')

prs.add_argument("-ped", type = str,
                 required = True,
                 help = 'A comma-separated file with a single family per line. 1st col: Pedigree ID. 2nd col: Parent1. 3rd col: Parent3. 4th col: Parent4. 5th - N cols: Every F2 individual, one per column. Missing data must be coded as NA',
                 metavar = 'CSV')

prs.add_argument("-err", type = str,
                 required = False,
                 help = 'Name of the output CSV file with the number and percentage of family errors for each F2.',
                 metavar = 'CSV')

# prs.print_help()

# ar = prs.parse_args('-vcf toy_test.vcf.gz -out toy_phased_F2s.vcf -ped pedigrees_toy.csv'.split())
# ar = prs.parse_args('-vcf Report_DCob20-5075_GBSext.vcf.gz -out phased_F2s.vcf -ped pedigrees.csv'.split())
ar = prs.parse_args('-vcf Report_DCob20-5075_GBSext_imputed.vcf.gz -out phased_F2s.vcf -ped pedigrees.csv -err family_errors.csv'.split())

# print(ar)


# %% Functions


def hyLen (spGT):
    """

    Parameters
    ----------
    spGT : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    up = sorted( list( {i for i in re.split('\W+', '/'.join(spGT))} ) )

    if len(up) == 3:
        del up[1]

    elif len(up) == 1:
        up = up + up

    return '/'.join(up)



def hyPar (pGT):
    """

    Parameters
    ----------
    pGT : TYPE
        DESCRIPTION.

    Returns
    -------
    pGT : TYPE
        DESCRIPTION.

    """

    if len(pGT) == 4:

        up1 = hyLen(pGT[:2])
        up2 = hyLen(pGT[2:])
        pGT[0] , pGT[1] = up1 , up2
        del pGT[2:4]

    elif len(pGT) == 3:

        up2 = hyLen(pGT[1:])
        pGT[1] = up2
        del pGT[2]

    return pGT


def readPeds (peds):
    """

    Parameters
    ----------
    peds : TYPE
        DESCRIPTION.

    Returns
    -------
    fs : TYPE
        DESCRIPTION.

    """

    with open(peds, 'r') as peds:

        i = 0
        fs = dict()
        cross = dict()

        for ped in peds:

            if ped.startswith('#'):
                continue

            else:

                ped = ped.strip().split(',')
                fs[str(i)] = {'p' : [i for i in ped[1:5] if i != 'NA' ]}
                fs[str(i)]['off'] = [i for i in ped[5:] if i != 'NA']

                for v in fs[str(i)]['off']:
                    cross[v] = ped[0]

                i += 1

        return fs , cross


def phase (pGT, oGT, h, off, eoff):
    """

    Parameters
    ----------
    pGT : TYPE
        DESCRIPTION.
    oGT : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    off : TYPE
        DESCRIPTION.
    eoff : TYPE
        DESCRIPTION.

    Returns
    -------
    ph : TYPE
        DESCRIPTION.
    eoff : TYPE
        DESCRIPTION.

    """

    ph = list()

    for n,i in enumerate(oGT):

        if (len({i for i in pGT}) < 2):

            upGT = {i for i in re.split('\W+', '/'.join(pGT))}

            con = [all([(j in upGT) for j in i.split('/')]),
                   upGT == {'.'},
                   i == './.' and len(upGT) == 2 ]

            if any(con):
                ph.append(i.replace('/','|'))

            elif i == './.' and len(upGT) == 1:
                ph.append(pGT[0].replace('/','|'))

            else:
                eoff[h[off[n]]] += 1
                ph.append(pGT[0].replace('/','|'))

        else:

            if i == '0/1':

                if pGT == ['0/1','1/1']:
                    ph.append(i.replace('/','|'))

                elif '1' in {i for i in re.split('\W+', pGT[0])}:
                    ph.append('1|0')

                else:
                    ph.append(i.replace('/','|'))

            else:
                ph.append(i.replace('/','|'))

    return ph , eoff


def procVCF(VCF, peds, oVCF):
    """

    Parameters
    ----------
    VCF : str
        DESCRIPTION.
    peds : str
        DESCRIPTION.
    oVCF : str
        DESCRIPTION.

    Returns
    -------
    eoff : dict
        DESCRIPTION.

    """

    with gzip.open(VCF, 'r') as VCF, open(oVCF, "a") as oVCF:

        fs , cross = readPeds(peds)
        l = 0

        for line in VCF:

            line = line.decode().strip().split('\t')

            if len(line) < 2:
                None

            elif line[0].startswith('#CHROM'):

                h = line
                eoff = {i: 0 for i in h[9:]}
                u = {h.index(j) for i in fs.values() for j in i['off'] }
                u = [i for i in range(9)] + sorted(list(u))
                line = [line[i] for i in u]

            else:

                l += 1
                GT = line[8].split(':').index('GT')

                for f in fs.keys():

                    p = [h.index(i) for i in fs[f]['p']]
                    off = [h.index(i) for i in fs[f]['off']]

                    oGT = [line[i].split(':')[GT] for i in off]
                    pGT = [line[i].split(':')[GT] for i in p]

                    pGT = hyPar(pGT)
                    ph , eoff = phase(pGT, oGT, h, off, eoff)

                    for n,i in enumerate(off):
                        line[i] = ph[n]

                line = [line[i] for i in u]

            oVCF.write('\t'.join(line) + '\n')

        h = [h[i] for i in range(len(h)) if i in u]

        eoff = [
            [ k, cross[k], str( round(v * 100 / l, 3) ), str(v) ]
             for k,v in eoff.items() if k in h
             ]

        return eoff


# %% Run

if __name__ == '__main__':

    eoff = procVCF(VCF = ar.vcf, peds = ar.ped, oVCF = ar.out)

    if ar.err:
        with open(ar.err, 'w') as fe:
            fe.write('Line,Pedigree,Percentage_of_errors,Number_of_errors\n')
            fe.writelines([','.join(v) + '\n' for v in eoff])


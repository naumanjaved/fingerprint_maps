import subprocess
import os
import itertools as it
import numpy as np
import sys
import argparse
import traceback
import time

def extract_similar_SNPs(chrom, VCF_file, int_directory, SIM):
    '''
    Writes to a new file the SNPs in an input VCF with population specific
    minor allele fractions(MAFs) that do not differ by more than SIM. For
    example,this function would reject a SNP with EUR_AF - AMR_AF  = 0.15 if
    SIM = 0.1, but would accept if SIM = 0.15.

    Parameters
    ----------
    chrom : string
        chromosome number
    VCF_directory : string
        directory path where original VCFs from 1000genomes phase 3 are stored
    int_directory : string
        directory path where intermediate and output files are stored
    SIM : float
        similarity value - the maximum pairwise difference in population
        specific MAF above which a SNP is rejected

    Returns
    -------
    None
    '''
    # loop through each line of the vcf, ignoring the 253 line header
    vcf = open(VCF_file, 'r')
    SNPs = []
    with vcf as variants:
        for k, line in enumerate(variants):
            if k > 252:
                parsed = line.rstrip('\n').split('\t')
                # exclude indels, multi-allelic SNPs
                if len(str(parsed[4])) != 1:
                    continue
                if len(str(parsed[3])) != 1:
                    continue
                SNP = parsed[2]
                # create string for each MAF values
                AFs =parsed[7].split(';')
                overallval = 0.0
                # initialize lists holding AFs
                values = []
                # initialize lists holding differences between population
                # AFs
                differences = []
                # loop through AF values and append to values
                for info in AFs:
                    if len(info.split('=')) != 2:
                        continue
                    entry = info.split('=')[0]
                    value = info.split('=')[1]
                    if entry == "AF":
                        overallval = float(value)
                    if entry == 'AFR_AF':
                        values.append(float(value))
                    if entry == 'AMR_AF':
                        values.append(float(value))
                    if entry == 'EUR_AF':
                        values.append(float(value))
                    if entry == 'EAS_AF':
                        values.append(float(value))
                    if entry == 'SAS_AF':
                        values.append(float(value))
                # compute pairwise differences for all pairs in values list
                differences = [abs(y-x) for x,y in it.combinations(values,2)]
                diffbool = True
                # if any pairwise difference exceeds 0.10 or one population
                # is not represented, then skip current SNP
                for diff in differences:
                    if diff > SIM:
                        diffbool = False
                        break
                    if not diffbool:
                        continue
                if len(values) != 5:
                        continue
                if SNP in SNPs or SNP == ".":
                    continue
                # append to list of "similar" SNPs
                with open(int_directory + "chr_" + chrom \
                    + "-common-SNPs.list", 'a') as similar_SNPs:
                    similar_SNPs.write(SNP + '\n')

def create_VCFs(chrom, VCF_file, int_directory, min_MAF):
    '''
    Script that uses vcftools to extract SNPs from each chromosome
    in the 1000 genomes phase 3 folder. These SNPs are:
    1) biallelic
    2) have population specific MAF differences no greater than 10% as
       determined by being listed in chr_i-common-SNPs.list
    3) are phased

    Parameters
    ----------
    chrom : string
        chromosome number
    VCF_file : string
        path to 1000 genomes VCF file
    int_directory : string
        directory path where intermediate and output files are stored
    min_MAF : float
        filter out all SNPs with population averaged MAF less than this value
    Returns
    -------
    None
    '''
    subprocess.check_call("vcftools --vcf " + VCF_file + " " \
    + "--maf " + str(min_MAF) + " " \
    + "--phased " \
    + "--remove-indels " \
    + "--min-alleles 2 --max-alleles 2 " \
    + "--recode --recode-INFO-all " \
    + "--snps " + int_directory + "chr_" + chrom + "-common-SNPs.list " \
    + "--out " + int_directory + "chr_" + chrom, shell=True)


def sort_VCF(chrom, int_directory):
    '''
    Sometimes VCFs have duplicate IDs or misformed ID names. This script
    simply sorts and rewrites the new vcf.

    Parameters
    ----------
    chrom : string
        chromosome number
    VCF_file : string
        path to 1000 genomes VCF file
    int_directory : string
        directory path where intermediate and output files are stored
    min_MAF : float
        filter out all SNPs with population averaged MAF less than this value

    Returns
    -------
    None
    '''
    # sort
    subprocess.check_call("tail -n+253 " + int_directory + "chr_" + chrom + ".recode.vcf" \
                            + " | awk -F'\t' '{print $3}' | sort | uniq -d > " \
                            + int_directory + chrom + ".duplicate", shell=True)
    # new file to write duplicate IDs to
    duplicates = open(int_directory + chrom + ".duplicate", 'r')
    snplist = []
    # determine and write duplicate IDs
    with duplicates as d:
        for snp in d:
            snplist.append(snp.rstrip('\n'))
    vcf = open(int_directory + "chr_" + chrom + ".recode.vcf", 'r')
    newvcf = open(int_directory + chrom + ".recode_u.vcf", 'w')
    # rewrite unique IDs to new file
    with vcf as f:
        for k, line in enumerate(f):
            if k <= 252:
                newvcf.write(line)
            if k > 252:
                if line.split('\t')[2] not in snplist:
                    newvcf.write(line)
    newvcf.close()

def create_PLINK_binary(chrom, int_directory, recomb_directory):
    '''
    Crates PLINK binary files(https://www.cog-genomics.org/plink2/formats) from
    input VCF for use with LDSC script(https://github.com/bulik/ldsc)

    Parameters
    ----------
    chrom : string
        chromosome number
    VCFs_directory : string
        directory path where original VCFs from 1000genomes phase 3 are stored
    int_directory : string
        directory path where intermediate and output files are stored
    recomb_directory : string
        directory path where 1000genomes recombination map files are stored

    Returns
    -------
    None
    '''
    recomb_file = recomb_directory \
                    + "genetic_map_chr" \
            + chrom \
            + "_combined_b37.txt* " \
            + chrom
    # if autosome, use 1000 G recombination map to write centimorgan
    # positions for the .bim file
    if chrom != "X":
        subprocess.check_call("plink --vcf " \
        + int_directory + chrom + ".recode_u.vcf" \
        + " --cm-map " + recomb_file \
    + " --make-bed" \
        + " --out " + int_directory + chrom, shell=True)
    # if x chromsome then just use basepair position instead of centimorgans
    else:
        subprocess.check_call("plink --vcf " \
    + int_directory + chrom + ".recode_u.vcf" \
    + " --make-bed " \
        + " --out " + int_directory + chrom, shell=True)

def LD_score(chrom, int_directory, LD_script, window_cm, window_kb):
    '''
    Calculates LDscore of all variants for each chromosome
    using LDSC script(https://github.com/bulik/ldsc)

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    LD_path : string
        directory containing python LD script
    window_cm : float
        for autosomes, the window in centimorgans over which to calculate LD
        scores
    window_kb : int
        for sex chromosome, the window in kilobases over which to calculate LD
        score

    Returns
    -------
    None
    '''
    if chrom != "X":
        subprocess.check_call("python " + LD_script \
        + " --bfile " + int_directory + chrom \
        + " --ld-wind-cm " + str(window_cm) \
        + " --out " + int_directory + "LD-" + chrom \
        + " --l2 --yes-really", shell=True)

    # X chromosome doesn't have a recombination map from 1000 genomes
    # so just set LD score to be calculate within an approximate
    # centimorgan window
    else:
        subprocess.check_call("python " + LD_path + "ldsc.py " \
        + "--bfile " + int_directory + chrom \
        + " --ld-wind-kb " + str(window_kb) \
        + " --out " + int_directory + "LD-" + chrom \
        + " --l2 --yes-really", shell=True)

def prune(chrom, int_directory, window, slide, cutoff):
    '''
    Take in all SNPs found on the chromosome and prune them to a set of
    independent variants. PLINK looks at all pairs of variants within an X
    variant window. If a pair has an R^2 correlation greater than some
    threshold, then one SNP within this pair is removed. After pruning within
    this window, PLINK slides along the chromosome by some number of SNPs.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    window : int
        size of window in SNPs over which to compute all pairwise correlations
    slide : int
        number of SNPs to slide over after each iteration
    cutoff : float
        maximum R^2 correlation between two SNPs, above which one SNP in the
        pair is pruned.

    Returns
    -------
    None
    '''
    if chrom != "X":
        subprocess.check_call("plink --bfile " + int_directory  + chrom \
        + " --indep-pairwise " + str(window) + " " + str(slide) + " " + str(cutoff) \
        + " --r" \
        + " --out " + int_directory + chrom, shell=True)

    else:
        subprocess.check_call("plink --bfile " + int_directory  + chrom \
        + " --indep-pairwise " + str(window) + " " + str(slide) + " " + str(cutoff) \
        + " --r" \
        + " --ld-xchr 1" \
        + " --out " + int_directory + chrom, shell=True)
#    subprocess.call("gzip -d " + int_directory + "LD-" + chrom + ".l2.ldscore.gz", shell=True)

def order(chrom, int_directory):
    '''
    Here we reassign LDscore_i_new = (1.0 - LD_score_i / max(LD_score)).
    This gives the highest scoring SNP the lowest p value, and thus the
    greatest priority when clumping.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    LD_file : string
        path of LD score file

    Returns
    -------
    None
    '''
    subprocess.call("gzip -d " + int_directory + "LD-" + chrom + ".l2.ldscore.gz", shell=True)
    # unzip files calculated by LDScore regression
    assoc = open(int_directory + "LD-" + chrom + ".l2.ldscore", 'r')
    # create a new .p file to write results to"
    p_file = open(int_directory + chrom + ".p", 'w')
    # write header
    p_file.write("SNP" + "\t" + 'P\n')
    p_dictionary = {}
    size_dict = {}
    # loop over the association file and for each SNP assign its LDscore
    # in a dictionary
    for line in assoc.readlines()[1:]:
        SNP = line.split('\t')[1]
        LD = float(line.split('\t')[3].rstrip('\n'))
        p_dictionary[SNP] = LD
        assoc.close()
    # calculate the max LDscore for the chromosome
    max_LD = float(max(p_dictionary.values()))
    # for each SNP in the sorted dictionary, map the LDscore
    # to LD_score_new_i = (1.0 - LD_score_i / max(LD_score))
    # Add 1 x 10^-14 so that no SNP has a p value of 0.0
    for SNP in sorted(p_dictionary, key=p_dictionary.get, reverse=True):
        p_val = (1.0 - (float(p_dictionary[SNP])/max_LD)) + 0.00000000000001
        if p_val < 1.0:
            p_file.write(SNP + '\t' + str(p_val) + '\n')
    del p_dictionary
    p_file.close()

def LD_separate(chrom, int_directory):
    '''
    Read in the association files created by order function
    and separate it out into independent(pruned.in) SNPs and
    dependent(pruned.out) SNPs.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored

    Returns
    -------
    None
    '''
    LD_file = open(int_directory + chrom + ".p", 'r')
    prune_file = open(int_directory + chrom + ".prune.in", 'r')
    in_file = open(int_directory + chrom + "_in.p", 'w')
    in_file.write('SNP' + '\t' + 'P\n')
    out_file = open(int_directory + chrom + "_out.p", 'w')
    out_file.write('SNP' + '\t' + 'P\n')
    out_SNPs = []
    index_SNPs = []
    with prune_file as prune:
        for k, line in enumerate(prune):
            index_SNPs.append(line.rstrip('\n'))

    with LD_file as LD:
        for k, line in enumerate(LD):
            if k > 0:
                SNP = line.rstrip('\n').split('\t')[0]
                if SNP in index_SNPs:
                    in_file.write(line)
                else:
                    out_file.write(line)

    in_file.close()
    out_file.close()

def clump(chrom, int_directory, cutoff, max_size):
    '''
    Clumps dependent variants(_out.p) around independent variants(_in.p) in a
    greedy manner.Each SNP "clumped" around an index variant should be within
    a certain window size(specified with clump-kb parameter) and have an r^2
    correlation of atleast "cutoff"(specified with --clump-r2) with the index
    variant. p1 and p2 represent the upper thresholds of signifigance, above
    which SNPs should be excluded. Here they are set to 1 to include all SNPs.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clumped
    max_size : int
        max distance in kb a SNP can be from an index variant in order to be
        included in a clump

    Returns
    -------
    None
    '''
    subprocess.check_call("plink --bfile " + int_directory + chrom \
                + " --clump " + int_directory + chrom + "_in.p " \
                + int_directory + chrom + "_out.p " \
                + "--clump-index-first " \
                + "--clump-p1 1.0 " \
                + "--clump-p2 1.0 " \
                + "--clump-r2 " + str(cutoff) + " " \
                + "--clump-kb " + str(max_size) + " " \
                + "--out " + int_directory + chrom, shell=True)


def reformat_clumps(chrom, int_directory):
    '''
    Script to read in clumps and write them into a haplotype map file readable
    by picard. This file format has recently been deprecated in favor of a VCF
    format file so this should be updated at some point.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clumped
    max_size : int
        max distance in kb a SNP can be from an index variant in order to be
        included in a clump

    Returns
    -------
    None
    '''
    # Sort the clumps by base pair position of index variant
    # and write to a new file
    subprocess.check_call("tail -n+2 " + int_directory + chrom + ".clumped"\
                            + "| sort -k4,4 > " + int_directory + chrom 
			    + "_sorted.clumped", shell=True)

    # create a dictionary of blocks where each key is a SNP
    # and is assigned a value equal to it's index variant
    block_dict = {}
    anchor_file = open(int_directory + chrom + "_sorted.clumped", 'r')
    anchors = []
    for line in anchor_file.readlines()[2:]:
        anchor = line.rstrip('\n').split()[2]
        anchors.append(anchor)
        block_dict[anchor] = ""
        snps = line.rstrip('\n').split()[11].split(',')
        for snp in snps:
            if snp.rstrip('(w)') not in anchors:
                variant = snp.split('(')[0]
                block_dict[variant] = anchor
    anchor_file.close()

    # read in old map
    old_map = open(int_directory + chrom + ".recode_u.vcf", 'r')
    new_map = open(int_directory + chrom + ".map", 'w')
    # exclude header
    # parse each line in old vcf file to extract the SNP
    # name, base pair position, minor and major alleles, MAF
    # and look up its index variant in the block dictionary
    # defined above. If the variant IS an index variant
    # then simply write its index variant as ''
    for line in old_map.readlines()[253:]:
        parsed = line.split('\t')
        for item in parsed[7].split(";"):
            entryname = item.split('=')[0]
            if entryname == "AF":
                AF = item.split('=')[1]
        tab = '\t'
        if parsed[2] in block_dict.keys():
            # order of fields is:
            # chr, pos, name, maj, min,
            # MAF, anchor
            newline = tab.join(parsed[:5]) + '\t' + \
                        AF + '\t' + \
                        block_dict[parsed[2]] + '\t' + \
                        "" + '\t' + '\n'
            new_map.write(newline)
    new_map.close()
    old_map.close()

def detect_negative_LD(chrom, int_directory):
    '''
    Reads in data from clumping procedure and writes out pairs of SNPs in the
    newly created map file that are in high negative LD that give a high r^2
    correlation.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clumped
    max_size : int
        max distance in kb a SNP can be from an index variant in order to be
        included in a clump

    Returns
    -------
    None
    '''

    r_file = open(int_directory + chrom + ".ld", 'r')
    # read in the .ld file created by the clumping procedure
    # and create a dictionary where each key is a SNP pair
    # in tuple format and the value is the r^2 correlation
    r_dict = {}
    with r_file as r:
        for k, line in enumerate(r):
            if k > 0:
                r_dict[(line.split()[2], line.split()[5])] = \
                line.split()[6].rstrip('\n')
    # loop over the newly created map file and for each
    # pair of SNP/anchor SNP, determine whether the
    # the pair's r^2 correlation is negative by checking
    # its value in the dictionary above. If it is, write
    # the pair to a new ".negLD" file

    map_file = open(int_directory + chrom + ".map", 'r')
    with map_file as maps:
        for num, line in enumerate(maps):
    #       if num > (end - 1) :
            if len(line.split('\t')) > 6:
                snp1 = line.split('\t')[2]
                snp2 = line.split('\t')[6]
                if (snp1, snp2) in r_dict.keys():
                    if float(r_dict[(snp1, snp2)]) < -0.50:
                        with open(int_directory + chrom + ".negLD", 'a') as negLDfile:
                            negLDfile.write(snp1 + '\t' + snp2 + '\t' \
                                            + r_dict[(snp1, snp2)] + '\n')
                    del r_dict[(snp1, snp2)]
                if (snp2, snp1) in r_dict.keys():
                    if float(r_dict[(snp2, snp1)]) < -0.50:
                        with open(int_directory + chrom + ".negLD", 'a') as negLDfile:
                            negLDfile.write(snp1 + '\t' + snp2 + '\t' \
                                            + r_dict[(snp2, snp1)] + '\n')
                    del r_dict[(snp2, snp1)]

def switch_alleles(chrom, cwd, int_directory):
    '''
    For each chromosomal map file, read in the pairs of SNPs found to be in
    negative LD and switch the major and minor alleles in the map file. Keep
    original map file and write the new map file to *.filtered.map

    Parameters
    ----------
    chrom : string
        chromosome number
    int_directory : string
        directory path where intermediate and output files are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clumped
    max_size : int
        max distance in kb a SNP can be from an index variant in order to be
        included in a clump

    Returns
    -------
    None
    '''

    tab = '\t'
    negfile = open(int_directory + chrom + ".negLD", 'r')
    neglist = []
    with negfile as neg:
        for line in neg:
            neglist.append(line.split()[0])
    newmapfile =open(cwd + chrom + ".filtered.map", 'w')
    anchor_maf_problems = []
    oldmapfile = open(int_directory + chrom + ".map", 'r')
    with oldmapfile as old:
        for k, line in enumerate(old):
            chrom =line.split('\t')[0]
            bp = line.split('\t')[1]
            maf = float(line.split('\t')[5])
            snp = line.split('\t')[2]
            minor = line.split('\t')[4]
            maj = line.split('\t')[3]
            anch = line.split('\t')[6]
            # switch alleles if found in the list of SNPs in negative LD
            if snp in neglist:
                newline = [chrom, bp, snp, minor, maj, str(maf), anch, "", "\n"]
                newmapfile.write(tab.join(newline))
            else:
                newmapfile.write(line)
    newmapfile.close()

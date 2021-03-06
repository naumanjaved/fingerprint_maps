import subprocess
import os
import itertools as it
import sys
import argparse
import traceback
import time
import vcf

def filter_VCFs(chrom, VCF_file, int_dir, min_MAF):

    '''
    Extract SNPs from each chromosome that are:
    1) biallelic
    2) are phased

    Parameters
    ----------
    chrom : string
        chromosome number
    VCF_file : string
        path to 1000 genomes VCF file
    int_dir : string
        directory where intermediate files are stored
    min_MAF : float
        filter out all SNPs with population averaged MAF less than this value
    Returns
    -------
    None
    '''

    name_prefix = int_dir + "chr_" + chrom
    command = "vcftools --vcf " + VCF_file + " " \
            + "--maf " + str(min_MAF) + " " \
            + "--phased " \
            + "--remove-indels " \
            + "--min-alleles 2 --max-alleles 2 " \
            + "--recode --recode-INFO-all " \
            + "--out " + name_prefix
    subprocess.call(command, shell=True)

def extract_similar_SNPs(chrom, int_dir, SIM, min_MAF):
    '''
    Writes to a new file the SNPs in an input VCF with population specific
    minor allele fractions(MAFs) that do not differ by more than SIM. For
    example,this function would reject a SNP with EUR_AF - AMR_AF  = 0.15 if
    SIM = 0.1.

    Parameters
    ----------
    chrom : string
        chromosome number
    VCF_directory : string
        directory where original VCFs from 1000genomes phase 3 are stored
    int_dir : string
        directory where intermediate files are stored
    SIM : float
        similarity value - the maximum pairwise difference in population
        specific MAF above which a SNP is rejected

    Returns
    -------
    None

    '''
    VCF_file = int_dir + 'chr_' + chrom + '.recode.vcf'
    vcf_reader = vcf.Reader(open(VCF_file, 'r'))  # Loop over VCF records

    sim_SNPs_file = int_dir + "chr_" + chrom + "-common-SNPs.list"
    similar_SNPs = open(sim_SNPs_file, 'w')

    for index, record in enumerate(vcf_reader):
        # store record number

        info = record.INFO

        SNP = str(record.ID)

        if isinstance(SNP,basestring):
            SNP = 'chr' + str(chrom) + ':' + str(index)

        values = []  # List to hold population allele frequencies(AFs)
        differences = []  # List to hold pairwise AF differences

        for field in info.keys():  # Load population AFs
            if '_AF' in field:
                values.append(float(info[field][0]))
                differences = [abs(y-x) for x, y in it.combinations(values, 2)]
        if max(differences) > SIM or len(values) != 5 or min(values) < 0.10:
            continue  # If similarity threshold is exceeded, skip SNP

        similar_SNPs.write(str(index) + '-' + SNP + '\n')
    similar_SNPs.close()

def keep_common_SNPs(chrom, int_dir):
    '''
    Create VCF containing SNPs identified by extract_similar_SNPs fxn

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored
    Returns
    -------
    None
    '''
    sim_SNPs_file = \
            open(int_dir + "chr_" + chrom + "-common-SNPs.list", 'r')
    VCF_file = \
            open(int_dir + "chr_" + chrom + ".recode.vcf", 'r')

    output_file = \
            open(int_dir + "chr_" + chrom + ".common.recode.vcf", 'w')

    SNP_dict = {}

    VCF_lines = VCF_file.readlines()

    VCF_snps = []

    for line in VCF_lines:
        if '#' not in line.split('\t')[0]:
            VCF_snps.append(line)
        else:
            output_file.write(line)
    del VCF_lines

    for line in sim_SNPs_file:
        snp_num, name = line.rstrip('\n').split('-')

        snp_num = int(snp_num)

        SNP_dict[snp_num] = name

    sim_SNPs_file.close()

    for entry in sorted(SNP_dict):

        newline = VCF_snps[int(entry)].split('\t')
        newline[2] = SNP_dict[entry]
        newline = '\t'.join([str(x) for x in newline])
        output_file.write(newline)

    VCF_file.close()
    output_file.close()

def create_PLINK_binary(chrom, int_dir, recomb_file):
    '''
    Creates PLINK binary files from input VCF for use with
    LDSC script(https://github.com/bulik/ldsc).

    Parameters
    ----------
    chrom : string
        chromosome number
    VCFs_directory : string
        directory path where original VCFs from 1000genomes phase 3 are stored
    int_dir : string
        directory path where intermediate files are stored
    recomb_file : string
        1000genomes recombination map file

    Returns
    -------
    None
    '''

    # Use recombination maps to write plink binary files with centimorgans
    command = "plink --vcf " \
                    + int_dir + "chr_" + chrom + ".common.recode.vcf" \
                    + " --cm-map " + recomb_file + " " + chrom \
                    + " --make-bed" \
                    + " --out " + int_dir + chrom

    subprocess.call(command, shell=True)

def LD_score(chrom, int_dir, LD_script, window_cm):
    '''
    Calculates LDscore of all variants for each chromosome
    using LDSC script(https://github.com/bulik/ldsc)

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored
    LD_script : string
        directory containing python LD script
    window_cm : float
        for autosomes, the window in centimorgans over which to calculate LD
        scores

    Returns
    -------
    None
    '''

    command = "python " + LD_script \
                 + " --bfile " + int_dir + chrom \
                 + " --ld-wind-cm " + str(window_cm) \
                 + " --out " + int_dir + "LD-" + chrom \
                 + " --l2 --yes-really"

    subprocess.call(command, shell=True)

def prune(chrom, int_dir, window, slide, cutoff):
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
    int_dir : string
        directory path where intermediate files are stored
    window : int
        size of window in kb over which to compute all pairwise correlations
    slide : int
        number of SNPs to slide over after each iteration
    cutoff : float
        maximum R^2 correlation between two SNPs, above which one SNP in the
        pair is pruned.

    Returns
    -------
    None
    '''
    prune_params = str(window) + " " \
                   + str(slide) + " " \
                   + str(cutoff)

    command = "plink --bfile " + int_dir + chrom \
                    + " --indep-pairwise " + prune_params \
                    + " --r" \
                    + " --out " + int_dir + chrom
    if chrom == 'X':
        command = "plink --bfile " + int_dir + chrom \
                            + " --indep-pairwise " + prune_params \
                            + " --r" \
                            + " --ld-xchr 1" \
                            + " --out " + int_dir + chrom

    subprocess.call(command, shell=True)

def order(chrom, int_dir):
    '''
    Here we reassign LDscore_i_new = (1.0 - LD_score_i / max(LD_score)).
    This gives the highest scoring SNP the lowest p value, and thus the
    greatest priority when clumping.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored
    LD_file : string
        path of LD score file

    Returns
    -------
    None
    '''

    association_file_name = int_dir + "LD-" + chrom + ".l2.ldscore"
    # unzip LD score files
    subprocess.call("gzip -d " + association_file_name + ".gz", shell=True)

    assoc = open(association_file_name, 'r')

    # write a new association file to rank SNPs by LDscore
    p_file = open(int_dir + chrom + ".p", 'w')

    p_file.write("SNP" + "\t" + 'P\n')
    p_dictionary = {}  # Create a SNP:LDscore dictionary
    size_dict = {}

    for line in assoc.readlines()[1:]:
        SNP = line.split('\t')[1]
        LD = float(line.split('\t')[3].rstrip('\n'))
        p_dictionary[SNP] = LD
        assoc.close()


    # Calculate the max LDscore for the chromosome
    max_LD = float(max(p_dictionary.values()))
    '''
    For each SNP in the sorted dictionary, map the LDscore
    to LD_score_new_i = (1.0 - LD_score_i / max(LD_score))
    Add 1e-6 to set the minimum "signifigance value"
    '''
    for SNP in sorted(p_dictionary, key=p_dictionary.get, reverse=True):

        p_val = (1.0 - float(p_dictionary[SNP])/(max_LD)) + 0.000001
        if p_val < 1.0:
            p_file.write(SNP + '\t' + str(p_val) + '\n')

    p_file.close()

def LD_separate(chrom, int_dir):

    '''
    Read in the association files created by order function
    and separate it out into independent(pruned.in) SNPs and
    dependent(pruned.out) SNPs.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored

    Returns
    -------
    None
    '''

    LD_file = open(int_dir + chrom + ".p", 'r')
    prune_file = open(int_dir + chrom + ".prune.in", 'r')
    in_file = open(int_dir + chrom + "_in.p", 'w')
    in_file.write('SNP' + '\t' + 'P\n')
    out_file = open(int_dir + chrom + "_out.p", 'w')
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

def clump(chrom, int_dir, cutoff, max_size):
    '''
    Clumps dependent variants(_out.p) around independent variants(_in.p)..Each
    SNP "clumped" around an index variant should be within a certain window
    size(specified with clump-kb parameter) and have an r^2 correlation of
    atleast "cutoff"(specified with --clump-r2) with the index variant. p1 and p2
    represent the upper thresholds of signifigance, above which SNPs should be
    excluded. Here they are set to 1 to include all SNPs.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clumped
    max_size : int
        max distance in kb a SNP can be from an index variant in order to be
        included in a clump

    Returns
    -------
    None
    '''
    command = "plink --bfile " + int_dir + chrom \
                     + " --clump " + int_dir + chrom + "_in.p " \
                     + int_dir + chrom + "_out.p " \
                     + "--clump-index-first " \
                     + "--clump-p1 1.0 " \
                     + "--clump-p2 1.0 " \
                     + "--clump-r2 " + str(cutoff) + " " \
                     + "--clump-kb " + str(max_size) + " " \
                     + "--out " + int_dir + chrom

    subprocess.call(command, shell=True)

def reformat_clumps(chrom, int_dir):
    '''
    Script to read in clumps and write them into a haplotype map file readable
    by picard. This file format has recently been deprecated in favor of a VCF
    format file so this should be updated at some point.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored

    Returns
    -------
    None
    '''
    # Sort the clumps by base pair position of index variant
    command = "tail -n+2 " + int_dir + chrom + ".clumped" \
                            + " | sort -k4,4 > " + int_dir + chrom \
                            + "_sorted.clumped"
    subprocess.call(command, shell=True)

    block_dict = {}  # Block dictionary with variant:index variant pairs
    anchor_file = open(int_dir + chrom + "_sorted.clumped", 'r')
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


    bim_dict = {}
    # Read in plink .bim file to recover major and minor alleles
    with open(int_dir + chrom + '.bim') as bim_file:
        for line in bim_file:
            line = line.rstrip('\n')
            entries = line.split('\t')
            variant =  entries[1]
            minor = entries[4]
            major = entries[5]
            bim_dict[variant] = (minor, major)

    old_map = vcf.Reader(open(int_dir + "chr_" + chrom + ".common.recode.vcf", 'r'))

    new_map = open(int_dir + chrom + ".map", 'w')
    '''
    parse each line in old vcf file to extract the SNP name, base
    pair position, minor and major alleles, MAF and look up its index
    variant in the block dictionary defined above. If the variant IS
    an index variant then write its index variant as ''
    '''
    for record in old_map:
        if record.ID in block_dict.keys():
            ref = record.REF[0]
            alt = record.ALT[0]
            maf = record.INFO['AF'][0]
            minor, major = bim_dict[record.ID]
            switched_maf = str(maf)
            if str(alt) == str(major):
                switched_maf = str(1.0 - float(maf))
            ref = major
            alt = minor

            newline = str(record.CHROM) + '\t' \
                      + str(record.POS) + '\t' \
                      + str(record.ID) + '\t' \
                      + str(ref) + '\t' \
                      + str(alt) + '\t' \
                      + str(switched_maf) + '\t' \
                      + str(block_dict[record.ID]) \
                      + '\t' + "" + '\t' + '\n'
            new_map.write(newline)
    new_map.close()

def detect_negative_LD(chrom, int_dir):
    '''
    Reads in data from clumping procedure and writes out pairs of SNPs in the
    newly created map file that are in high negative LD that give a high r^2
    correlation.

    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate and output files are stored

    Returns
    -------
    None

    '''
    r_file = open(int_dir + chrom + ".ld", 'r')
    r_dict = {}  # Create a dictionary of SNPs:r^2 pairs
    with r_file as r:
        for k, line in enumerate(r):
            if k > 0:
                r_dict[(line.split()[2], line.split()[5])] = \
                        line.split()[6].rstrip('\n')
    '''
    Loop over the newly created map file and for each
    pair of SNP/anchor SNP, determine whether the
    the pair's r^2 correlation is negative by checking
    its value in the dictionary above. If it is, write
    the pair to a new ".negLD" file
    '''

    map_file = open(int_dir + chrom + ".map", 'r')
    neg_LD_file_name = int_dir + chrom + ".negLD"
    with map_file as maps:
        for num, line in enumerate(maps):
            if len(line.split('\t')) > 6:
                snp1 = line.split('\t')[2]
                snp2 = line.split('\t')[6]
                if (snp1, snp2) in r_dict.keys():
                    if float(r_dict[(snp1, snp2)]) < 0.0:
                        with open(neg_LD_file_name, 'a') as negLDfile:
                            negLDfile.write(snp1 + '\t' + snp2 + '\t'
                                            + r_dict[(snp1, snp2)] + '\n')
                            del r_dict[(snp1, snp2)]
                if (snp2, snp1) in r_dict.keys():
                    if float(r_dict[(snp2, snp1)]) < 0.0:
                        with open(neg_LD_file_name, 'a') as negLDfile:
                            negLDfile.write(snp1 + '\t' + snp2 + '\t'
                                            + r_dict[(snp2, snp1)] + '\n')
                            del r_dict[(snp2, snp1)]

def switch_alleles(chrom, int_dir, out_dir):
    '''
    For each chromosomal map file, read in the pairs of SNPs found to be in
    negative LD, switch the major and minor alleles in the map file, and
    switch the MAF to 1-MAF.Keep original map file and write the new map file
    to *.filtered.map
    Parameters
    ----------
    chrom : string
        chromosome number
    int_dir : string
        directory path where intermediate files are stored

    Returns
    -------
    None

    '''

    # If no pairs were in negative linkage, copy to output folder and skip
    if not os.path.isfile(int_dir + chrom + ".negLD"):
        command = "cp " + int_dir + chrom + ".map " \
                + out_dir + chrom + ".filtered.map"
        subprocess.call(command, shell=True)
        return None

    tab = '\t'
    negfile = open(int_dir + chrom + ".negLD", 'r')
    neglist = []
    with negfile as neg:
        for line in neg:
            neglist.append(line.split()[0])
    newmapfile =open(out_dir + chrom + ".filtered.map", 'w')
    oldmapfile = open(int_dir + chrom + ".map", 'r')
    with oldmapfile as old:
        for k, line in enumerate(old):
            chrom =line.split('\t')[0]
            bp = line.split('\t')[1]
            maf = line.split('\t')[5]
            snp = line.split('\t')[2]
            minor = line.split('\t')[4]
            maj = line.split('\t')[3]
            anch = line.split('\t')[6]
            # Switch alleles and replace maf with 1-maf if found in the list of SNPs in negative LD
            if snp in neglist:
                newline = [chrom, bp, snp, minor, maj, str(1.0 - float(maf)), anch, "", "\n"]
                newmapfile.write(tab.join(newline))
            else:
                newmapfile.write(line)
    newmapfile.close()

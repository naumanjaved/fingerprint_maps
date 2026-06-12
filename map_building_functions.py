import subprocess
import sys
import os
import gzip
import shutil
import itertools as it
import vcf


def filter_VCFs(chrom, VCF_file, int_dir, min_MAF):
    '''
    Extract SNPs from each chromosome that are:
    1) biallelic
    2) are phased

    Parameters
    ----------
    chrom : str
        chromosome number
    VCF_file : str
        path to 1000 genomes VCF file
    int_dir : str
        directory where intermediate files are stored
    min_MAF : float
        filter out all SNPs with population averaged MAF less than this value
    '''
    output_vcf = os.path.join(int_dir, "chr_" + chrom + ".recode.vcf")
    command = ("bcftools view " + VCF_file + " "
               + "--types snps "
               + "--min-alleles 2 --max-alleles 2 "
               + "--phased "
               + "-e 'MAF < " + str(min_MAF) + "' "
               + "-Ov -o " + output_vcf)
    subprocess.run(command, shell=True, check=True)


def extract_similar_SNPs(chrom, int_dir, SIM, min_MAF):
    '''
    Writes to a new file the SNPs in an input VCF with population specific
    minor allele fractions(MAFs) that do not differ by more than SIM.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory where intermediate files are stored
    SIM : float
        similarity value - the maximum pairwise difference in population
        specific MAF above which a SNP is rejected
    min_MAF : float
        minimum minor allele frequency threshold
    '''
    vcf_path = os.path.join(int_dir, 'chr_' + chrom + '.recode.vcf')
    sim_snps_path = os.path.join(int_dir, "chr_" + chrom + "-common-SNPs.list")

    vcf_reader = vcf.Reader(filename=vcf_path)

    with open(sim_snps_path, 'w') as similar_SNPs:
        for index, record in enumerate(vcf_reader):
            info = record.INFO

            if record.ID is None:
                SNP = 'chr' + str(chrom) + ':' + str(index)
            else:
                SNP = str(record.ID)

            values = []
            for field in info:
                if '_AF' in field:
                    values.append(float(info[field][0]))

            if len(values) != 5:
                continue

            differences = [abs(y - x) for x, y in it.combinations(values, 2)]

            if max(differences) > SIM or min(values) < min_MAF:
                continue

            similar_SNPs.write(str(index) + '-' + SNP + '\n')


def keep_common_SNPs(chrom, int_dir):
    '''
    Create VCF containing SNPs identified by extract_similar_SNPs fxn

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    '''
    sim_snps_path = os.path.join(int_dir, "chr_" + chrom + "-common-SNPs.list")
    vcf_path = os.path.join(int_dir, "chr_" + chrom + ".recode.vcf")
    output_path = os.path.join(int_dir, "chr_" + chrom + ".common.recode.vcf")

    SNP_dict = {}
    with open(sim_snps_path, 'r') as sim_SNPs_file:
        for line in sim_SNPs_file:
            snp_num, name = line.rstrip('\n').split('-', 1)
            SNP_dict[int(snp_num)] = name

    VCF_snps = []
    with open(vcf_path, 'r') as VCF_file, \
         open(output_path, 'w') as output_file:
        for line in VCF_file:
            if '#' in line.split('\t')[0]:
                output_file.write(line)
            else:
                VCF_snps.append(line)

        for entry in sorted(SNP_dict):
            newline = VCF_snps[entry].split('\t')
            newline[2] = SNP_dict[entry]
            output_file.write('\t'.join(newline))


def create_PLINK_binary(chrom, int_dir, recomb_file):
    '''
    Creates PLINK binary files from input VCF for use with
    LDSC script(https://github.com/bulik/ldsc).

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    recomb_file : str
        1000genomes recombination map file
    '''
    vcf_path = os.path.join(int_dir, "chr_" + chrom + ".common.recode.vcf")
    out_prefix = os.path.join(int_dir, chrom)

    command = ("plink --vcf " + vcf_path
               + " --cm-map " + recomb_file + " " + chrom
               + " --make-bed"
               + " --out " + out_prefix)
    subprocess.run(command, shell=True, check=True)


def LD_score(chrom, int_dir, LD_script, window_cm):
    '''
    Calculates LDscore of all variants for each chromosome
    using LDSC script(https://github.com/bulik/ldsc)

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    LD_script : str
        path to ldsc.py
    window_cm : float
        for autosomes, the window in centimorgans over which to calculate LD scores
    '''
    bfile = os.path.join(int_dir, chrom)
    out_prefix = os.path.join(int_dir, "LD-" + chrom)

    command = (sys.executable + " " + LD_script
               + " --bfile " + bfile
               + " --ld-wind-cm " + str(window_cm)
               + " --out " + out_prefix
               + " --l2 --yes-really")
    subprocess.run(command, shell=True, check=True)


def prune(chrom, int_dir, window, slide, cutoff):
    '''
    Prune SNPs to a set of independent variants. PLINK looks at all pairs of
    variants within a window. If a pair has an R^2 correlation greater than
    some threshold, one SNP is removed. After pruning, PLINK slides along the
    chromosome by some number of SNPs.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    window : int
        size of window in kb over which to compute all pairwise correlations
    slide : int
        number of SNPs to slide over after each iteration
    cutoff : float
        maximum R^2 correlation between two SNPs, above which one is pruned
    '''
    bfile = os.path.join(int_dir, chrom)
    prune_params = str(window) + " " + str(slide) + " " + str(cutoff)

    command = ("plink --bfile " + bfile
               + " --indep-pairwise " + prune_params
               + " --r"
               + " --out " + bfile)
    if chrom == 'X':
        command += " --ld-xchr 1"

    subprocess.run(command, shell=True, check=True)


def order(chrom, int_dir):
    '''
    Reassign LDscore_i_new = (1.0 - LD_score_i / max(LD_score)).
    This gives the highest scoring SNP the lowest p value, and thus the
    greatest priority when clumping.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    '''
    gz_path = os.path.join(int_dir, "LD-" + chrom + ".l2.ldscore.gz")

    p_dictionary = {}
    with gzip.open(gz_path, 'rt') as assoc:
        next(assoc)
        for line in assoc:
            fields = line.split('\t')
            SNP = fields[1]
            LD = float(fields[3].rstrip('\n'))
            p_dictionary[SNP] = LD

    max_LD = max(p_dictionary.values())

    p_file_path = os.path.join(int_dir, chrom + ".p")
    with open(p_file_path, 'w') as p_file:
        p_file.write("SNP\tP\n")
        for SNP in sorted(p_dictionary, key=p_dictionary.get, reverse=True):
            p_val = (1.0 - p_dictionary[SNP] / max_LD) + 0.000001
            if p_val < 1.0:
                p_file.write(SNP + '\t' + str(p_val) + '\n')


def LD_separate(chrom, int_dir):
    '''
    Read in the association files created by order function
    and separate into independent(pruned.in) and dependent(pruned.out) SNPs.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    '''
    ld_path = os.path.join(int_dir, chrom + ".p")
    prune_path = os.path.join(int_dir, chrom + ".prune.in")
    in_path = os.path.join(int_dir, chrom + "_in.p")
    out_path = os.path.join(int_dir, chrom + "_out.p")

    with open(prune_path, 'r') as prune_file:
        index_SNPs = set(line.rstrip('\n') for line in prune_file)

    with open(ld_path, 'r') as LD_file, \
         open(in_path, 'w') as in_file, \
         open(out_path, 'w') as out_file:

        in_file.write('SNP\tP\n')
        out_file.write('SNP\tP\n')

        next(LD_file)
        for line in LD_file:
            SNP = line.rstrip('\n').split('\t')[0]
            if SNP in index_SNPs:
                in_file.write(line)
            else:
                out_file.write(line)


def clump(chrom, int_dir, cutoff, max_size):
    '''
    Clumps dependent variants around independent variants. Each SNP "clumped"
    around an index variant should be within a certain window size and have an
    r^2 correlation of at least cutoff with the index variant.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediates are stored
    cutoff : float
        minimum r^2 correlation required to become part of a clump
    max_size : int
        max distance in kb a SNP can be from an index variant
    '''
    bfile = os.path.join(int_dir, chrom)
    in_p = os.path.join(int_dir, chrom + "_in.p")
    out_p = os.path.join(int_dir, chrom + "_out.p")

    command = ("plink --bfile " + bfile
               + " --clump " + in_p + " " + out_p + " "
               + "--clump-index-first "
               + "--clump-p1 1.0 "
               + "--clump-p2 1.0 "
               + "--clump-r2 " + str(cutoff) + " "
               + "--clump-kb " + str(max_size) + " "
               + "--out " + bfile)
    subprocess.run(command, shell=True, check=True)


def reformat_clumps(chrom, int_dir):
    '''
    Read in clumps and write them into a haplotype map file readable by picard.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    '''
    clumped_path = os.path.join(int_dir, chrom + ".clumped")
    sorted_path = os.path.join(int_dir, chrom + "_sorted.clumped")

    if not os.path.isfile(clumped_path):  # plink writes no .clumped when there are no input snps
        raise ValueError('no clumps produced for chrom ' + chrom
                         + ' - check that ' + chrom + '_in.p contains snps')

    command = "tail -n+2 " + clumped_path + " | sort -k4,4 > " + sorted_path
    subprocess.run(command, shell=True, check=True)

    block_dict = {}
    anchors = []
    with open(sorted_path, 'r') as anchor_file:
        for line in anchor_file.readlines()[2:]:
            anchor = line.rstrip('\n').split()[2]
            anchors.append(anchor)
            block_dict[anchor] = ""
            snps = line.rstrip('\n').split()[11].split(',')
            for snp in snps:
                if snp.rstrip('(w)') not in anchors:
                    variant = snp.split('(')[0]
                    block_dict[variant] = anchor

    bim_dict = {}
    bim_path = os.path.join(int_dir, chrom + '.bim')
    with open(bim_path) as bim_file:
        for line in bim_file:
            entries = line.rstrip('\n').split('\t')
            variant = entries[1]
            minor = entries[4]
            major = entries[5]
            bim_dict[variant] = (minor, major)

    vcf_path = os.path.join(int_dir, "chr_" + chrom + ".common.recode.vcf")
    old_map = vcf.Reader(filename=vcf_path)

    map_path = os.path.join(int_dir, chrom + ".map")
    with open(map_path, 'w') as new_map:
        for record in old_map:
            if record.ID in block_dict:
                alt = record.ALT[0]
                maf = record.INFO['AF'][0]
                minor, major = bim_dict[record.ID]
                switched_maf = str(maf)
                if str(alt) == str(major):
                    switched_maf = str(1.0 - float(maf))

                newline = '\t'.join([
                    str(record.CHROM), str(record.POS), str(record.ID),
                    str(major), str(minor), switched_maf,
                    str(block_dict[record.ID]), "", ""
                ]) + '\n'
                new_map.write(newline)


def detect_negative_LD(chrom, int_dir):
    '''
    Reads in data from clumping procedure and writes out pairs of SNPs
    in negative LD that give a high r^2 correlation.

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate and output files are stored
    '''
    ld_path = os.path.join(int_dir, chrom + ".ld")
    r_dict = {}
    with open(ld_path, 'r') as r_file:
        for k, line in enumerate(r_file):
            if k > 0:
                fields = line.split()
                r_dict[(fields[2], fields[5])] = fields[6].rstrip('\n')

    map_path = os.path.join(int_dir, chrom + ".map")
    neg_LD_path = os.path.join(int_dir, chrom + ".negLD")

    neg_pairs = []
    with open(map_path, 'r') as map_file:
        for line in map_file:
            parts = line.split('\t')
            if len(parts) > 6:
                snp1 = parts[2]
                snp2 = parts[6]
                for pair in [(snp1, snp2), (snp2, snp1)]:
                    if pair in r_dict and float(r_dict[pair]) < 0.0:
                        neg_pairs.append(
                            snp1 + '\t' + snp2 + '\t' + r_dict[pair] + '\n')
                        del r_dict[pair]

    if neg_pairs:
        with open(neg_LD_path, 'w') as neg_LD_file:
            for pair_line in neg_pairs:
                neg_LD_file.write(pair_line)


def switch_alleles(chrom, int_dir, out_dir):
    '''
    For each chromosomal map file, read in the pairs of SNPs found to be in
    negative LD, switch the major and minor alleles in the map file, and
    switch the MAF to 1-MAF. Keep original map file and write the new map file
    to *.filtered.map

    Parameters
    ----------
    chrom : str
        chromosome number
    int_dir : str
        directory path where intermediate files are stored
    out_dir : str
        directory path where output files are stored
    '''
    map_path = os.path.join(int_dir, chrom + ".map")
    output_path = os.path.join(out_dir, chrom + ".filtered.map")
    neg_ld_path = os.path.join(int_dir, chrom + ".negLD")

    if not os.path.isfile(neg_ld_path):
        shutil.copy2(map_path, output_path)
        return

    neg_snps = set()
    with open(neg_ld_path, 'r') as negfile:
        for line in negfile:
            neg_snps.add(line.split()[0])

    with open(map_path, 'r') as oldmapfile, \
         open(output_path, 'w') as newmapfile:
        for line in oldmapfile:
            parts = line.split('\t')
            snp = parts[2]
            if snp in neg_snps:
                chr_val = parts[0]
                bp = parts[1]
                maf = parts[5]
                minor = parts[4]
                maj = parts[3]
                anch = parts[6]
                newline = '\t'.join([chr_val, bp, snp, minor, maj,
                                     str(1.0 - float(maf)), anch, "", "\n"])
                newmapfile.write(newline)
            else:
                newmapfile.write(line)

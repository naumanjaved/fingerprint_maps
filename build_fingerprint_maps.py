"""

@author: Nauman Javed

"""
import map_building_functions as build
import os
import argparse
import traceback
import time
import datetime
from multiprocessing import Pool

__version__ = '1.0.1'
MASTHEAD = "---------------------------------------------------------------\n"
MASTHEAD += "* Fingerprint_Maps\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "--------------------------------------------------------------\n"


class Logger:
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        self.log_fh.write(msg + '\n')
        self.log_fh.flush()
        print(msg)


def create_maps(chrom, vcf_file, recomb_file, args, log, int_dir, out_dir):

    log.log('Pre-processing VCFs.')
    build.filter_VCFs(chrom, vcf_file, int_dir, args.min_MAF)
    log.log('Finished processing VCFs.')

    log.log('Creating list of SNPs with similar MAFs across populations...')
    build.extract_similar_SNPs(chrom, int_dir, args.similarity, args.min_MAF)
    log.log('Finished creating lists of similar SNPs...')

    log.log('Creating VCFs with common variants')
    build.keep_common_SNPs(chrom, int_dir)
    log.log('Finished creating VCFs with common variants')

    log.log('Creating PLINK binary files...')
    build.create_PLINK_binary(chrom, int_dir, recomb_file)
    log.log('Finished creating binary files.')

    log.log('Calculating LDScores...')
    build.LD_score(chrom, int_dir, args.LD_script,
                   args.LDScore_window_autosome)
    log.log('Finished calculating LDScores.')

    log.log('Creating PLINK association files...')
    build.order(chrom, int_dir)
    log.log('Finished creating PLINK association files.')

    log.log('Pruning SNPs...')
    build.prune(chrom, int_dir,
                args.prune_window, args.prune_slide, args.prune_cutoff)
    log.log('Finished pruning SNPs.')

    log.log('Separating out LDScores into dependent and independent SNP files...')
    build.LD_separate(chrom, int_dir)
    log.log('Finished separating LDScores.')

    log.log('Clumping SNPs...')
    build.clump(chrom, int_dir,
                args.clump_cutoff, args.max_distance_clump)
    log.log('Finished clumping SNPs.')

    log.log('Building map file...')
    build.reformat_clumps(chrom, int_dir)
    log.log('Finished building map files.')

    log.log('Detecting negative LD...')
    build.detect_negative_LD(chrom, int_dir)
    log.log('Finished recording negative LD.')

    log.log('Switching alleles...')
    build.switch_alleles(chrom, int_dir, out_dir)
    log.log('Finished switching negative LD alleles.')


def run_chromosome(chrom, args, int_dir, out_dir):
    vcf_file = args.VCF_file.format(chrom=chrom)
    recomb_file = args.recomb_file.format(chrom=chrom)

    if not os.path.isfile(vcf_file):
        raise ValueError('VCF file not found for chrom ' + chrom + ': ' + vcf_file)
    if not os.path.isfile(recomb_file):
        raise ValueError('recomb file not found for chrom ' + chrom + ': ' + recomb_file)

    log = Logger(os.path.join(out_dir, chrom + '.log'))

    start_time = time.time()
    try:
        opts = vars(args)
        header = MASTHEAD
        header += 'build_fingerprint_maps.py (chrom ' + chrom + ') \\\n'
        options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\'
                   for x in opts]
        header += '\n'.join(options).replace('True', '').replace('False', '')
        header = header[0:-1] + '\n'
        log.log(header)
        log.log('Beginning analysis of chromosome {C} at {T}'.format(
            C=chrom, T=time.ctime()))
        create_maps(chrom, vcf_file, recomb_file, args, log, int_dir, out_dir)
    except Exception:
        log.log(traceback.format_exc())
        raise
    finally:
        log.log('Analysis of chromosome {C} finished at {T}'.format(
            C=chrom, T=time.ctime()))
        time_elapsed = round(time.time() - start_time, 2)
        log.log('Total time elapsed: {T}'.format(
            T=str(datetime.timedelta(seconds=time_elapsed))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    accepted_chromosomes = [str(i) for i in range(1, 23)] + ['X']

    parser.add_argument('--int_dir', required=True, type=str,
        help='directory to store intermediate files')
    parser.add_argument('--out_dir', required=True, type=str,
        help='directory to store output files')
    parser.add_argument('--recomb_file', required=True, type=str,
        help='Shapeit recombination file (use {chrom} placeholder for multi-chrom runs)')
    parser.add_argument('--chrom', required=True, type=str, nargs='+',
        help='chromosome(s) to process (1-22, X, or "all")')
    parser.add_argument('--VCF_file', required=True, type=str,
        help='VCF file path (use {chrom} placeholder for multi-chrom runs)')
    parser.add_argument('--LD_script', required=True, type=str,
        help='path to ldsc.py')

    parser.add_argument('--similarity', default=0.10, type=float,
        help='maximum difference in population specific MAF allowed')
    parser.add_argument('--min_MAF', default=0.10, type=float,
        help='minimum minor allele fraction of SNPs included in map')
    parser.add_argument('--LDScore_window_autosome', default=1.0, type=float,
        help='window size (in centimorgans) over which to calculate LDScore')
    parser.add_argument('--prune_window', default=10000, type=int,
        help='window size in kb for PLINK prune function')
    parser.add_argument('--prune_slide', default=5, type=int,
        help='number of SNPs to slide over on each iteration of PLINK prune')
    parser.add_argument('--prune_cutoff', default=0.10, type=float,
        help='maximum r^2 correlation allowed between pruned SNPs')
    parser.add_argument('--clump_cutoff', default=0.85, type=float,
        help='minimum r^2 correlation required for SNPs to be clumped together')
    parser.add_argument('--max_distance_clump', default=10000, type=int,
        help='maximum distance in kb a SNP can be from index SNP when forming clump')
    parser.add_argument('--jobs', default=1, type=int,
        help='number of chromosomes to process in parallel (default: 1)')

    args = parser.parse_args()

    if not os.path.isfile(args.LD_script):
        raise ValueError('--LD_script not found: ' + args.LD_script)

    if not (0.0 <= args.similarity <= 1.0):
        raise ValueError('--similarity must be between 0.0 and 1.0')
    if not (0.0 <= args.min_MAF <= 1.0):
        raise ValueError('--min_MAF must be between 0.0 and 1.0')
    if args.LDScore_window_autosome < 0.0:
        raise ValueError('--LDScore_window_autosome must be >= 0.0')
    if args.prune_window <= 0:
        raise ValueError('--prune_window must be > 0')
    if args.prune_slide <= 0:
        raise ValueError('--prune_slide must be > 0')
    if not (0.0 <= args.prune_cutoff <= 1.0):
        raise ValueError('--prune_cutoff must be between 0.0 and 1.0')
    if not (0.0 <= args.clump_cutoff <= 1.0):
        raise ValueError('--clump_cutoff must be between 0.0 and 1.0')
    if args.max_distance_clump <= 0:
        raise ValueError('--max_distance_clump must be > 0')
    if args.jobs < 1:
        raise ValueError('--jobs must be >= 1')

    if args.chrom == ['all']:
        chromosomes = accepted_chromosomes
    else:
        chromosomes = args.chrom
        for c in chromosomes:
            if c not in accepted_chromosomes:
                raise ValueError('--chrom values must be 1-22, X, or "all"')

    if len(chromosomes) > 1:
        if '{chrom}' not in args.VCF_file:
            raise ValueError(
                '--VCF_file must contain {chrom} placeholder for multi-chrom runs')
        if '{chrom}' not in args.recomb_file:
            raise ValueError(
                '--recomb_file must contain {chrom} placeholder for multi-chrom runs')

    int_dir = args.int_dir
    os.makedirs(int_dir, exist_ok=True)

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    if len(chromosomes) == 1 or args.jobs == 1:
        for chrom in chromosomes:
            run_chromosome(chrom, args, int_dir, out_dir)
    else:
        with Pool(args.jobs) as pool:
            pool.starmap(run_chromosome,
                         [(c, args, int_dir, out_dir) for c in chromosomes])

"""

@author: Nauman Javed

"""
from __future__ import division
import map_building_functions as build
import subprocess
import os
import itertools as it
#import numpy as np
import sys
import argparse
import traceback
import time
import datetime

__version__ = '1.0.1'
MASTHEAD = "---------------------------------------------------------------\n"
MASTHEAD += "* Fingerprint_Maps\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "--------------------------------------------------------------\n"

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print >>self.log_fh, msg
        print msg

def create_maps(args, log):



    log.log('Pre-processing VCFs.')
    build.filter_VCFs(args.chrom, args.VCF_file, int_dir, args.min_MAF)
    log.log('Finished processing VCFs.')

    log.log('Creating list of SNPs with similar MAFs across populations...')
    build.extract_similar_SNPs(args.chrom,
            int_dir, args.similarity, args.min_MAF)
    log.log('Finished creating lists of similar SNPs...')

    log.log('Creating VCFs with common variants')
    build.keep_common_SNPs(args.chrom, int_dir)
    log.log('Finished creating VCFs with common variants')

    log.log('Creating PLINK binary files...')
    build.create_PLINK_binary(str(args.chrom), int_dir,
                              args.recomb_file)
    log.log('Finished creating binary files.')

    log.log('Calculating LDScores...')
    build.LD_score(str(args.chrom), int_dir, args.LD_script,
                    args.LDScore_window_autosome)
    log.log('Finished calculating LDScores.')

    log.log('Creating PLINK association files...')
    build.order(args.chrom, int_dir)
    log.log('Finished creating PLINK association files.')

    log.log('Pruning SNPs...')
    build.prune(args.chrom, int_dir,
                int(args.prune_window), args.prune_slide, args.prune_cutoff)
    log.log('Finished pruning SNPs.')

    log.log('Separating out LDScores into dependent and independent SNP files...')
    build.LD_separate(args.chrom, int_dir)
    log.log('Finished separating LDScores.')

    log.log('Clumping SNPs...')
    build.clump(args.chrom, int_dir,
               args.clump_cutoff, args.max_distance_clump)
    log.log('Finished clumping SNPs.')

    log.log('Building map file...')
    build.reformat_clumps(args.chrom, int_dir)
    log.log('Finished building map files.')

    log.log('Detecting negative LD...')
    build.detect_negative_LD(args.chrom, int_dir)
    log.log('Finished recording negative LD.')

    log.log('Switching alleles...')
    build.switch_alleles(args.chrom, int_dir, out_dir)
    log.log('Finished switching negative LD alleles.')


# Argument parsing
parser = argparse.ArgumentParser()

# Directory and reference file specifications'
parser.add_argument('--int_dir', default=None, type=str,
    help='directory to store intermediate files')
parser.add_argument('--out_dir', default=None, type=str,
    help='directory to store output files')
parser.add_argument('--recomb_file', default=None, type=str,
    help='Shapeit recombination file')
parser.add_argument('--chrom', default=None, type=str,
    help='Chromosome for which to calculate.')
parser.add_argument('--VCF_file', default=None, type=str,
    help='VCF file name for the selected chromosome')
parser.add_argument('--LD_script', default=None, type=str,
    help='Directory containing ldsc.py ')

# SNP filtering and map calculation parameters
parser.add_argument('--similarity', default=0.10, type=float,
    help='Maximum difference in population specific MAF allowed')
parser.add_argument('--min_MAF', default=0.10, type=float,
    help='Minimum minor allele fraction of SNPs that will be included in map')
parser.add_argument('--LDScore_window_autosome', default=1.0, type=float,
    help='Window size(in centimorgans) over which to calculate LDScore.')
parser.add_argument('--prune_window', default=10000, type=int,
    help='Window size in kb for PLINK prune function')
parser.add_argument('--prune_slide', default=5, type=int,
    help='Number of SNPs to slide over on each iteration of PLINK prune')
parser.add_argument('--prune_cutoff', default=0.10, type=float,
    help='Maximum r^2 correlation allowed between pruned SNPs')
parser.add_argument('--clump_cutoff', default=0.85, type=float,
    help='Minimum r^2 correlation required for SNPs to be clumped together')
parser.add_argument('--max_distance_clump', default=10000, type=int,
    help='Maximum distance in kb a SNP can be from index SNP when forming clump')

if __name__ == "__main__":
    args = parser.parse_args()
    accepted_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8',
                            '9', '10', '11', '12', '13', '14', '15', '16',
                            '17', '18', '19', '20', '21', '22', 'X']

    # argument error handling

    if args.int_dir is None:
        raise ValueError('--int_dir is required.')

    if args.out_dir is None:
        raise ValueError('--out_dir is required.')

    if args.recomb_file is None:
        raise ValueError('--recomb_file is required.')
    if not os.path.isfile(args.recomb_file):
        raise ValueError('--recomb_file not found.')

    if args.chrom is None :
        raise ValueError('--chrom is required.')
    if args.chrom not in accepted_chromosomes:
        raise ValueError('--chrom must be between 1-22 or X.')

    if args.VCF_file is None:
        raise ValueError('--VCF_file is required.')
    if not os.path.isfile(args.VCF_file):
        raise ValueError('--VCF_file not found.')

    if args.LD_script is None:
        raise ValueError('--LD_script is required.')
    if not os.path.isfile(args.LD_script):
        raise ValueError('--LD_script not found.')

    if args.similarity is None:
        raise ValueError('--similarity is required.')

    if (args.similarity < 0.0) or (args.similarity > 1.0):
        raise ValueError('--similarity must be float > 0.0 and < 1.0')

    if args.min_MAF is None:
        raise ValueError('--min_MAF is required.')

    if (args.min_MAF < 0.0) or (args.min_MAF > 1.0):
        raise ValueError('--min_MAF  must be float > 0.0 and < 1.0')

    if args.LDScore_window_autosome is None:
        raise ValueError('--LDScore_window_autosome is required.')
    if args.LDScore_window_autosome < 0.0:
        raise ValueError('--LDScore_window_autosome must be float > 0.0')

    if args.prune_window is None:
        raise ValueError('--prune_window is required.')
    if args.prune_window <= 0:
        raise TypeError('--prune_window must be an integer > 0')

    if args.prune_slide is None:
        raise ValueError('--prune_slide is required.')
    if args.prune_slide <= 0:
        raise ValueError('--prune_slide must be an int > 0')

    if args.prune_cutoff is None:
        raise ValueError('--prune_cutoff is required')
    if (args.prune_cutoff < 0.0) or (args.prune_cutoff > 1.0):
        raise TypeError('--prune_cutoff must be a float > 0.0 and < 1.0')

    if args.clump_cutoff is None:
        raise ValueError('--clump_cutoff is required.')
    if (args.clump_cutoff < 0.0) or (args.clump_cutoff > 1.0):
        raise ValueError('--clump_cutoff must be a float between 0.0 and 1.0')

    if args.max_distance_clump is None:
        raise ValueError('--max_distance_clump is required.')
    if args.max_distance_clump <= 0:
        raise TypeError('--max_distance_clump must be an int > 0')

    int_dir = args.int_dir

    if not os.path.exists(int_dir):
        os.makedirs(int_dir)

    if int_dir[-1] != '/':
        int_dir += '/'

    out_dir = args.out_dir

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if out_dir[-1] != '/':
        out_dir += '/'


    log = Logger(os.path.join(args.out_dir, args.chrom + '.log'))

    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        header = MASTHEAD
        header += 'build_fingerprint_maps.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in opts.keys()]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()
        create_maps(args, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=\
    str(datetime.timedelta(seconds=time_elapsed))))

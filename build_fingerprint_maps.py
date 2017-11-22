"""
Created on Tue Nov 21 15:12:16 2017

@author: Nauman Javed
"""
from __future__ import division
import map_building_functions as build
import subprocess
import os
import itertools as it
import numpy as np
import sys
import argparse
import traceback
import time

__version__ = '1.0.0'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* Fingerprint_Maps\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "*********************************************************************\n"

class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module
    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print >>self.log_fh, msg
        print msg

def create_map(args, log):

    full_path = os.getcwd() + "/"
    intermediate_directory = full_path + "intermediates/"
#    recomb_file = "1000GP_Phase3/genetic_map_chr" + chrom + "_combined_b37.txt*"

    build.extract_similar_SNPs(args.chromosome, VCF_file,
                                intermediate_directory, args.similarity)

    build.create_VCFs(args.chromosome, args.VCF_file,
                        intermediate_directory, args.min_MAF)

    build.sort_VCF(args.chromosome, intermediate_directory)

    build.create_PLINK_binary(args.chromosome, intermediate_directory,
                              args.recomb_directory)

    build.LD_score(args.chromosome, intermediate_directory, args.LD_script,
                    args.LDScore_window_autosome, args.LDScore_window_X)

    build.prune(args.chromosome, intermediate_directory,
                args.prune_window, args.prune_slide, args.prune_cutoff)

    build.LD_separate(chrom, intermediate_directory)

    build.clump(args.chromosome, intermediate_directory,
                args.clump_cutoff, args.max_distance_clump)

    build.reformat_clumps(args.chromosome, intermediate_directory)

    build.detect_negative_LD(args.chromosome, intermediate_directory)

    build.switch_alleles(args.chromosome, full_path, intermediate_directory)

parser = argparse.ArgumentParser()
# Directory specifications'
parser.add_argument('--recomb_directory', default=None, type=str,
    help='Directory pointing to where shapeit recombination files are located.')
parser.add_argument('--chromosome', default=None, type=str,
    help='Chromosome for which to calculate.')
parser.add_argument('--VCF_file', default=None, type=str,
    help='VCF file name for the selected chromosome')
parser.add_argument('--LD_script', default=None, type=str,
    help='Directory containing ldsc.py ')

# SNP filtering and map calculation parameters
parser.add_argument('--similarity', default=0.10, type=float,
    help='Maximum difference in population specific MAF allowed')
parser.add_argument('--min_MAF', default=0.05, type=float,
    help='Minimum minor allele fraction of SNPs that will be included in map')
parser.add_argument('--LDScore_window_autosome', default=1.0, type=float,
    help='Window size(in centimorgans) over which to calculate LDScore.')
parser.add_argument('--LDScore_window_X', default=871, type=int,
    help='Window size(kb for X chromosome) over which to calculate LDScore.')
parser.add_argument('--prune_window', default=1000, type=int,
    help='Window size in # SNPs for PLINK prune function')
parser.add_argument('--prune_slide', default=5, type=int,
    help='Number of SNPs to slide over on each iteration of PLINK prune')
parser.add_argument('--prune_cutoff', default=0.1, type=float,
    help='Maximum r^2 correlation allowed between pruned SNPs')
parser.add_argument('--clump_cutoff', default=0.9, type=float,
    help='Minimum r^2 correlation required for SNPs to be clumped together')
parser.add_argument('--max_distance_clump', default=1500, type=float,
    help='Maximum distance in kb a SNP can be from index SNP when forming clump')

if __name__ == "__main__":
    args = parser.parse_args()
    accepted_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8',
                            '9', '10', '11', '12', '13', '14', '15', '16',
                            '17', '18', '19', '20', '21', '22', 'X']

    if not isinstance(args.recomb_directory, str):
    	raise TypeError('--recomb_directory must be string pointing to where recomb files are stored')
    if args.recomb_directory is None:
        raise ValueError('--recomb_directory is required.')

    if not isinstance(args.chromosome, str):
        raise TypeError('--chromosome must be a string 1-22 or X')
    if args.chromosome is None or args.chromosome not in accepted_chromosomes:
        raise ValueError('--chromosome is required.')

    if not isinstance(args.VCF_file, str):
            raise TypeError('--VCF_file must be full path pointing to VCF file')
    if args.VCF_file is None:
        raise ValueError('--VCF_file_name is required.')

    if not isinstance(args.LD_script, str):
        raise TypeError('--LD_script must be full path pointing to ldsc.py')
    if args.LD_script is None:
        raise ValueError('--LD_script is required.')

    if not isinstance(args.similarity, float):
        raise TypeError('--similarity must be a float between 0.0 and 1.0')
    if args.similarity is None or args.similarity < 0.0 or args.similarity > 1.0:
        raise ValueError('--similarity is required - must be float between 0.0 and 1.0')

    if not isinstance(args.min_MAF, float):
        raise TypeError('--min_MAF must be a float between 0.0 and 1.0')
    if args.min_MAF is None or args.min_MAF < 0.0 or args.min_MAF> 1.0:
        raise ValueError('--min_MAF is required - must be float between 0.0 and 1.0')

    if not isinstance(args.LDScore_window_autosome, float):
        raise TypeError('--LDScore_window_autosome must be a float > 0.0')
    if args.LDScore_window_autosome is None or args.LDScore_window_autosome < 0.0:
        raise ValueError('--LDScore_window_autosome is required - must be float > 0.0')

    if not isinstance(args.LDScore_window_X, int):
        raise TypeError('--LDScore_window_X must be an int > 0')
    if args.LDScore_window_X is None or args.LDScore_window_X <=0:
        raise ValueError('--LDScore_window_X is required - must be int > 0')

    if not isinstance(args.prune_window, int):
        raise TypeError('--prune_window must be an int > 0')
    if args.prune_window is None or args.prune_window <=0:
        raise ValueError('--prune_window is required.')

    if not isinstance(args.prune_slide, int):
        raise TypeError('--prune_slide must be an int > 0')
    if args.prune_slide is None or args.prune_slide <=0:
        raise ValueError('--prune_slide is required.')

    if not isinstance(args.prune_cutoff, float):
        raise TypeError('--prune_cutoff must be a float between 0.0 and 1.0')
    if args.prune_cutoff is None or args.prune_cutoff < 0.0 or args.prune_cutoff > 1.0:
        raise ValueError('--prune_cutoff is required - float between 0.0 and 1.0')

    if not isinstance(args.clump_cutoff, float):
        raise TypeError('--clump_cutoff must be a float between 0.0 and 1.0')
    if args.clump_cutoff is None or args.clump_cutoff < 0.0 or args.clump_cutoff > 1.0:
        raise ValueError('--clump_cutoff is required - float between 0.0 and 1.0')

    if not isinstance(args.max_distance_clump, int):
        raise TypeError('--max_distance_clump must be an int > 0')
    if args.max_distance_clump is None or args.max_distance_clump <=0:
        raise ValueError('--max_distance_clump is required - must be int > 0')


    log = Logger(args.chromsome+'.log')

    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += './build_fingerprint_maps.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
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
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

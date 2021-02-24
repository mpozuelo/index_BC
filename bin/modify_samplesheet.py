#!/usr/bin/env python


import pandas as pd
import argparse
import sys

def parse_args(args=None):
    Description = 'Check samplesheet and add files to bed file'
    Epilog = """Example usage: python modify_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet.")
    parser.add_argument('FILE_OUT', help="Output samplesheet.")
    return parser.parse_args(args)


def add_bed_file(FileIn,FileOut):
    #Open input file
    fi = open(FileIn, 'r')

    # Load mosdepth thresholds.bed.gz into a pandas dataframe
    input = pd.read_csv(fi, delimiter=',', index_col=False, low_memory=False)

    # Open output file
    fo = open(FileOut, 'w')

    basefolder = '/datos/ngs/dato-activo/data/02_rfastq/'
    # Dictionary for bed files
        # Write header
    #fo.write("%s\n" %('\t'.join(l_th[1:])))

    # Compute percentages
    input['fastq1'] = basefolder + input['platform'] + '/' + input['run'] + '/' + input['lane'] + '/' + input['run'] + '_' + input['lane'] + '_read_1.fq.gz'
    input['fastq2'] = basefolder + input['platform'] + '/' + input['run'] + '/' + input['lane'] + '/' + input['run'] + '_' + input['lane'] + '_read_2.fq.gz'


    input.to_csv(fo, index = False)
    fi.close()


def main(args=None):
    args = parse_args(args)
    add_bed_file(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())

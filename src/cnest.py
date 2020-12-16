#!/usr/bin/env python3.8

import argparse
import os
import subprocess
import logging
import sys
import shutil

logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s | %(levelname)s: %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S')
logger = logging.getLogger('CNest')

def get_args():
    parser = argparse.ArgumentParser(description='CNest wrapper')
    subparsers = parser.add_subparsers(help='Steps of CNest', dest='step')
    # step 1
    parser_1 = subparsers.add_parser('step1', help='Create project directory')
    parser_1.add_argument('--project', dest='project', required=True, type=str, help='project name')
    parser_1.add_argument('--bed', dest='bed_file', required=True, type=str, help='input sorted BED3 file name')
    parser_1.add_argument('--debug', dest='debug', default=False, action='store_true', help='Show debug messages')
    # step 2
    parser_2 = subparsers.add_parser('step2', help='BAM/CRAM to binary')
    parser_2.add_argument('--project', dest='project', required=True, type=str, help='project name')
    parser_2.add_argument('--sample', dest='sample_id', required=True, type=str, help='sample name')
    parser_2.add_argument('--input', dest='input_file', required=True, type=str, help='input BAM/CRAM file name')
    parser_2.add_argument('--debug', dest='debug', default=False, action='store_true', help='Show debug messages')
    # step 3
    parser_3 = subparsers.add_parser('step3', help='Gender QC')
    parser_3.add_argument('--project', dest='project', required=True, type=str, help='project name')
    parser_3.add_argument('--debug', dest='debug', default=False, action='store_true', help='Show debug messages')
    args = parser.parse_args()
    return args


def step1(project, bed_path):
    accepted_chr = [str(i) for i in range(1,23)] + ['X', 'Y']
    # make project (sub-)directories if not existed
    dirs = [f'{project}', f'{project}/txt',
    f'{project}/bin', f'{project}/tmp']
    for mydir in dirs:
        if not os.path.exists(mydir):
            os.mkdir(mydir)
    index_path = f'{project}/index.txt'
    index_tab_path = f'{project}/index_tab.txt'
    with open(bed_path) as fin, open(index_path, 'w') as fix, open(index_tab_path, 'w') as ftab:
        for line in fin:
            elements = line.split('\t')
            # index.txt uses the same chrom names as BED (BAM/CRAM)
            # index_tab.txt uses 1-22, 23 for X and 24 for Y
            chrom = elements[0].replace('chr', '')
            assert chrom in accepted_chr, f'Only autosomes and X/Y are allowed. Your chromosome is {chrom}'
            if chrom == 'X':
                chrom = '23'
            elif chrom == 'Y':
                chrom = '24'
            fix.write(f'{elements[0]}:{elements[1]}-{elements[2]}')
            ftab.write(f'{chrom}\t{elements[1]}\t{elements[2]}')


def step2(project_root, sample_id, input_file):
    """
    	project_root=/output_location/${project_name}
		echo $project_name $input_file_name $sample_id
		/software/applications/ngscnv/ngs ${project_root} bam-to-rd /input_location/$input_file_name $sample_id
		/software/applications/ngscnv/ngs ${project_root} rd-dump $sample_id > ${project_root}/txt/$sample_id
		Rscript /resources/run.R processtobin ${project_root} $sample_id
		rm ${project_root}/$sample_id
		rm ${project_root}/txt/$sample_id

        TODO:
        1. Add BAM/CRAM index file check
        2. 
    """
    logger.info('Start step2')
    cmd1 = ['/software/applications/ngscnv/ngs', project_root, 'bam-to-rd', input_file, sample_id]
    cmd2 = ['/software/applications/ngscnv/ngs', project_root, 'rd-dump', sample_id]
    cmd3 = ['Rscript', '/resources/run.R', 'processtobin', project_root, sample_id]
    logger.debug('CMD=' + " ".join(cmd1))
    process1 = subprocess.run(cmd1, capture_output=True)
    if process1.returncode == 0:
        logger.info('bam-to-rd done.')
        logger.debug('CMD=' + " ".join(cmd2))
        outpath = os.path.join(project_root, 'txt', sample_id)
        logger.debug('OUT=' + outpath)
        with open(outpath, 'w') as f:
            process2 = subprocess.run(cmd2, stdout=f)
    else:
        logger.error(f'bam-to-rd failed with exit code {process1.returncode}.\n{process1.stderr}')
    if process2.returncode == 0:
        logger.info('rd-dump done.')
        logger.debug('CMD=' + " ".join(cmd3))
        process3 = subprocess.run(cmd3, capture_output=True)
    if process3.returncode == 0:
        logger.info('processtobin done.')
    # clean temp files
    os.remove(os.path.join(project_root, sample_id))
    os.remove(os.path.join(project_root, 'txt', sample_id))
    logger.info('Step 2 done!')


def step3(project):
    """
    	project_root=/output_location/${project_name}
		index_file=${project_root}/index_tab.txt
		qc_file=${project_root}/gender_qc.txt
		gender_file=${project_root}/gender_classification.txt
		Rscript /resources/run.R classify_gender ${project_root} ${index_file} ${qc_file} ${gender_file}
    """
    project_root = '/output_location/' + project
    index_file = project_root + '/index_tab.txt'
    qc_file = project_root + '/gender_qc.txt'
    gender_file = project_root + '/gender_classification.txt'
    cmd = ['Rscript', '/resources/run.R', 'classify_gender', project_root, index_file, qc_file, gender_file]
    process = subprocess.run(cmd, capture_output=True)
    if process.returncode == 0:
        logger.info('classify_gender done.')


if __name__ == '__main__':
    args = get_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    if args.step == 'step1':
        # step 1
        step1(args.project, args.bed_file)
    elif args.step == 'step2':
        # step 2
        step2(args.project, args.sample_id, args.input_file)
    elif args.step == 'step3':
        # step 3
        step3(args.project)
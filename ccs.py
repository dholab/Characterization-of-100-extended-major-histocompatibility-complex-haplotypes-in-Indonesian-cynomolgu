import subprocess
import os
import utils
import status
import logging
import sys

# initialize logger
log = logging.getLogger(__name__)

def makeCcs(subreads, out_dir, minPredictedAccuracy='0.9', minLength='1000', maxLength='1500'):
    '''use smrtlink ccs to produce consensus sequence'''

    # path to smrtlink ccs
    smrtlink_ccs_path = '/slipstream/SMRT4/SMRT/smrtcmds/bin/ccs'

    # check that subreads file exists
    if os.path.exists(subreads) == False:
        status.printStatus('Error: Specified subread file does not exist. Check your file path and try again.')
        return

    # filename of input file
    subreads_basename = os.path.splitext(os.path.basename(subreads))[0]
    print(subreads_basename)

    # create output directory if it doesn't exist
    utils.createOutputFolder(out_dir)

    # call ccs
    cmd = [smrtlink_ccs_path,
           '--minPredictedAccuracy',
           minPredictedAccuracy,
           '--minLength',
           minLength,
           '--maxLength',
           maxLength,
           subreads,
           out_dir + '/' + subreads_basename + '.ccs.bam']

    status.printStatus('CCS command: ' + ' '.join(cmd))
    status.printStatus('CCS processing of ' + subreads + ' started')
    subprocess.call(cmd)
    status.printStatus('CCS processing of ' + subreads + ' completed')
    status.printStatus('Output CCS file saved to ' + out_dir + '/' + subreads_basename + '.ccs.bam')

    # create fastq file
    fastq_path = makeFastq(out_dir + '/' + subreads_basename + '.ccs.bam')
    return fastq_path

def makeFastq(ccs_bam):
    '''use smrtlink bam2fq to produce gzip compressed FASTQ file from CCS bam'''

    # path to smrtlink bam2fastq
    smrtlink_bam2fastq_path = '/slipstream/SMRT4/SMRT/smrtcmds/bin/bam2fastq'

    # create fastq output file name
    ccs_basename = os.path.splitext(os.path.basename(ccs_bam))[0]
    fastq_output = os.path.dirname(ccs_bam) + '/' + ccs_basename

    # call bam2fastq
    cmd = [smrtlink_bam2fastq_path,
           ccs_bam,
           '-o',
           fastq_output,
           ]

    status.printStatus('bam2fastq command: ' + ' '.join(cmd))
    status.printStatus('bam2fastq processing of ' + ccs_bam + ' started')
    subprocess.call(cmd)
    status.printStatus('bam2fastq processing of ' + ccs_bam + ' completed')
    status.printStatus('gzip compressed FASTQ file saved to ' + fastq_output + '.fastq.gz')

    # return path to output fastq file
    return fastq_output + '.fastq.gz'

# if __name__ == '__main__':  # if run directly from the command line
#     # command line parameters
#     import argparse
#     parser = argparse.ArgumentParser()
#     parser.add_argument("out_dir", help='Folder that will store all output files')
#     parser.add_argument("--subreads", required=True,
#                         help='Path to file of PacBio subreads. Will be converted to CCS file.')
#     parser.add_argument("--ccsMinAccuracy", required=False,
#                         help='Set minPredicted accuracy (from 0-1) for retaining CCS reads. Default=0.9. Recommend 0.999 for de novo allele discovery.')
#     parser.add_argument("--ccsMinLength", required=False,
#                         help='Set minLength in bp for retaining CCS reads. Default=1000. Set to minimum expected amplicon size.')
#     parser.add_argument("--ccsMaxLength", required=False,
#                         help='Set maxLength in bp for retaining CCS reads. Default=1500. Set to minimum expected amplicon size.')
#     args = parser.parse_args()
#
#     # make output folder if it doesn't exist
#     utils.createOutputFolder(args.out_dir)
#
#     # configure log to stdout
#     logging.basicConfig(filename=args.out_dir + '/log.txt', filemode='w', level=logging.DEBUG,
#                         format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
#
#     # run with command line parameters
#     d = {}
#     d['subreads'] = args.subreads
#     d['out_dir'] = args.out_dir
#     if args.ccsMinAccuracy is not None: d['minPredictedAccuracy'] = args.ccsMinAccuracy
#     if args.ccsMinLength is not None: d['minLength'] = args.ccsMinLength
#     if args.ccsMaxLength is not None: d['maxLength'] = args.ccsMaxLength
#
#     # log command line
#     status.printStatus('Command line statment: ' + ' '.join(sys.argv))
#
#     # run makeCcs function
#     makeCcs(**d)

# # test invokation
#     d = {}
#     d['subreads'] = '/slipstream/pacbio/pacbio_raw/pacbio48/3_C01/m54178_170519_124037.subreads.bam'
#     d['out_dir'] = '/slipstream/shared_data/19070/pacbio48-default-minQuality/11'
#     d['minPredictedAccuracy'] = '0.9'
#     d['minLength'] = '1000'
#     d['maxLength'] = '1500'
#
#     # log command line
#     status.printStatus('CCS paramaeters: ' + d)
#
#     # run makeCcs function
#     makeCcs(**d)

# test bam2fastq
# makeFastq('/slipstream/shared_data/19070/pacbio48//20170604082124/ccs//m54178_170519_124037.subreads.ccs.bam')
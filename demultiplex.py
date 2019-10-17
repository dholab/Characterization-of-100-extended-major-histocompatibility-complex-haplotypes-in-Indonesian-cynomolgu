import labkeyInteract
import re
from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
import gzip
import logging
import utils
import status
import shutil
import time

# specify a Pacbio run identifier (e.g., PacBio48)
# get run_id corresponding to this identifier from 'runs' list
# then select records from this run from 'Samples list'
# use labkeyInteract class to communicate with Labkey

# initialize logger
log = logging.getLogger(__name__)

def getRunId(pacbio_id):
    '''inherit pacbio_id (e.g., PacBio48) from parent function and retrieve run identifier'''

    # get Pacbio run ID for specified run
    # necessary because Samples table stores run_id as foreign key lookup to runs table
    # debug modification by JRC 09202018

    pacbio_run_id = labkeyInteract.LabkeySelectRows()
    pacbio_run_id.serverContext('dho/pacbio')
    pacbio_run_id.set_filters('pacbio_id', pacbio_id)
    result = pacbio_run_id.selectRows(labkey_schema='lists', labkey_table='runs')

    # debug modification by JRC 09202018
    print('result is')
    print(result)
    time.sleep(5)
    # extract run number from result
    runNumber = result['rows'][0]['run_num']

    # log whether pacbio_id corresponding to run_id is found
    if runNumber != '':
        status.printStatus(pacbio_id + ' found in dholk.primate.wisc.edu')

    return runNumber

def getSamples(pacbio_id):
    '''retrieve sample information from genotyping Samples table'''

    # get runId corresponding to pacbio_id
    runId = getRunId(pacbio_id)

    # get samples from specified PacBio run
    pacbio_samples = labkeyInteract.LabkeySelectRows()
    pacbio_samples.serverContext('dho/pacbio')
    pacbio_samples.set_filters('run_id', runId)
    result = pacbio_samples.selectRows(labkey_schema='genotyping', labkey_table='Samples')

    # log count of samples in pacbio_id
    status.printStatus(str(result['rowCount']) + ' samples detected in ' + pacbio_id)
    status.printStatus('Barcode configuration')
    # log information on each sample
    print('OC_ID\tForward Barcode\tReverse Barcode')

    samples = {} # initialize samples dictionary

    for i in result['rows']:
        # use oc_id if it exists, otherwise use animal_id to identify sample name
        if i['oc_animal_id'] == None:
            sample_name = i['animal_id']
        else:
            sample_name = i['oc_animal_id']

        # run normalizeBarcodes to create PacBio standard identifiers
        renamed_barcodes = normalizeBarcodes(i)

        # print samples
        print(sample_name + '\t' + renamed_barcodes[0]+ '\t' + renamed_barcodes[1])

        # create dictionary with sample name and barcodes
        samples[sample_name] = [renamed_barcodes[0], renamed_barcodes[1]]

    return samples

def normalizeBarcodes(row):
    '''standardize barcode names so they match pacbio_barcodes_384.fasta barcode names'''

    # find forward and reverse barcode identifiers
    fwd = re.search('[0-9]*', str(row['forward_barcode']))
    rev = re.search('[0-9]*', str(row['reverse_barcode']))

    # rename barcodes to use PacBio standard identifiers
    renamed_fwd_barcode = fwd.group().zfill(4) + '_Forward'
    renamed_rev_barcode = rev.group().zfill(4) + '_Forward'

    return [renamed_fwd_barcode, renamed_rev_barcode]

def pacbioBarcodeDict():
    '''make dictionary from PacBio barcodes'''

    # make dictionary
    pacbioBarcodes = {}

    print('parsing barcodes')
    for record in SeqIO.parse('/slipstream/oc/pacbio/pacbioSequel/ref/pacbio_barcodes_384.fasta', "fasta"):
        pacbioBarcodes[record.id] = str(record.seq)

    return pacbioBarcodes

def makeBarcodeManifest(samples, out_fasta):
    '''create FASTA file containing serialized barcodes for use with PacBio SMRT LINK'''

    # bug hunting 09202018 JRC
    # print('makeBarcodeMAnifest variables are ')
    # print(samples)
    # print(out_fasta)

    # create PacBio barcode dictionary to lookup against
    pacbioLookup = pacbioBarcodeDict()

    # loop over all samples

    # counter for uniquifying barcodes
    ct = 1

    # write output to file
    with open(out_fasta, 'w') as the_file:
        for seq_name, barcode_seqs in samples.items():
            # print FASTA header
            the_file.write('>' + str(ct).rjust(3, '0') + '-' + seq_name + '_' + barcode_seqs[0] + '\n')

            # get barcode name
            forward_sample_barcode = (barcode_seqs[0])

            # print barcode sequence
            the_file.write(pacbioLookup[forward_sample_barcode] + '\n')

            # increment serial
            ct += 1

            # print FASTA header
            the_file.write('>' + str(ct).rjust(3, '0') + '-' + seq_name + '_' + barcode_seqs[1] + '\n')

            # get barcode name
            reverse_sample_barcode = (barcode_seqs[1])

            # print barcode sequence
            the_file.write(pacbioLookup[reverse_sample_barcode] + '\n')

            # increment serial
            ct += 1

    # log creation of FASTA file
    print('--[' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '] ' + 'FASTA file written to ' + out_fasta + ' --')

def parseBarcodes(samples, input_ccs_fastq, out_dir):
    '''parse barcodes from gzip-compressed FASTQ of PacBio CCS reads'''

    # create output directory if it doesn't exist
    utils.createOutputFolder(out_dir)

    # create PacBio barcode dictionary to lookup against
    pacbioLookup = pacbioBarcodeDict()

    # create dictionary of sample IDs and barcode sequences
    searchDict = {}
    for seq_name, barcode_seqs in samples.items():
        searchDict[seq_name] = [pacbioLookup[barcode_seqs[0]], pacbioLookup[barcode_seqs[1]]]

    # open gzip-compressed FASTQ
    with gzip.open(input_ccs_fastq, "rt") as handle:

        # make dictionary to hold barcode-split seq records
        perBarcodeDict = {}

        # initialize dictionary with names of each sample
        for j in searchDict:
            perBarcodeDict[j]=[]

        # log every 1000 sequences processed
        log_every_n = 1000

        # iterate through generator containing FASTQ sequences
        for idx, i in enumerate(SeqIO.parse(handle, "fastq")):

            # print status message every 1000 sequences processed
            if (idx % log_every_n) == 0:
                status.printStatus(str(idx) + ' FASTQ reads demultiplexed')

            # for each sequence, look for the presence of barcodes at the start and end
            for j in searchDict:
                # redo to use re.search to find barcodes not at very end of sequence
                # if i.seq.startswith(searchDict[j][0]) and i.seq.endswith(searchDict[j][1]):

                # regular expression to find barcodes in forward orientation
                prog = re.compile(searchDict[j][0] + ('.*') + searchDict[j][1])

                # test if regular expression is found in sequence
                # need to cast i.seq to string to use re.search

                if prog.search(str(i.seq)):
                    # write matching barcodes to perBarcodeDict - store in memory
                    x = perBarcodeDict[j]
                    x.append(i)
                    perBarcodeDict[j]= x

                # handle inserts in the opposite orientation
                # create Biopython sequence object containing barcode sequences
                forward_seq = Seq(searchDict[j][0])
                reverse_seq = Seq(searchDict[j][1])

                # reverse complement
                forward_seq_rc = forward_seq.reverse_complement()
                reverse_seq_rc = reverse_seq.reverse_complement()

                # find FASTQ sequences matching reverse complemented barcodes
                # if i.seq.startswith(forward_seq_rc) and i.seq.endswith(reverse_seq_rc):

                # because of the SMRTBell orientation, second barcode gets listed first in reverse complement orientation
                prog = re.compile(str(reverse_seq_rc) + '.*' + str(forward_seq_rc))

                # need to cast i.seq to string to use re.search
                if prog.search(str(i.seq)):
                    # store matches in dictionary
                    x = perBarcodeDict[j]
                    x.append(i)
                    perBarcodeDict[j]= x

        # write output files containing reads matching each barcode
        for i in perBarcodeDict:
            count = SeqIO.write(perBarcodeDict[i], out_dir + '/' + i + '.fastq', 'fastq')

            # compress fastq file and remove uncompressed version

            with open(out_dir + '/' + i + '.fastq', 'rb') as f_in:
                with gzip.open(out_dir + '/' + i + '.fastq.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.remove(out_dir + '/' + i + '.fastq') # remove uncompressed

            # log
            status.printStatus(str(count) + ' barcoded reads saved from sample ' + i )
            status.printStatus('gzip-compressed demultipled FASTQ file saved to ' + out_dir + '/' + i + '.fastq.gz')

def makeSmrtlinkFasta(pacbioRunId, makeSmrtlinkBarcodes):
    '''get sample information from DHOLK then make barcode manifest'''

    # retrieve sample information from dholk
    samples = getSamples(pacbioRunId)

    # make smrtlink compatible FASTA barcode file
    makeBarcodeManifest(samples, makeSmrtlinkBarcodes)

def demultiplexFastq(ccs, pacbioRunId, demultiplexedFastqPath):
    '''Using barcode information from dholk, use python to demultiplex barcoded reads'''

    # print('PacBio Run ID is ' + str(pacbioRunId))
    # retrieve sample information from dholk
    samples = getSamples(pacbioRunId)
    # print(samples)

    # print(ccs)

    # split FASTQ by barcodes
    parseBarcodes(samples, ccs, demultiplexedFastqPath)

# test demultiplexing
# demultiplexFastq('/slipstream/shared_data/19070/pacbio48/20170604130117/ccs/m54178_170519_124037.subreads.ccs.fastq.gz', 'PacBio48', '/slipstream/shared_data/19364/09/')

# if __name__ == '__main__':  # if run directly from the command line
#     # command line parameters
#     import argparse
#     parser = argparse.ArgumentParser()
#     parser.add_argument("out_dir", help='Folder that will store all output files')
#     parser.add_argument("--ccs", required=False,
#                         help='Path to file of PacBio CCS reads')
#     parser.add_argument("--pacbioRunId", required=False,
#                         help='(e.g., PacBio48 - from dholk PacBio sample database')
#     parser.add_argument("--makeSmrtlinkBarcodes", required=False,
#                         help='Name of smrtlink-compatible barcode FASTA file')
#     parser.add_argument("--demultiplexedFastqPath", required=False,
#                         help='Path to FASTQ files split by barcodes specified in dholk')
#     args = parser.parse_args()
#
#     # make output folder if it doesn't exist
#     utils.createOutputFolder(args.out_dir)
#
#     # configure log to stdout
#     logging.basicConfig(filename=args.out_dir + '/log.txt', filemode='w', level=logging.DEBUG,
#                         format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
#
#     # log command line
#     status.printStatus('Command line statment: ' + ' '.join(sys.argv))
#
#     # make smrtlink-compatible barcode file if makeSmrtLinkBarcodes specified
#     if args.makeSmrtlinkBarcodes is not None:
#         makeSmrtlinkFasta(args.pacbioRunId, args.makeSmrtlinkBarcodes)
#
#     # make demultiplexed FASTQ if demultiplexedFastqPath specified
#     if args.demultiplexedFastqPath is not None:
#         demultiplexFastq(args.ccs, args.pacbioRunId, args.demultiplexedFastqPath)

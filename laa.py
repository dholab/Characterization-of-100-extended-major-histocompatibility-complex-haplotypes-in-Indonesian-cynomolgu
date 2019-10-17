import subprocess
import os
import sys
import status
import utils
import fileinput
from Bio import SeqIO

# these are functions for running long-amplicon analysis using SMRT Link v5.1 on dhogal
# expects demultiplexed FASTQ input files

def extractSequenceNames(gzip_fastq):
    '''convert FASTQ to FASTA and then extract sequence names to new file'''

    # path to reformat.sh, update as needed
    bbmap_reformat_sh = '/slipstream/oc/jrcanalysisdata/mhcAnalysis/bbmap/reformat.sh'

    # create temporary whitelist sequence path
    whitelist_sequences = gzip_fastq + '.whitelist.tmp.txt'
    print(whitelist_sequences)

    # create reformat.sh command to convert fastq to fasta
    cmd = [bbmap_reformat_sh ,
           'in=' + gzip_fastq,
           'out=' + whitelist_sequences + '.tmp.fasta']

    # print bbmap command
    status.printStatus(' '.join(cmd))

    # call reformat.sh
    subprocess.call(cmd)

    # need to remove trailing /ccs from FASTA file
    # use code from https://stackoverflow.com/questions/17140886/how-to-search-and-replace-text-in-a-file-using-python

    with fileinput.FileInput(whitelist_sequences + '.tmp.fasta', inplace=True) as file:
        for line in file:
            print(line.replace('/ccs', ''), end='')

    # extract sequence names to new file
    with open(whitelist_sequences, 'w') as the_file:
        for seq_record in SeqIO.parse(whitelist_sequences + '.tmp.fasta', "fasta"):
            the_file.write(seq_record.id + '\n')

    # return path to fasta_output
    return whitelist_sequences

def runLongAmpliconAnalysis(subreadsetXML, whitelistSequences, outputPrefix, minLength='1000', maxLength='1500', maxReads='20000', maxClusteringReads='5000'):
    '''run SMRT Link v5 long amplicon analysis'''

    # runs LAA to generate amplicon sequences from PacBio Sequel data
    # subreadsetXML can be from a single dataset, or merged datasets where new XML files are created using dataset create
    # whitelistFasta is a file containing sequences that will be analyzed by LAA, typically sequences from a single sample
    # defaults are set for typical MHC class I genotyping and should be adjusted depending on target
    # note: LAA default minLength=3000 will cause most of our analyses to fail so minLength should almost always be set
    # increasing maxClusteringReads will allow more alleles to be detected at the expense of speed:
    # LAA default of 500 clustering reads runs each sample in ~2 minutes, MHC class I default of 10000 takes ~30 minutes
    # but detects more alleles. Setting even higher values like 100,000 clustering reads causes runtimes of several hours.
    # maxReads can be set very high to ensure that all reads are used to accurately define clusters. This doesn't significantly
    # impact runtime.

    # use outputPrefix to specify the folder and prefix for output files
    # eg '/slipstream/shared_data/19364/09/'
    # eg '/slipstream/shared_data/19364/09/BM115.

    # path to SMRT Link v6.0 LAA
    laa_path = '/slipstream/oc/pacbio/smrtlink_v6/smrtcmds/bin/laa'

    # create output folder if it doesn't exist
    utils.createOutputFolder(os.path.dirname(outputPrefix))

    # create laa command
    laa_cmd = [laa_path,'--whitelist=' + whitelistSequences,'--logFile=' + outputPrefix + '.log.txt','--resultFile=' + outputPrefix + '.amplicon_analysis.fastq','--junkFile=' + outputPrefix + '.amplicon_analysis_chimeras_noise.fastq','--reportFile=' + outputPrefix + '.amplicon_analysis_summary.csv','--inputReportFile=' + outputPrefix + '.amplicon_analysis_input.csv','--subreadsReportPrefix=' + outputPrefix + '.amplicon_analysis_subreads',subreadsetXML]
    print(laa_cmd)
#               '--minLength=' + minLength,
#               '--maxLength=' + maxLength,
#               '--maxReads=' + maxReads,
#               '--maxClusteringReads=' + maxClusteringReads,'--whitelist=' + whitelistSequences,
#               '--logFile=' + outputPrefix + '.log.txt',
#               '--resultFile=' + outputPrefix + '.amplicon_analysis.fastq',
#               '--junkFile=' + outputPrefix + '.amplicon_analysis_chimeras_noise.fastq',
#               '--reportFile=' + outputPrefix + '.amplicon_analysis_summary.csv',
#               '--inputReportFile=' + outputPrefix + '.amplicon_analysis_input.csv',
#               '--subreadsReportPrefix=' + outputPrefix + '.amplicon_analysis_subreads',
#               subreadsetXML]
#    laa_cmd = [laa_path,
#               '--minLength=' + minLength,
#               '--maxLength=' + maxLength,
#               '--maxReads=' + maxReads,
#               '--maxClusteringReads=' + maxClusteringReads,
#               '--whitelist=' + whitelistSequences,
#               '--logFile=' + outputPrefix + '.log.txt',
#               '--resultFile=' + outputPrefix + '.amplicon_analysis.fastq',
#               '--junkFile=' + outputPrefix + '.amplicon_analysis_chimeras_noise.fastq',
#               '--reportFile=' + outputPrefix + '.amplicon_analysis_summary.csv',
#               '--inputReportFile=' + outputPrefix + '.amplicon_analysis_input.csv',
#               '--subreadsReportPrefix=' + outputPrefix + '.amplicon_analysis_subreads',
#               subreadsetXML]


    # print laa command
    status.printStatus(' '.join(laa_cmd))

    # call laa
    subprocess.call(laa_cmd)

    # return path to LAA fastq output
    return outputPrefix + 'amplicon_analysis.fastq'

def cleanupTempFiles(f):
    # remove temporary files
    os.remove(f)

def generateWhitelistAndRunLaa(gzip_fastq, subreadsetXML, output_folder, minLength='1000', maxLength='1500', maxReads='20000', maxClusteringReads='5000'):
    '''workflow to generate whitelist sequences from demultiplexed FASTQ and use as input for LAA'''

    print(gzip_fastq)
    # generate whitelist sequences
    whitelist_path = extractSequenceNames(gzip_fastq)

    # extract sample name from .fastq.gz and use to create full output file location and prefix
    gzip_fastq_basename = os.path.splitext(os.path.basename(gzip_fastq))[0]
    outputPrefix = output_folder + '/' + gzip_fastq_basename

    # run LAA
    runLongAmpliconAnalysis(subreadsetXML,
                            whitelist_path,
                            outputPrefix,
                            minLength,
                            maxLength,
                            maxReads,
                            maxClusteringReads)

    # remove temp files
    cleanupTempFiles(gzip_fastq + '.whitelist.tmp.txt')
    cleanupTempFiles(gzip_fastq + '.whitelist.tmp.txt' + '.tmp.fasta')

### tests ###
# test fastq to fasta
# extractSequenceNames('/slipstream/shared_data/19364/09/BM115.fastq.gz', '/slipstream/shared_data/19364/09/BM115.txt')

# test amplicon analysis function
# runLongAmpliconAnalysis('/slipstream/pacbio/pacbio_raw/UWBC/2017-05-23/3_C01/m54178_170519_124037.subreadset.xml',
#                         '/slipstream/shared_data/19364/09/BM115.txt',
#                        '/slipstream/shared_data/19364/10/BM115/BM115.')


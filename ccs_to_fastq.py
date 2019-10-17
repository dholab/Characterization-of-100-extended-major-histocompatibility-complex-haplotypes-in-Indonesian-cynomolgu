import os
import subprocess
import status

def makeFastq(ccs_bam):
    '''use smrtlink bam2fq to produce gzip compressed FASTQ file from CCS bam'''
    '''/slipstream/oc/pacbio/smrtlink_v6/smrtcmds/bin/'''

    # path to smrtlink bam2fastq
    smrtlink_bam2fastq_path = '/slipstream/oc/pacbio/smrtlink_v6/smrtcmds/bin/bam2fastq'

    # create fastq output file name
    ccs_basename = os.path.splitext(os.path.basename(ccs_bam))[0]
    fastq_output = os.path.dirname(ccs_bam) + '/' + ccs_basename
    print(fastq_output)

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

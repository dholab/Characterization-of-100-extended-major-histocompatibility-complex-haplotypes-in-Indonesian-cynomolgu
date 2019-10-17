import subprocess
import os
import status
import utils
import pandas as pd
import logging
from shutil import copyfile
import sys
import labkeyInteract

# initialize logger
log = logging.getLogger(__name__)

def mapReads(in_fastq, ref_fasta, out_dir, experiment):
    '''use mapPacBio.sh from bbmap to identify reference sequenecs matched by one or more PacBio reads with no substitutions (indels allowed)'''

    # mapPacBio path (first part gets path to folder running script)
    bbmap_pacbio = (os.path.dirname(os.path.realpath(__file__))) + '/bbmap_37_28/mapPacBio.sh'

    # get sample name from input file
    # need to strip off .gz and .fastq extensions sequentially

    sample_name = os.path.splitext(os.path.splitext(os.path.basename(in_fastq))[0])[0]
    print('Sample name: ' + sample_name)

    # create output genotyping folder if it doesn't exist
    sample_dir = utils.createOutputFolder(out_dir + '/genotyping/' + sample_name)

    # create bbmap command
    cmd = [bbmap_pacbio,
           'in=' + in_fastq,
           'ref=' + ref_fasta,
           'covstats=' + sample_dir + '/' + sample_name + '.covstats.tmp.txt',
           'outm=' + sample_dir + '/' + sample_name + '.mapped.bam',
           'outu=' + sample_dir + '/' + sample_name + '.unmapped.fastq.gz',
           'statsfile=' + sample_dir + '/' + sample_name + '.mapping_stats.txt',
           'subfilter=0',
           'nzo=t',
           'ambiguous=all',
           'maxlen=1500',
           'minid=0.9',
           'maxindel=10',
           'minratio=0.8',
           'twocolumn=t',
           'ow=t']

    # print bbmap command
    status.printStatus(' '.join(cmd))

    # call bbmap
    # suppress stats output (saved to file, no need to clutter stderr)
    # FNULL = open(os.devnull, 'w')
    subprocess.call(cmd)
    # FNULL.close()

    # add descriptors to covstats output
    with open(sample_dir + '/' + sample_name + '.covstats.tmp.txt', 'r') as f:
        with open(sample_dir + '/' + sample_name + '.covstats.txt', 'w') as g:
            for idx, line in enumerate(f):
                # print header in first line, otherwise value of sample_name
                if idx == 0:
                    g.write('sample_name' + '\t' + line.rstrip('\n') + '\t' + 'ref_fasta\tanalysis_path\texperiment\n')
                else:
                    g.write(sample_name + '\t' + line.rstrip('\n') + '\t' + ref_fasta + '\t' + out_dir + '\t' + experiment + '\n')

    # remove temporary covstats.tmp.txt file after covstats.txt with sample ID prepared
    if os.path.exists(sample_dir + '/' + sample_name + '.covstats.tmp.txt'):
        os.remove(sample_dir + '/' + sample_name + '.covstats.tmp.txt')

    # copy reference file to output folder
    copyfile(ref_fasta, out_dir + '/genotyping/' + os.path.basename(ref_fasta))

    # return covstats file
    return sample_dir + '/' + sample_name + '.covstats.txt'

def pivotTable(covstats_files, out_dir):
    '''make pivot-table display from bbmap genotyping results'''

    # expect list of covstats.txt files from mapReads

    # make concatentated dataframe
    df = pd.concat((pd.read_csv(f, sep = '\t', dtype={'Avg_fold':'float'}) for f in covstats_files))

    # change Avg_fold data to float and round to one decimal
    df.Avg_fold = df.Avg_fold.round(1)

    # rename columns in dataframe
    df.rename(index=str, columns={'#ID' : 'allele', 'Avg_fold' : 'read_count'}, inplace=True)

    # pivot data
    pivoted = pd.pivot_table(df, values='read_count', index=['allele'], columns='sample_name').reset_index()

    # create output genotyping folder if it doesn't exist
    genotyping_dir = utils.createOutputFolder(out_dir + '/genotyping/')

    # create CSV output files
    df.to_csv(genotyping_dir + '/genotyping_list.tsv', sep='\t')
    pivoted.to_csv(genotyping_dir + '/genotyping_pivot.tsv', sep='\t')

    # import data to labkey
    importLabkey(df)

def mapReadsFolder(fastq_folder, ref_fasta, out_dir, experiment):
    '''map FASTQ reads to reference for all files in folder and make pivottable from results'''

    # create list to store covstats paths
    covstats = []

    # count number of fastq files that will be processed
    fastq_count = 0
    for filename in os.listdir(fastq_folder):
        if filename.endswith(".fastq.gz"):
            fastq_count += 1

    # run mapReads on FASTQ files in specific folder

    for idx, filename in enumerate(os.listdir(fastq_folder)):
        if filename.endswith(".fastq.gz"):
            # run mapReads for each file
            # return covstats file path - add to covstats list
            status.printStatus('Genotyping FASTQ file ' + str(idx + 1) + ' of ' + str(fastq_count))
            covstats.append(mapReads(fastq_folder + '/' + filename, ref_fasta, out_dir, experiment))
            continue
        else:
            continue

    status.printStatus('Make pivot table from: ' + ', '.join(covstats))

    # create pivottable
    pivotTable(covstats, out_dir)

def importLabkey(df):
    '''import tabular genotypes into https://dholk.primate.wisc.edu/list/dho/gs/grid.view?listId=1630'''

    # make list of records from tabular dataframe
    labkey_data = df.to_dict('records')

    # add to labkey
    x = labkeyInteract.LabkeyInsertRows()
    x.serverContext('/dho/gs/')
    x.labkey_schema = 'lists'
    x.labkey_table = 'pacbio_genotypes'
    x.insertRows(x.labkey_schema, x.labkey_table, labkey_data)

    # log
    status.printStatus(str(len(labkey_data)) + ' sample genotypes added to dholk')

if __name__ == '__main__':  # if run directly from the command line
    # command line parameters
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", help='Folder that will store all output files')
    parser.add_argument("fastq_folder", help='Path to folder containing FASTQ files to genotype')
    parser.add_argument("ref_fasta", help='Path to reference FASTA file to map reads against')
    parser.add_argument("experiment", help='Experiment number')
    args = parser.parse_args()

    # make output folder if it doesn't exist
    utils.createOutputFolder(args.out_dir)

    # configure log to stdout
    logging.basicConfig(filename=args.out_dir + '/log.txt', filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S')

    # log command line
    status.printStatus('Command line statment: ' + ' '.join(sys.argv))

    # map reads and summarize results
    mapReadsFolder(args.fastq_folder, args.ref_fasta, args.out_dir, args.experiment)

    # example invokation
    # anaconda3/bin/python /slipstream/shared_data/pycharm/dhogal/19070/genotyping.py /slipstream/shared_data/19070/pacbio48-default-minQuality/16/ /slipstream/shared_data/19070/pacbio48-default-minQuality/12/fastq/ /slipstream/shared_data/19070/pacbio48-default-minQuality/ipd-mhc-20170523.fasta
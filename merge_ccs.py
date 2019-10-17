import os
import subprocess

# purpose: merge CCS datasets when the same library is run on multiple SMRTCells

# overview
# 1. Merge BAM files with samtools merge
# 2. Create pacbio bam index with smrtlink
# The indices are necessary for other smrtlink tools (like bam2fq) to work correctly

def samtoolsMerge(samtools_path, output_merged_bam_path, *input_bamfiles):
    '''input two or more BAM file paths to merge. Merge wtih samtools merge'''

    # path to samtools, merge command, and output_merged_bam_path
    merge_cmd = [samtools_path, 'merge', output_merged_bam_path]

    # iterate over BAM files and add to merge_cmd
    for i in input_bamfiles:
        merge_cmd.append(i)

    # run samtools merge
    subprocess.call(merge_cmd)

    # index merged BAM
    pacbio_index_bam(output_merged_bam_path)

    # return path to merged BAM file
    return output_merged_bam_path

def pacbio_index_bam(merged_bam):
    '''input merged CCS file and make pacbio BAM index'''

    # run pbindex on merged ccs
    # index saved in same directory as merged bam

    # path to pbindex
    pbindex_path = '/slipstream/SMRT4/SMRT/smrtcmds/bin/pbindex'

    subprocess.call([pbindex_path, merged_bam])


## merge and index PacBio 48

samtoolsMerge('/usr/local/bin/samtools', '/slipstream/pacbio/pacbio_processed/19243/PacBio48/PacBio48.merged.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio48/19243-20170624223255/ccs/m54178_170519_124037.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio48/19243-20170624225635/ccs/m54178_170525_213850.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio48/19243-20170625021449/ccs/m54178_170526_070329.subreads.ccs.bam')

## merge and index PacBio 49

samtoolsMerge('/usr/local/bin/samtools', '/slipstream/pacbio/pacbio_processed/19243/PacBio49/PacBio49.merged.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio49/19243-20170625033514/ccs/m54178_170518_161609.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio49/19243-20170625065217/ccs/m54178_170519_022739.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio49/19243-20170625161508/ccs/m54178_170622_142521.subreads.ccs.bam')

## merge and index PacBio 51

samtoolsMerge('/usr/local/bin/samtools', '/slipstream/pacbio/pacbio_processed/19243/PacBio51/PacBio51.merged.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio51/19243-20170626013948/ccs/m54178_170620_213502.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio51/19243-20170626033640/ccs/m54178_170621_074622.subreads.ccs.bam')

## merge and index PacBio 52

samtoolsMerge('/usr/local/bin/samtools', '/slipstream/pacbio/pacbio_processed/19243/PacBio52/PacBio52.merged.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio52/19243-20170626050830/ccs/m54178_170621_175917.subreads.ccs.bam',
              '/slipstream/pacbio/pacbio_processed/19243/PacBio52/19243-20170626054333/ccs/m54178_170622_041226.subreads.ccs.bam')


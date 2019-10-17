import os
import logging

# initialize logger
log = logging.getLogger(__name__)

def createOutputFolder(out_dir):
    '''Test if output folder and create if it doesn't exist'''

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        log.info('Output directory ' + out_dir + ' created')

    return out_dir
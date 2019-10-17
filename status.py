from datetime import datetime
import logging

# import logger
log = logging.getLogger(__name__)

# print status update with timestamp

def printStatus(status):
    print('--[' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '] ' + status + '--')
    log.info(status)

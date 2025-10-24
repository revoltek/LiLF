from dask_jobqueue import SLURMCluster
from dask.distributed import Client 
import time, sys, glob
import numpy as np
from dask import delayed
from dask.distributed import Client, as_completed

import os

from LiLF import lib_ms, lib_util, lib_log

os.makedirs('logs/dask-test', exist_ok=True)
logger_obj = lib_log.Logger('pipeline-dask-test')
logger = lib_log.logger


s = lib_util.Scheduler(log_dir=logger_obj.log_dir, backend='slurm', container_path=os.path.expanduser('~') + "/pill.simg")

#s._cluster.scale(cores=1)

logger.debug(s._client.scheduler_info())
logger.debug(s._client.run(lambda: os.uname()))


MSs = lib_ms.AllMSs(glob.glob('/beegfs/lofar/boxelaar/deepfields/Elais-N1/mss-dask-test/TC*.MS'), s,)
logger.info(f"Found {len(MSs.getListObj())} MSs to process.")

lib_util.check_rm('/beegfs/lofar/boxelaar/deepfields/Elais-N1/mss-dask-test/TC*.MS-test')
MSs.run(f'DP3 msin=$pathMS msout=$pathMS-test steps=[count]', log='$nameMS_count.log', commandType='DP3')







'''
def MSs_check(MS: lib_ms.MS):
    command = f"DP3 msin={MS.pathMS} msout=. steps=[count]"
    logger.info(f"Generated command for {MS.pathMS}: {command}")
    os.system(command)
    return command

mode = 'future'
if mode == 'future':
    commands_ran = []
    for future in as_completed(s._client.map(MSs_check, MSs.getListObj())): 
        command = future.result()
        commands_ran.append(command)

    print(commands_ran)
    print("Done")
    
elif mode == 'submit':
    futures = s._client.map(MSs_check, MSs.getListObj(), key=MSs.getListStr()) 
    for future in as_completed(futures):
        MS_name = future.key
        command = future.result()
        logger.info(f"Finished processing {MS_name} with command: {command}")
        
    print("Done")




s._client.close()
#'''

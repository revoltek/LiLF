from dask_jobqueue import SLURMCluster
from dask.distributed import Client 
import time, sys, glob
import numpy as np
from dask import delayed
from dask.distributed import Client, as_completed

import os

from LiLF import lib_ms, lib_util, lib_log

logger_obj = lib_log.Logger('pipeline-dask-test')
logger = lib_log.logger


s = lib_util.SLURMScheduler(log_dir=logger_obj.log_dir, container_path=os.path.expanduser('~') + "/pill.simg", walltime='02:00:00', max_cpus_per_node=16)



data_dir = 'mss'
MSs = lib_ms.AllMSs(glob.glob(data_dir + '/TC*.MS'), s,)
logger.info(f"Found {len(MSs.getListObj())} MSs to process.")

lib_util.check_rm(data_dir + '/TC*.MS-test')
MSs.run(f'DP3 msin=$pathMS msout=. steps=[count]', log='$nameMS_count.log')







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

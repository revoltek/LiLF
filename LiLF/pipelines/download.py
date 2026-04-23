# Pipeline to download LOFAR measurement sets via authenticated wget.

import glob
import os
import re
from LiLF import lib_scheduler, lib_log, lib_walker

def run(step):
    
    lib_log.Logger(f'pipeline-{step.kind}-{step.name}')
    s = lib_scheduler.Scheduler(dry_run=False)
    w = lib_walker.Walker(f'pipeline-{step.kind}-{step.name}.walker')
    logger = lib_log.logger

    logger.info(f"Started step: {step.kind} - {step.name}")

    output_dir = step['output']
    download_file = step['download_file']
    macaroon_file = step['macaroon_file']

    os.makedirs(output_dir, exist_ok=True)

    if download_file is not None:
        with w.if_todo('download'):
            with open(download_file,'r') as df:
                logger.info(f'Downloading into {output_dir}...')
                downloaded = glob.glob(os.path.join(output_dir, '*MS'))
                downloaded = [os.path.basename(p) for p in downloaded]
                # add renamed files
                if os.path.exists('renamed.txt'):
                    with open('renamed.txt','r') as flog:
                        downloaded += [line.rstrip('\n') for line in flog]
                with open(macaroon_file,'r') as mf:
                    macaroon = mf.readlines()[-2].strip('\n')

                for i, line in enumerate(df):
                    if line == "\n": continue
                    ms = re.findall(r'L[0-9]*.*_SB[0-9]*_uv', line)[0]
                    if ms+'.MS' in downloaded or ms+'.dppp.MS' in downloaded: continue
                    s.add(f'wget -nv --header "Authorization: bearer {macaroon}" --no-check-certificate "{line[:-1]}" -O - | tar -x -C {output_dir}', log=ms+'_download.log', commandType='general')
                    logger.debug('Queue download of: '+line[:-1])
                s.run(check=True, max_proc=4)
    else:
        logger.warning('No download file specified, finishing.')

    w.alldone()
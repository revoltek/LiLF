#!/usr/bin/env python3

import argparse, sys, glob, os
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5, lib_dd_parallel, lib_cat
import lsmtool

logger_obj = lib_log.Logger('pipeline-simul')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-simul.walker')

parset_dir_para = os.path.dirname(os.path.realpath(__file__))+'/../parsets/LOFAR_ddparallel'
parset_dir_ser = os.path.dirname(os.path.realpath(__file__))+'/../parsets/LOFAR_ddserial'

def main(ms, skymodel, h5parm):
    logger.info("Simulation script initialized.")
    logger.info(f"Input ms: {ms}, skymodel: {skymodel}, h5parm: {h5parm}")
    MSs = lib_ms.AllMSs( ms, s, check_flags=True, check_consistency=True)
    MSs.run('taql "UPDATE $pathMS SET DATA=0"', log='$nameMS_taql.log', commandType='general')
    phase_center = MSs.getListObj()[0].getPhaseCentre()
    # load skymodel and fill the model_column
    skymodel_basename = skymodel.split('/')[-1]
    for MS in MSs.getListStr():
        lib_util.check_rm(MS + '/' + skymodel_basename)
        logger.debug('Copy: ' + skymodel + ' -> ' + MS)
        os.system('cp -r ' + skymodel + ' ' + MS)
    sm = lsmtool.load(skymodel, beamMS=ms[0])
    patches = sm.getPatchNames()
    for patch in patches:
        logger.info(f"Processing patch: {patch}")
        MSs.run(f'DP3 {parset_dir_para}/DP3-predict-beam.parset msin=$pathMS pre.sourcedb=$pathMS/{skymodel_basename} pre.sources={patch} \
                msout.datacolumn=MODEL_DATA pre.beammode=array_factor', log='$nameMS_pre.log', commandType='DP3')
        MSs.run(f'DP3 {parset_dir_para}/DP3-cor.parset msin=$pathMS cor.parmdb={h5parm} msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                       cor.correction=phase000 cor.invert=False cor.direction=Dir{patch[-2:]}', log='$nameMS_cor.log', commandType="DP3")     
        MSs.run('taql "UPDATE $pathMS SET DATA=DATA+MODEL_DATA"', log='$nameMS_taql.log', commandType='general')

    logger.info("Set CORRECTED_DATA = DATA")    
    MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA"',
                    log='$nameMS_taql.log', commandType='general')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation script")
    parser.add_argument("--ms", type=str, nargs='+', help="Input strings for ms (space-separated for multiple values)")
    parser.add_argument("--skymodel", type=str, help="Input string for skymodel") # already divided into patches
    parser.add_argument("--h5parm", type=str, help="Input string for h5parm") # with a direction per patch
    args = parser.parse_args()
    
    ms = args.ms
    skymodel = args.skymodel
    h5parm = args.h5parm
    if ms and skymodel and h5parm:
        main(ms, skymodel, h5parm)  
    else:
        print("No input string provided for ms.")
        sys.exit(1)
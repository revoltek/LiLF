#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa

This pipeline chunk
1. Performs the first chunk of flagging.
'''

import glob, logging

import lib_ms, lib_util


def main():
    #pathDirectoryMain = "/disks/strw3/oei/uGMRTCosmosCut-PiLF"
    pathDirectoryMain = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsTarget/P149.7+03.4/MSs"
    
    pathDirectoryLog  = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/logs/"
    nameFileLog       = "pipeline_uGMRT_flag.log"
    pathFileLog       = pathDirectoryLog + nameFileLog
    
    
    lib_util.printLineBold("Starting log at '" + pathFileLog + "'...")
    
    
    logging.basicConfig(filename = pathFileLog, level = logging.DEBUG)
    logging.info("Started 'pipeline_uGMRT_flag.py'!")
    
    for pathMS in glob.glob(pathDirectoryMain + "/*MS"):
        MSObject = lib_ms.Ms(pathMS)
        print (MSObject.find_nchan())
        print (MSObject.find_chanband())


if (__name__ == "__main__"):
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (C) 2018 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
addcol2ms.py
Create a new column in a measumrent set
"""

import optparse, logging
import casacore.tables as pt
import numpy
import logging

logging.basicConfig(level=logging.DEBUG)

def main(options):
    ms = options.ms
    if ms == '':
            logging.error('You have to specify an input MS, use -h for help')
            return
    cols = options.cols
    incol = options.incol
    dysco = options.dysco
    
    t = pt.table(ms, readonly=False, ack=False)

    for col in cols.split(','):
        if col not in t.colnames():
            logging.info('Adding the output column '+col+' to '+ms+'.')
            if incol == '':
                # prepare col metadata
                cd = t.getcoldesc('DATA')
                coldmi = t.getdminfo('DATA')
                if dysco:
                    cd['dataManagerType'] = 'DyscoStMan'
                    cd['dataManagerGroup'] = 'DyscoData'
                    coldmi = {'NAME': col,'SEQNR': 3,'SPEC': {'dataBitCount': 10,'distribution': 'TruncatedGaussian','distributionTruncation': 2.5,'normalization': 'AF','studentTNu': 0.0,'weightBitCount': 12},'TYPE': 'DyscoStMan'}
                # not as performing as standard DATA
                else:
                    coldmi["NAME"] = col
                    cd['dataManagerType'] = 'StandardStMan'
                    cd['dataManagerGroup'] = 'SSMVar'
                    coldmi = {'NAME': col,'SEQNR': 0,'SPEC': {'ActualCacheSize': 2,'BUCKETSIZE': 32768,'IndexLength': 799832,'PERSCACHESIZE': 2},'TYPE': 'StandardStMan'}

                cd['comment'] = 'Added by addcol2ms'
                t.addcols(pt.makecoldesc(col, cd), coldmi)

                # if non dysco is done by default
                if options.dysco:
                    logging.warning('Setting '+col+' = 0')
                    pt.taql("update $t set "+col+"=0")

            else:
                # prepare col metadata
                coldmi = t.getdminfo(incol)
                coldmi['NAME'] = col
                cd = t.getcoldesc(incol)

                cd['comment'] = 'Added by addcol2ms'
                t.addcols(pt.makecoldesc(col, cd), coldmi)

                logging.warning('Setting '+col+' = '+incol)
                pt.taql("update $t set "+col+"="+incol)

        else:
            logging.warning('Column '+col+' already exists.')

    t.close()
        
opt = optparse.OptionParser()
opt.add_option('-m','--ms',help='Input MS [no default].',default='')
opt.add_option('-c','--cols',help='Output column, comma separated if more than one [no default].',default='')
opt.add_option('-i','--incol',help='Input column to copy in the output column, otherwise it will be set to 0 [default set to 0].',default='')
opt.add_option('-d','--dysco',help='Enable dysco dataManager for new columns (copied columns always get the same dataManager of the original)',action="store_true",default=False)
options, arguments = opt.parse_args()
main(options)


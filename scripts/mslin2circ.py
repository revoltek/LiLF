#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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

import optparse
import numpy
import sys
import casacore.tables as pt
from casacore.quanta import quantity

def checkfile(inms):
  if inms == '':
     print('Error: give an input MS')
     sys.exit()
  #Check feed polarization
  t = pt.table(inms,ack=False)
  feed = pt.table(t.getkeyword('FEED'),ack=False)
  poltyp=feed.getcell('POLARIZATION_TYPE',0) 
  feed.close()
  t.close() 
  if options.reverse == True:
    if poltyp[0] != 'R' and poltyp[0] != 'L':
       print("WARNING: Data is not from circularly polarized feed but I'm converting a column from circular")
  else: 
    if poltyp[0] != 'X' and poltyp[0] != 'Y':
       print("WARNING: Data is not from linearly polarized feed but I'm converting a column from linear.")

def setupiofiles(inms, outms, incolumn, outcolumn):
  """
  if inms!=outms copy it and then work only on outms
  """
  if outms == None:
     outms = inms

  if inms != outms :
     t = pt.table(inms)
     t.copy(outms, True, True)
     t.close()
     print("Finished copy.")

  # create output column if doesn't exist
  to = pt.table(outms, readonly=False)
  if not outcolumn in to.colnames():
      print("Add column %s" % outcolumn)
      ti = pt.table(inms)
      coldmi = ti.getdminfo(incolumn)
      coldmi['NAME'] = outcolumn
      to.addcols(pt.makecoldesc(outcolumn, ti.getcoldesc(incolumn)), coldmi)
      pt.taql("update $to set "+outcolumn+"="+incolumn)
      ti.close()
  to.close()
  return outms


def mslin2circ(incol, outcol, outms, skipmetadata):
  tc = pt.table(outms, readonly=False, ack=False)
  dataXY = tc.getcol(incol)
  I=numpy.complex(0.0,1.0)
  dataRL = 0.5* numpy.transpose(numpy.array([
           +dataXY[:,:,0]-I*dataXY[:,:,1]+I*dataXY[:,:,2]+dataXY[:,:,3],
           +dataXY[:,:,0]+I*dataXY[:,:,1]+I*dataXY[:,:,2]-dataXY[:,:,3],
           +dataXY[:,:,0]-I*dataXY[:,:,1]-I*dataXY[:,:,2]-dataXY[:,:,3],
           +dataXY[:,:,0]+I*dataXY[:,:,1]-I*dataXY[:,:,2]+dataXY[:,:,3]]),
           (1,2,0))
  tc.putcol(outcol,dataRL)

  #Change metadata information to be circular feeds
  if not skipmetadata:
    feed = pt.table(tc.getkeyword('FEED'),readonly=False,ack=False)
    for tpart in feed.iter('ANTENNA_ID'):
        tpart.putcell('POLARIZATION_TYPE',0,['R','L'])

    polariz = pt.table(tc.getkeyword('POLARIZATION'),readonly=False,ack=False)
    polariz.putcell('CORR_TYPE',0,[5,6,7,8])

  tc.close()

def mscirc2lin(incol, outcol, outms, skipmetadata):
  tc = pt.table(outms,readonly=False, ack=False)
  dataRL = tc.getcol(incol)
  I=numpy.complex(0.0,1.0)
  dataXY = 0.5* numpy.transpose(numpy.array([
              +dataRL[:,:,0]+dataRL[:,:,1]+dataRL[:,:,2]+dataRL[:,:,3],
           I*(+dataRL[:,:,0]-dataRL[:,:,1]+dataRL[:,:,2]-dataRL[:,:,3]),
           I*(-dataRL[:,:,0]-dataRL[:,:,1]+dataRL[:,:,2]+dataRL[:,:,3]),
              +dataRL[:,:,0]-dataRL[:,:,1]-dataRL[:,:,2]+dataRL[:,:,3] ]),
           (1,2,0))
  tc.putcol(outcol,dataXY)

  #Change metadata information to be circular feeds
  if not skipmetadata:
    feed = pt.table(tc.getkeyword('FEED'),readonly=False, ack=False)
    for tpart in feed.iter('ANTENNA_ID'):
        tpart.putcell('POLARIZATION_TYPE',0,['X','Y'])

    polariz = pt.table(tc.getkeyword('POLARIZATION'),readonly=False, ack=False)
    polariz.putcell('CORR_TYPE',0,[9,10,11,12])

  tc.close()


def mergeweights(outms):
  """
  Merge weights (weights become the average across the 4 polarizations)
  """
  print("WARNING: updating weights, cannot reverse to original.")
  tc = pt.table(outms,readonly=False, ack=False)
  weights = tc.getcol('WEIGHT_SPECTRUM')
  shape = weights.shape
  # find the mean along the pol axis and then expand the array
  weights = numpy.repeat(numpy.mean(weights, axis=2), 4, axis=1).reshape(shape)
  tc.putcol('WEIGHT_SPECTRUM',weights)
  tc.close()


def mergeflags(outms):
  """
  Merge flags (if a pol is flagged, flag everything)
  """
  #taql("UPDATE $outms set FLAG=True where any(FLAG)")
  tc = pt.table(outms,readonly=False, ack=False)
  flag = tc.getcol('FLAG')
  print("Initial flags:", numpy.count_nonzero(flag))
  shape = flag.shape
  # find if any data is flagged along the pol axis and then expand the array
  flag = numpy.repeat( numpy.any(flag, axis=2), 4, axis=1).reshape(shape)
  print("Final flags:", numpy.count_nonzero(flag))
  tc.putcol('FLAG',flag)
  tc.close()


def updatehistory(outms):
  """
  Update history to show that this script has modified original data
  """
  tc = pt.table(outms,readonly=False)
  th = pt.table(tc.getkeyword('HISTORY'), readonly=False, ack=False)
  nr=th.nrows()
  th.addrows(1)
  tr=th.row()
  tr.put(nr,{'TIME': quantity('today').get('s').get_value(), 'OBSERVATION_ID':0,'MESSAGE': ' ', 'PRIORITY': ' ', 'ORIGIN': ' ','OBJECT_ID':0, 'APPLICATION':'mslin2circ','CLI_COMMAND':[''],'APP_PARAMS': ['']})


opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS (format: ms:COLUMN, default column: DATA)',default='')
opt.add_option('-o','--outms',help='Output MS (format: ms:COLUMN, default ms: InputMS, default column: DATA)',default='')
opt.add_option('-r','--reverse',action="store_true",default=False,help='Convert from circular to linear')
opt.add_option('-s','--skipmetadata',action="store_true",default=False,help='Skip setting the metadata correctly')
opt.add_option('-w','--weights',action="store_true",default=False,help='Weights are updated to reflect the combined polarization (cannot be undone with -r)')
options, arguments = opt.parse_args()

if options.outms == '':
    options.outms = options.inms.split(':')[0]

if len( options.inms.split(':') ) == 2:
    incolumn = options.inms.split(':')[1]
else:
    incolumn = 'DATA'

if len( options.outms.split(':') ) == 2:
    outcolumn = options.outms.split(':')[1]
else:
    outcolumn = 'DATA'

inms = options.inms.split(':')[0]
outms = options.outms.split(':')[0]
checkfile(inms)
outms = setupiofiles(inms, outms, incolumn, outcolumn)

print("INFO: inms: "+inms+" (column: "+incolumn+")")
print("INFO: outms: "+outms+" (column: "+outcolumn+")")

if options.reverse == True:
   mscirc2lin(incolumn, outcolumn, outms, options.skipmetadata)
else:
   mslin2circ(incolumn, outcolumn, outms, options.skipmetadata)
if options.weights: mergeweights(outms)
mergeflags(outms)
updatehistory(outms)

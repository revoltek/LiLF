#!/usr/bin/env python

import os, sys, glob
import numpy as np
from astropy.io import fits

def make_beam(imageroot, outfile='beam.fits'):

    imagefiles = sorted(glob.glob(imageroot+'-[0-9]*-image.fits'))
    beamfiles = sorted(glob.glob(imageroot+'-[0-9]*-beam-I.fits'))

    if len(beamfiles) < 2 or len(imagefiles) < 2:
        print('Need at least 2 beam images.')
        sys.exit()

    # get weight
    weights = []
    for imagefile in imagefiles:
        with fits.open(imagefile) as imagefits:
            header = imagefits[0].header
            weights.append(header['WSCIMGWG'])
    #print('Weights:', weights)

    # combine beams
    with fits.open(beamfiles[0]) as beamfits:
        beam = beamfits[0].data * weights[0]
        header = beamfits[0].header
    for i, beamfile in enumerate(beamfiles[1:]):
        with fits.open(beamfile) as beamfits:
            beam += beamfits[0].data * weights[i+1]

    beam /= np.sum(weights)
    print("Writing: "+outfile)
    fits.writeto(outfile, beam, header, overwrite=True)

if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] image-root-name \n Francesco de Gasperin', version='1.0')
    opt.add_option('-o', '--outfile', help='Filename for the combined beam fits (default=beam.fits)', default='beam.fits')
    opt.add_option('-i', '--infile', help='Root name of the beam images (default=None)', default=None)
    (options, args) = opt.parse_args()

    outfile = options.outfile
    infile = options.infile

    if infile is None:
        print('Missing filename.')
        sys.exit()

    make_beam(infile, outfile)

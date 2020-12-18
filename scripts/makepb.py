#!/usr/bin/env python

# use idg to create a primary beam file for LOFAR
# the output of wsclean is re-arranged to have a single beam file even if more channels are used

import os, sys, glob
from astropy.io import fits

def make_beam(mss, outfile='beam.fits', pixscale=10, size=5, nchans=1):

    for ms in mss:
        os.system('DPPP msin=%s msin.datacolumn=CORRECTED_DATA msout=. msout.datacolumn=DATA steps=[setbeam] \
        setbeam.type=setbeam setbeam.beammode=default' % ms)

    # make beam image with wsclean
    mss_string = ' '.join(mss)
    size = size*3600/pixscale
    cmd = 'wsclean -j 64 -data-column DATA -parallel-reordering 4 -name __beam -use-idg -idg-mode cpu -channels-out 1 \
            -grid-with-beam -use-differential-lofar-beam -beam-aterm-update 400 \
            -size %i %i -scale %farcsec -niter 0 -no-update-model-required -pol I %s' \
            % (int(size), int(size), float(pixscale), mss_string)
    print(cmd)
    os.system(cmd)

    if nchans == 1:
        print("Writing: "+outfile)
        os.system('mv __beam-beam-I.fits %s' % outfile)
        os.system('rm __beam*')
    else:
        imagefiles = sorted(glob.glob('__beam-[0-9]*-image.fits'))
        beamfiles = sorted(glob.glob('__beam-[0-9]*-beam-I.fits'))

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
        os.system('rm __beam*')


if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] image-root-name \n Francesco de Gasperin', version='1.0')
    opt.add_option('-o', '--outfile', help='Filename for the beam fits (default=beam.fits)', default='beam.fits')
    opt.add_option('-f', '--fromimg', help='Image to use to get size/pixscale (default=None)', default=None)
    opt.add_option('-s', '--size', help='Image size (deg) (default=5)', default=5)
    opt.add_option('-p', '--pixscale', help='Pixelscale (arcsec) (default=10)', default=10)
    opt.add_option('-n', '--nchans', help='Number of output channels to create the beam (default=1)', default=1)
    (options, args) = opt.parse_args()

    outfile = options.outfile
    fromimg = options.fromimg
    size = int(options.size)
    pixscale = float(options.pixscale)
    nchans = int(options.nchans)
    mss = args

    if len(mss) == 0:
        print('Missing measurementsets.')
        sys.exit()

    if fromimg is not None:
        header = fits.getheader(fromimg)
        size = header['NAXIS1']
        pixscale = header['CDELT1']*3600

    make_beam(mss, outfile, pixscale, size, nchans)

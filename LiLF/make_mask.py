#!/usr/bin/env python

# create a mask using bdsm of an image

def make_mask(image_name, mask_name=None, threshisl=5, atrous_do=False, rmsbox=(100,10), adaptive_thresh=50, mask_combine=None, write_srl=False):

    import sys, os
    import numpy as np
    from astropy.io import fits as pyfits
    import bdsf

    # wavelets are required to fit gaussians
    if atrous_do or write_srl: stop_at = None
    else: stop_at = 'isl'

    # DO THE SOURCE DETECTION
    img = bdsf.process_image(image_name, rms_box=rmsbox, frequency=54e6, \
        thresh_isl=float(threshisl), thresh_pix=float(threshisl*5/3.), rms_map=True, mean_map='zero', atrous_do=atrous_do, atrous_jmax=4, \
        adaptive_rms_box=True, adaptive_thresh=adaptive_thresh, rms_box_bright=(30,5), \
        flagging_opts=True, flag_maxsize_fwhm=0.5, stop_at=stop_at, quiet=True, debug=False)

    # WRITE THE MASK FITS
    if mask_name == None: mask_name = image_name+'.newmask'
    if os.path.exists(mask_name): os.system('rm -r ' + mask_name)
    img.export_image(img_type='island_mask', img_format='fits', outfile=mask_name, clobber=True)

    # WRITE CATALOGUE
    if write_srl:
        img.write_catalog(format='bbs', catalog_type='gaul', outfile=mask_name.replace('fits','skymodel'), clobber=True)

    del img

    # do an pixel-by-pixel "OR" operation with a given mask
    if not mask_combine is None:
        print("Doing a pix-by-pix OR with %s." % mask_combine)
        with pyfits.open(mask_combine) as fits:
            data_comb = fits[0].data
        with pyfits.open(mask_name) as fits:
            data = fits[0].data
            assert data.shape() == data_comb.shape()
            data[(data_comb == 1.)] = 1.
            fits[0].data = data
            fits.writeto(mask_name, overwrite=True)

    return mask_name

if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] imagename \n Francesco de Gasperin', version='1.0')
    opt.add_option('-i', '--threshisl', help='Threshold island (default=3)', type='int', default=3)
    opt.add_option('-t', '--atrous_do', help='BDSM extended source detection (default=False)', action='store_true', default=False)
    opt.add_option('-m', '--newmask', help='Mask name (default=imagename with mask in place of image)', default=None)
    opt.add_option('-c', '--combinemask', help='Mask name of a mask to add to the found one (default=None)', default=None)
    opt.add_option('-r', '--rmsbox', help='rms box size (default=100,10)', default='100,10')
    opt.add_option('-t', '--atrous_do', help='BDSM extended source detection (default=False)', action='store_true', default=False)
    opt.add_option('-s', '--write_srl', help='Write SRL (default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()
    
    rmsbox = (int(options.rmsbox.split(',')[0]),int(options.rmsbox.split(',')[1]))
    make_mask(args[0].rstrip('/'), options.newmask, options.threshisl, options.atrous_do, rmsbox, options.combinemask, options.write_srl)

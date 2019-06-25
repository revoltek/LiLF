import os, sys, glob
import numpy as np
import astropy.io.fits as pyfits
import lsmtool
import pyregion
from scipy.ndimage.measurements import label
from LiLF import make_mask, lib_util
from LiLF.lib_log import logger

class Image(object):
    def __init__(self, imagename, userReg = None, beamReg= None ):
        """
        userMask: keep this region when making masks
        BeamReg: ds9 region file of the beam
        """
        if 'MFS' in imagename: suffix = 'MFS-image.fits'
        else: suffix = 'image.fits'

        self.imagename    = imagename
        self.root         = imagename.replace(suffix, '')
        self.maskname     = imagename.replace(suffix, 'mask.fits')
        self.skymodel     = imagename.replace(suffix, 'sources.txt')
        self.skymodel_cut = imagename.replace(suffix, 'sources-cut.txt')
        self.skydb        = imagename.replace(suffix, 'sources-cut.skydb')
        self.userReg      = userReg
        self.beamReg      = beamReg

    def rescaleModel(self, funct_flux):
        """
        Rescale the model images to a certain total flux estimated using funct_flux(nu)
        nu is taken from the fits file.

        funct_flux: is a function of frequency (Hz) which returns the total expected flux (Jy) at that frequency.
        """
        for model_img in sorted(glob.glob(self.root+'*model*.fits')):
            fits = pyfits.open(model_img)
            # get frequency
            assert fits[0].header['CTYPE3'] == 'FREQ'
            nu = fits[0].header['CRVAL3']
            data = fits[0].data
            # find expected flux
            flux = funct_flux(nu)
            current_flux = np.sum(data)
            # rescale data
            scaling_factor = flux/current_flux
            logger.warning('Rescaling model %s by: %f' % (model_img, scaling_factor))
            data *= scaling_factor
            fits[0].data = data
            fits.writeto(model_img, overwrite=True)
            fits.close()


    def makeMask(self, threshisl=5, atrous_do=True, remove_extended_cutoff=0., maskname=None):
        """
        Create a mask of the image where only believable flux is

        remove_extended_cutoff: if >0 then remove all islands where sum(brightness_pixels)/(#pixels^2) < remove_extended_cutoff
        this is useful to remove extended sources from the mask. This higher this number the more compact must be the source.
        A good value is 0.001 for DIE cal images.

        maskname: if give, then use a specific maskname
        """
        if maskname is None: maskname = self.maskname

        if not os.path.exists(maskname):
            logger.info('%s: Making mask...' % self.imagename)
            make_mask.make_mask(image_name=self.imagename, mask_name=maskname, threshisl=threshisl, atrous_do=atrous_do)
        if self.userReg is not None:
            logger.info('%s: Adding user mask (%s)...' % (self.imagename, self.userReg))
            blank_image_reg(maskname, self.userReg, inverse=False, blankval=1)

        if remove_extended_cutoff > 0:

            # get data
            with pyfits.open(self.imagename) as fits:
                data = fits[0].data
            # get mask
            with pyfits.open(maskname) as fits:
                mask = fits[0].data
                # for each island calculate the catoff
                blobs, number_of_blobs = label(mask.astype(int).squeeze(), structure=[[1,1,1],[1,1,1],[1,1,1]])
                for i in range(1,number_of_blobs):
                    this_blob = blobs == i
                    max_pix = np.max(data[0,0,this_blob])
                    ratio = np.sum(data[0,0,this_blob])/np.sum(mask[0,0,this_blob])**2
                    if max_pix < 1. and ratio < remove_extended_cutoff:
                        mask[0,0,this_blob] = False
                    #mask[0,0,this_blob] = ratio # debug

                # write mask back
                fits[0].data = mask
                fits.writeto(maskname, overwrite=True)


    def selectCC(self, keepInBeam=True):
        """
        remove cc from a skymodel according to masks
        keepInBeam: if beamReg is present and is True: remove sources outside beam
                    if beamReg is present and is False: remove source inside beam
        """
        self.makeMask()

        if self.beamReg is not None:
            logger.info('Predict (apply beam reg %s)...' % self.beamReg)
            blank_image_reg(self.maskname, self.beamReg, inverse=keepInBeam, blankval=0) # if keep_in_beam set to 0 everything outside beam.reg

        # apply mask
        logger.info('%s: Apply mask on skymodel...' % self.imagename)
        lsm = lsmtool.load(self.skymodel)
        lsm.select('%s == True' % self.maskname)
        lsm.group('single') # group to 1 patch
        lsm.write(self.skymodel_cut, format = 'makesourcedb', clobber=True)
        del lsm

        # convert from txt to blob
        logger.info('%s: Make skydb...' % self.imagename)
        lib_util.check_rm(self.skydb)
        os.system('makesourcedb outtype="blob" format="<" in="'+self.skymodel_cut+'" out="'+self.skydb+'"')


    def getNoise(self, boxsize=None):
        """
        Return the rms of all the non-masked pixels in an image
        boxsize : limit to central box of this pixelsize
        """   
        self.makeMask()

        with pyfits.open(self.imagename) as fits:
            with pyfits.open(self.maskname) as mask:
                data = fits[0].data
                mask = mask[0].data
                if boxsize is not None:
                    if len(data.shape)==4:
                        _,_,ys,xs = data.shape
                        data = data[0,0,ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
                        mask = mask[0,0,ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
                    else:
                        ys,xs = data.shape
                        data = data[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
                        mask = mask[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
    
                return np.nanstd(data[mask==0])



def flatten(f, channel = 0, freqaxis = 0):
    """
    Flatten a fits file so that it becomes a 2D image. Return new header and data
    """
    from astropy import wcs

    naxis=f[0].header['NAXIS']
    if (naxis < 2):
        raise RadioError("Can\'t make map from this")
    if (naxis == 2):
        return f[0].header,f[0].data

    w               = wcs.WCS(f[0].header)
    wn              = wcs.WCS(naxis = 2)

    wn.wcs.crpix[0] = w.wcs.crpix[0]
    wn.wcs.crpix[1] = w.wcs.crpix[1]
    wn.wcs.cdelt    = w.wcs.cdelt[0:2]
    wn.wcs.crval    = w.wcs.crval[0:2]
    wn.wcs.ctype[0] = w.wcs.ctype[0]
    wn.wcs.ctype[1] = w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = f[0].header['NAXIS1']
    header["NAXIS2"] = f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r = f[0].header.get(k)
        if (r):
            header[k] = r

    slice = []
    for i in range(naxis,0,-1):
        if (i <= 2):
            slice.append(np.s_[:],)
        elif (i == freqaxis):
            slice.append(channel)
        else:
            slice.append(0)

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header, f[0].data[slice]


def blank_image_fits(filename, maskname, outfile = None, inverse = False, blankval = 0.):
    """
    Set to "blankval" all the pixels inside the given region
    if inverse=True, set to "blankval" pixels outside region.

    filename: fits file
    region: ds9 region
    outfile: output name
    inverse: reverse region mask
    blankval: pixel value to set
    """

    if (outfile == None):
        outfile = filename

    with pyfits.open(maskname) as fits:
        mask = fits[0].data
    
    if (inverse): mask = ~(mask.astype(bool))

    with pyfits.open(filename) as fits:
        data = fits[0].data

        assert mask.shape == data.shape # mask and data should be same shape

        sum_before = np.sum(data)
        data[mask] = blankval
        logger.debug("%s: Blanking (%s): sum of values: %f -> %f" % (filename, maskname, sum_before, np.sum(data)))
        fits.writeto(outfile, overwrite=True)

 
def blank_image_reg(filename, region, outfile = None, inverse = False, blankval = 0., op = "AND"):
    """
    Set to "blankval" all the pixels inside the given region
    if inverse=True, set to "blankval" pixels outside region.
    If a list of region is provided the operation is applied to each region one after the other

    filename: fits file
    region: ds9 region or list of regions
    outfile: output name
    inverse: reverse final *combined* mask
    blankval: pixel value to set
    op: how to combine multiple regions with AND or OR
    """

    if outfile == None: outfile = filename
    if not type(region) is list: region=[region]

    # open fits
    with pyfits.open(filename) as fits:
        origshape    = fits[0].data.shape
        header, data = flatten(fits)
        sum_before   = np.sum(data)
        if (op == 'AND'):
            total_mask = np.ones(shape = data.shape).astype(bool)
        if (op == 'OR'):
            total_mask = np.zeros(shape = data.shape).astype(bool)
        for this_region in region:
            # extract mask
            r    = pyregion.open(this_region)
            mask = r.get_mask(header=header, shape=data.shape)
            if (op == 'AND'):
                total_mask = total_mask & mask
            if (op == 'OR'):
                total_mask = total_mask | mask
        if (inverse):
            total_mask = ~total_mask
        data[total_mask] = blankval
        # save fits
        fits[0].data = data.reshape(origshape)
        fits.writeto(outfile, overwrite=True)

    logger.debug("%s: Blanking (%s): sum of values: %f -> %f" % (filename, region, sum_before, np.sum(data)))

def make_fits(filename, shape, fill_value=1):
    """
    Create a fits file
    """
    data = np.full(shape=shape, fill_value=fill_value)
    hdu = pyfits.PrimaryHDU(data)
    hdul = pyfits.HDUList([hdu])
    hdul.writeto(filename, overwrite=True)

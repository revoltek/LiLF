import os, sys
import numpy as np

# Load external script make_mask!

import logging
logger = logging.getLogger("PiLL")

class Image(object):
    def __init__(self, imagename, region_facet = None, user_mask = None):
        self.imagename    = imagename
        self.maskname     = imagename.replace('MFS-image.fits', 'mask.fits')
        self.skymodel     = imagename.replace('MFS-image.fits', 'sources.txt')
        self.skymodel_cut = imagename.replace('MFS-image.fits', 'sources-cut.txt')
        self.region_facet = region_facet
        self.user_mask    = user_mask

    def make_mask(self, threshisl=5):
        """
        Create a mask of the image where only believable flux is
        """
        logger.info('%s: Making mask...' % self.imagename)
        if not os.path.exists(self.maskname):
            make_mask(image_name=self.imagename, mask_name=self.maskname, threshisl=threshisl, atrous_do=True)
        if self.user_mask is not None:
            logger.info('%s: Adding user mask (%s)...' % (self.imagename, self.user_mask))
            blank_image_reg(self.maskname, self.user_mask, inverse=False, blankval=1)

    def select_cc(self):
        """
        remove cc from a skymodel according to masks
        """
        self.make_mask()

        if self.region_facet is not None:
            logger.info('Predict (apply facet mask %s)...' % self.region_facet)
            blank_image_reg(self.maskname, self.region_facet, inverse=True, blankval=0) # set to 0 pixels outside facet mask

        # apply mask
        logger.info('%s: Apply mask on skymodel...' % self.imagename)
        lsm = lsmtool.load(self.skymodel)
        lsm.select('%s == True' % self.maskname)
        lsm.write(self.skymodel_cut, format = 'makesourcedb', clobber=True)
        del lsm


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


def size_from_reg(filename, regions, coord):
    """
    find number of pixels necessary to cover the intersection of the regions with a square centered on coord
    filename : a fits file
    regions : ds9 regions
    coord : coordinate of the image center
    """
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import pyregion

    fits         = pyfits.open(filename)
    header, data = flatten(fits)
    w            = pywcs.WCS(header)

    # find max dist in pixel on reference image
    x_c, y_c     = w.all_world2pix(coord[0], coord[1], 0, ra_dec_order=True)
    #print x_c, y_c
    
    # ditance would overestimate, get max of x-x_c and y-y_c
    mask = np.ones(shape=data.shape, dtype=bool)
    for region in regions:
        r    = pyregion.open(region)
        mask = (mask & r.get_mask(header=header, shape=data.shape))
    y, x = mask.nonzero()

    if (len(x) == 0):
        return 0

    max_size = 2 * np.max([np.abs(np.array(y) - y_c), np.abs(np.array(x) - x_c)])

    return int(max_size)

 
def scale_from_ms(ms):
    """
    Get the pixel scale in arcsec for a full-res image.
    It is 1/3 of the max resolution assuming zenit observation.
    Completely flagged lines are removed
    """
    from pyrap.tables import table
    import numpy as np

    c          = 299792458.                                   # in metres per second

    t          = table(ms, ack = False).query('not all(FLAG)')
    col        = t.getcol('UVW')
    t.close()

    t          = table(ms+'/SPECTRAL_WINDOW', ack=False)
    wavelength = c / t.getcol('REF_FREQUENCY')[0]             # in metres
    t.close()
    #print 'Wavelength:', wavelength,'m (Freq: '+str(t.getcol('REF_FREQUENCY')[0]/1.e6)+' MHz)'

    maxdist    = np.max( np.sqrt(col[:,0] ** 2 + col[:,1] ** 2) )

    return int(round(wavelength / maxdist * (180 / np.pi) * 3600 / 3.)) # in arcseconds


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
    import astropy.io.fits as pyfits

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
        fits.writeto(outfile, clobber=True)

 
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
    import astropy.io.fits as pyfits
    import pyregion

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
        fits.writeto(outfile, clobber=True)

    logger.debug("%s: Blanking (%s): sum of values: %f -> %f" % (filename, region, sum_before, np.sum(data)))


def get_noise_img(filename, boxsize=None, niter=20, eps=1e-5):
    """
    Return the rms of all the pixels in an image
    boxsize : limit to central box of this pixelsize
    niter : robust rms estimation
    eps : convergency
    """
    import astropy.io.fits as pyfits
    
    with pyfits.open(filename) as fits:
        data = fits[0].data
        if boxsize is None:
            subim = data
        else:
           if len(data.shape)==4:
                _,_,ys,xs = data.shape
                subim = data[0,0,ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
           else:
                ys,xs = data.shape
                subim = data[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2].flatten()
        oldrms = 1.
        for i in range(niter):
            rms = np.nanstd(subim)
            #print len(subim),rms
            if np.abs(oldrms-rms)/rms < eps:
                return rms
            subim=subim[np.abs(subim)<5*rms]
            oldrms=rms
        raise Exception('Failed to converge')

#def nan2zeros(filename):
#    """
#    Replace NaNs to zeros in a fits file
#    """
#    import astropy.io.fits as pyfits
#    with pyfits.open(filename) as fits:
#        fits[0].data = np.nan_to_num(fits[0].data)
#        fits.writeto(filename, clobber=True)
#
#
#def get_coord_centroid(filename, region):
#    """
#    Get centroid coordinates from an image and a region
#    filename: fits file
#    region: ds9 region
#    """
#    import astropy.io.fits as pyfits
#    import astropy.wcs as pywcs
#    import pyregion
#    from scipy.ndimage.measurements import center_of_mass
#
#    fits = pyfits.open(filename)
#    header, data = flatten(fits)
#
#    # extract mask and find center of mass
#    r = pyregion.open(region)
#    mask = r.get_mask(header=header, shape=data.shape)
#    dec_pix, ra_pix = center_of_mass(mask)
#    
#    # convert to ra/dec in angle
#    w = pywcs.WCS(fits[0].header)
#    #w = w.celestial # needs newer astropy
#    ra, dec = w.all_pix2world(ra_pix, dec_pix, 0, 0, 0, ra_dec_order=True)
#
#    fits.close()
#    return float(ra), float(dec)
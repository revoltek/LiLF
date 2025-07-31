import os, glob
import numpy as np
import astropy.io.fits as pyfits
import casacore.images as pim
from casacore import quanta
import lsmtool
import pyregion
from scipy.ndimage.measurements import label
from astropy.stats import median_absolute_deviation
from LiLF import make_mask, lib_util
from LiLF.lib_log import logger
import astropy.io.fits as fits
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.visualization import AsymmetricPercentileInterval, SqrtStretch, LogStretch, ImageNormalize

class Image(object):
    def __init__(self, imagename, userReg = None, beamReg= None ):
        """
        userReg: keep this region when making masks
        BeamReg: ds9 region file of the beam
        """
        if 'MFS' in imagename: suffix = '-MFS-image.fits'
        elif 'image.fits' in imagename: suffix = '-image.fits'
        elif 'restored.fits' in imagename: suffix = '.app.restored.fits'
        else: suffix = '.fits'
        if userReg == '': userReg = None
        if beamReg == '': beamReg = None

        self.imagename    = imagename
        self.root         = imagename.replace(suffix, '')
        self.maskname     = imagename.replace(suffix, '-mask.fits')
        self.skymodel     = imagename.replace(suffix, '-sources.txt')
        self.skymodel_cut = imagename.replace(suffix, '-sources-cut.txt')
        self.skydb        = imagename.replace(suffix, '-sources-cut.skydb')
        self.userReg      = userReg
        self.beamReg      = beamReg

    def calc_flux(self, mask):
        """
        Get flux inside given region. Adapted from Martin Hardcastle's radiomap class
        """
        extract = mask

        with fits.open(self.imagename) as phdu:
            head, lhdu = flatten(phdu)
            gfactor = 2.0 * np.sqrt(2.0 * np.log(2.0))
            f = phdu[0]
            prhd = phdu[0].header
            units = prhd.get('BUNIT')
            if units is None:
                units = prhd.get('UNIT')
            if units != 'JY/BEAM' and units != 'Jy/beam':
                print('Warning: units are', units, 'but code expects JY/BEAM')
            bmaj = prhd.get('BMAJ')
            bmin = prhd.get('BMIN')

            bmaj = np.abs(bmaj)
            bmin = np.abs(bmin)

            w = self.getWCS()
            cd1 = -w.wcs.cdelt[0]
            cd2 = w.wcs.cdelt[1]
            if ((cd1 - cd2) / cd1) > 1.0001 and ((bmaj - bmin) / bmin) > 1.0001:
                print('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

            bmaj /= cd1
            bmin /= cd2
            area = 2.0 * np.pi * (bmaj * bmin) / (gfactor * gfactor)

            d = [lhdu]

            region = pyregion.open(extract).as_imagecoord(prhd)

            for i, n in enumerate(d):
                mask = region.get_mask(hdu=f, shape=np.shape(n))
                data = np.extract(mask, d)
                nndata = data[~np.isnan(data)]
                flux = np.sum(nndata) / area

            return flux

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

    def nantozeroModel(self):
        """
        Set nan/inf to 0 in all model images
        """
        for modelimage in sorted(glob.glob(self.root + '*model*.fits')):
            with pyfits.open(modelimage, mode='update') as fits:
                for hdu in fits:
                    if not np.isfinite(hdu.data).all():
                        logger.info(f"Model image '{modelimage}' has '{np.sum(~np.isfinite(hdu.data))}' NaN/inf values (NaN/inf -> zeros)")
                    hdu.data[~np.isfinite(hdu.data)] = 0

    # TODO: separate makemask (using breizorro) and makecat (using bdsf)
    def makeMask(self, threshpix=5, atrous_do=False, rmsbox=(100,10), remove_extended_cutoff=0., only_beam=False, maskname=None,
                 write_srl=False, write_gaul=False, write_ds9=False, mask_combine=None):
        """
        Create a mask of the image where only believable flux is

        remove_extended_cutoff: if >0 then remove all islands where sum(brightness_pixels)/(#pixels^2) < remove_extended_cutoff
        this is useful to remove extended sources from the mask. This higher this number the more compact must be the source.
        A good value is 0.001 for DIE cal images.

        maskname: if give, then use a specific maskname
        only_beam: set to 0 outside the beam
        """
        if maskname is None: maskname = self.maskname

        if not os.path.exists(maskname):
            logger.info('%s: Making mask (%s)...' % (self.imagename, maskname))
            make_mask.make_mask(image_name=self.imagename, mask_name=maskname, threshpix=threshpix, atrous_do=atrous_do,
                                rmsbox=rmsbox, write_srl=write_srl, write_gaul=write_gaul, write_ds9=write_ds9, mask_combine=mask_combine)

        if remove_extended_cutoff > 0:

            # get data
            with pyfits.open(self.imagename) as fits:
                data = np.squeeze(fits[0].data)

            # get mask
            with pyfits.open(maskname) as fits:
                mask = np.squeeze(fits[0].data)

                # for each island calculate the catoff
                blobs, number_of_blobs = label(mask.astype(int).squeeze(), structure=[[1,1,1],[1,1,1],[1,1,1]])
                for i in range(1,number_of_blobs):
                    this_blob = blobs == i
                    max_pix = np.max(data[this_blob])
                    ratio = np.sum(data[this_blob])/np.sum(mask[this_blob])**2
                    if max_pix < 1. and ratio < remove_extended_cutoff:
                        mask[this_blob] = False
                    #mask[0,0,this_blob] = ratio # debug

                # write mask back
                fits[0].data[0,0] = mask
                fits.writeto(maskname, overwrite=True)

        if self.userReg is not None:
            logger.info('%s: Adding user mask (%s)...' % (self.imagename, self.userReg))
            blank_image_reg(maskname, self.userReg, inverse=False, blankval=1)

        if only_beam and self.beamReg is not None:
            logger.info('%s: Restricting to the beam (%s)...' % (self.imagename, self.beamReg))
            blank_image_reg(maskname, self.beamReg, inverse=True, blankval=0)


    def selectCC(self, checkBeam=True, keepInBeam=True, maskname=None):
        """
        remove cc from a skymodel according to masks
        checkBeam: remove according to beam (see keepInBeam)
        keepInBeam: if beamReg is present and is True: remove sources outside beam
                    if beamReg is present and is False: remove source inside beam
        maskname: a possible mask, otherwise try to use standard
        """
        if maskname is None: maskname = self.maskname
        if not os.path.exists(maskname):
            raise("Missing mask in selectCC: %s." % maskname)

        if checkBeam:
            if self.beamReg is None:
                raise('Missing beam in selectCC.')
            logger.info('Predict (apply beam reg %s)...' % self.beamReg)
            blank_image_reg(maskname, self.beamReg, inverse=keepInBeam, blankval=0) # if keep_in_beam set to 0 everything outside beam.reg

        # apply mask
        logger.info('%s: Apply mask (%s) on skymodel...' % (self.imagename,maskname))
        lsm = lsmtool.load(self.skymodel)
        lsm.select('%s == True' % maskname)
        lsm.group('single') # group to 1 patch
        lsm.write(self.skymodel_cut, format = 'makesourcedb', clobber=True)
        del lsm

        # convert from txt to blob
        logger.info('%s: Make skydb...' % self.imagename)
        lib_util.check_rm(self.skydb)
        os.system('makesourcedb outtype="blob" format="<" in="'+self.skymodel_cut+'" out="'+self.skydb+'"')

    def getNoise(self, boxsize=None, useMask=True):
        """
        Return the rms of all the non-masked pixels in an image
        boxsize : limit to central box of this pixelsize
        """   
        if useMask:
            self.makeMask()

            with pyfits.open(self.imagename) as fits:
                with pyfits.open(self.maskname) as mask:
                    data = np.squeeze(fits[0].data)
                    mask = np.squeeze(mask[0].data)
                    if boxsize is not None:
                        ys,xs = data.shape
                        data = data[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2]
                        mask = mask[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2]
    
            return np.nanstd(data[mask==0])
        else:
            with pyfits.open(self.imagename) as fits:
                data = np.squeeze(fits[0].data)
                if boxsize is not None:
                    ys,xs = data.shape
                    data = data[ys/2-boxsize/2:ys/2+boxsize/2,xs/2-boxsize/2:xs/2+boxsize/2]
            
            eps = 1e-3; niter = 100; sigma = 5
            data = data[ ~np.isnan(data) & (data != 0) ] # remove nans and 0s
            initial_len = len(data)
            if initial_len == 0: return 0
            mad_old = 0.
            for i in range(niter):
                mad = median_absolute_deviation(data)
                logger.debug('%s: MAD noise: %.1f mJy on %f%% data' % (self.imagename, mad*1e3, 100*len(data)/initial_len))
                if np.isnan(mad): return 0
                if np.abs(mad_old-mad)/mad < eps:
                    rms = np.nanstd( data )
                    logger.debug('%s: Noise: %.1f mJy/b (data len: %i -> %i - %.2f%%)' % (self.imagename, rms*1e3, initial_len, len(data), 100*len(data)/initial_len))
                    return rms

                data = data[np.abs(data) < (sigma*mad)]
                mad_old = mad

    def getMaxMinRatio(self):
        """
        Return the ratio of the max over min in the image
        """   
        with pyfits.open(self.imagename) as fits:
            data = np.squeeze(fits[0].data)
            return np.abs(np.max(data)/np.min(data))

    def getBeam(self):
        """
        Return the beam size of the image
        """
        this_pim = pim.image(self.imagename)
        info_dict = this_pim.info()['imageinfo']['restoringbeam']
        # get beam info
        bpar_ma = quanta.quantity(info_dict['major']).get_value('arcsec')
        bpar_mi = quanta.quantity(info_dict['minor']).get_value('arcsec')
        bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
        #print('\n{0} - Beam: maj {1:0.3f} (arcsec), min {2:2.3f} (arcsec), pa {3:0.2f} (deg)'.format(img, bpar_ma, bpar_mi,bpar_pa))
        return (bpar_ma,bpar_mi,bpar_pa)

    def getFreq(self):
        """
        :return:
        The flux of the image
        """
        with pyfits.open(self.imagename) as fits:
            if fits[0].header['CTYPE3'] == 'FREQ':
                return fits[0].header['CRVAL3']
            elif fits[0].header['CTYPE4'] == 'FREQ':
                return fits[0].header['CRVAL4']
            else:
                raise RuntimeError('Cannot find frequency in image %s' % self.imagename)

    def getWCS(self, flat=True):
        """
        get the WCS from the fits image header
        Returns
        -------
        WCS
        """
        with fits.open(self.imagename) as phdu:
            if flat:
                hdr, data_res = flatten(phdu)
                return wcs.WCS(hdr)
            else:
                return wcs.WCS(phdu[0].header)

    def plotimage(self, outplotname, regionfile=None, minmax=None, \
                      cmap='gist_heat', regioncolor='yellow', stretch_type='sqrt', regionalpha=0.6):
        """
        minmax: in mJy/beam, if None then use 20-99.99 percentile
        """

        # image noise info
        head, data = flatten(self.imagename)
        imagenoiseinfo = self.getNoise(useMask=False)
    
        f = plt.figure(figsize=(8, 8), dpi=450)
        ax = f.add_subplot(111, projection=self.getWCS() ) #, slices=('x', 'y', 0, 0))

        data *= 1e3  # convert to mJy/beam
        if stretch_type == 'sqrt':
            interval = AsymmetricPercentileInterval(20, 99.99)  # 20 - 99.99 percentile
            stretch = SqrtStretch()
        elif stretch_type == 'log':
            interval = AsymmetricPercentileInterval(80, 99.999)  # 80 - 99.99 percentile
            stretch = LogStretch()
        int_min, int_max = interval.get_limits(data)
        if int_max > 20 * imagenoiseinfo*1e3: int_max = 20 * imagenoiseinfo*1e3
        if minmax != None: int_min, int_max = minmax
        norm = ImageNormalize(data, vmin=float(int_min), vmax=float(int_max), stretch=stretch)#, clip=True)
        img = ax.imshow(data, origin='lower', interpolation='nearest', cmap=cmap, norm=norm)
        #img = ax.imshow(data, cmap=cmap, vmax=minmax[1], vmin=minmax[0])
        
        ax.set_title(self.imagename + ' (noise = {} mJy/beam)'.format(round(imagenoiseinfo * 1e3, 3)),fontsize=6)
    
        ax.grid(False)
        ax.coords[0].set_axislabel_position('b') # for some reason this needs to be hardcoded, otherwise RA gets on top axis
        ax.coords[0].set_ticks_position('bt')
        ax.coords[0].set_ticklabel_position('b') # for some reason this needs to be hardcoded, otherwise RA gets on top axis
    
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')
        try:
            from astropy.visualization.wcsaxes import add_beam
            add_beam(ax, header=head, frame=False) 
        except Exception as e:
            print(f"Cannot plot beam on image, failed with error: {e}. Skipping.")

        cbar = plt.colorbar(img, shrink=0.8)
        cbar.set_label('Flux density (mJy beam$^{-1}$)')
        #ax.add_artist(_add_astropy_beam(fitsimagename))
 
        try: 
            if regionfile is not None:
                import regions
                ds9regions = regions.Regions.read(regionfile, format='ds9')
                for ds9region in ds9regions:
                    reg = ds9region.to_pixel(self.getWCS())
                    reg.plot(ax=ax, alpha=1)
        except Exception as e:
            print(f"Cannot overplot facets, failed with error: {e}. Skipping.")
        
        #if os.path.isfile(outplotname + '.png'):
        #    os.system('rm -f ' + outplotname + '.png')
        plt.savefig(outplotname, format='png', bbox_inches='tight')
        plt.close()


def flatten(f, channel = 0, freqaxis = 0):
    """
    Flatten a fits file so that it becomes a 2D image. Return new header and data
    f: fits file or fits HDUList
    """
    if isinstance(f, str):
        f = pyfits.open(f)

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

    slicing = []
    for i in range(naxis,0,-1):
        if (i <= 2):
            slicing.append(np.s_[:],)
        elif (i == freqaxis):
            slicing.append(channel)
        else:
            slicing.append(0)

    try:
        # if the beam is present, we need to keep it
        header["BMAJ"] = f[0].header["BMAJ"]
        header["BMIN"] = f[0].header["BMIN"]
        header["BPA"] = f[0].header["BPA"]
    except KeyError:
        # if the beam is not present, we do not need to keep it
        pass

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header, f[0].data[tuple(slicing)]


def select_connected_island(filename, coords, outfile=None):
    """
    Select only connected pixels starting from coords on a mask file.
    All other pixels will be set to 0.
    coords: SkyCoord with units
    """
    from scipy.ndimage import label

    if (outfile == None):
        outfile = filename
    
    with pyfits.open(filename) as fits:
        mask = fits[0].data
        w = wcs.WCS(fits[0].header).celestial
        selected_pixel = w.world_to_pixel(SkyCoord(coords[0], coords[1]))

        # Label connected components in the array
        labeled_array, __num_features = label(mask.byteswap().newbyteorder())
        # Get the label of the starting pixel
        #print(labeled_array.shape, selected_pixel)
        selected_label = labeled_array[0,0,int(selected_pixel[1]),int(selected_pixel[0])]
        # Extract only the pixels belonging to the same island as the starting pixel
        mask[:] = 0
        mask[labeled_array == selected_label] = 1

        fits[0].data = mask
        fits.writeto(outfile, overwrite=True)


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
        mask = fits[0].data.astype(bool)
    
    if (inverse): mask = ~(mask)

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


def regrid(image_in, header_from, image_out):
    """
    Regrid 'image_in' to the header of 'header_from' and write it in 'image_out'
    """
    from astropy.io import fits
    from reproject import reproject_interp, reproject_exact
    reproj = reproject_exact

    # get input and header for regridding
    header_rep, data_rep = flatten(fits.open(header_from))
    header_in, data_in = flatten(fits.open(image_in))

    # do the regrid
    logging.info('Regridding %s->%s' % (image_in, image_out))
    data_out, footprint = reproj((data_in, header_in), header_rep, parallel=True)

    # write output
    header_rep =  fits.open(header_from)[0].header
    hdu = fits.PrimaryHDU(header=header_rep, data=[[data_out]])
    hdu.writeto(image_out, overwrite=True)

def add_beam(imagefile, bmaj, bmin, bpa):
    """
    Add/change beam info to fits header
    """
    with pyfits.open(imagefile) as fits:
        header = fits[0].header
        header["BMAJ"] = bmaj
        header["BMIN"] = bmin
        header["BPA"] = bpa
        fits[0].header = header
        fits.writeto(imagefile, overwrite=True)

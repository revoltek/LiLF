import pyregion
from LiLF.lib_img import flatten
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.wcs.utils as utils
from astropy.io import fits as pyfits


def checkregion(pointing, fitsimage, region):
    """
    Check if target region is fully covered in the input image.
    Currently works for box and circles only.
    """

    hdul = fits.open(fitsimage)
    head, datafits = flatten(hdul)
    r = pyregion.open(region)

    if len(r[0].coord_list) == 5:
        region_shape = 1
        racen, deccen, xlen, ylen, angle = r[0].coord_list

        ramax = racen + (xlen)
        ramin = racen - (xlen)

        decmax = deccen + (ylen)
        decmin = deccen - (ylen)

    elif len(r[0].coord_list) == 3:
        region_shape = 2
        racen, deccen, radius = r[0].coord_list

        ramax = racen + (radius)
        ramin = racen - (radius)

        decmax = deccen + (radius)
        decmin = deccen - (radius)


    w = WCS(head)
    coords_cen = SkyCoord(racen*u.deg, deccen*u.deg, frame='fk5')
    coords_max = SkyCoord(ramax*u.deg, decmax*u.deg, frame='fk5')
    coords_min = SkyCoord(ramin*u.deg, decmin*u.deg, frame='fk5')

    rapix_cen, decpix_cen = utils.skycoord_to_pixel(coords_cen, wcs=w)
    rapix_max, decpix_max = utils.skycoord_to_pixel(coords_max, wcs=w)
    rapix_min, decpix_min = utils.skycoord_to_pixel(coords_min, wcs=w)

    im_xcen = (hdul[0].header['CRVAL1'])
    im_ycen = (hdul[0].header['CRVAL2'])
    xshift = (hdul[0].header['CDELT1'])
    yshift = (hdul[0].header['CDELT2'])
    pixnum_x = (hdul[0].header['NAXIS1'])
    pixnum_y = (hdul[0].header['NAXIS2'])

    minra_image = im_xcen - (abs(xshift)*pixnum_x)
    maxra_image = im_xcen + (abs(xshift)*pixnum_x)
    mindec_image = im_ycen - (abs(yshift)*pixnum_y)
    maxdec_image = im_ycen + (abs(yshift)*pixnum_y)


    if minra_image < 0:
        minra_image = 360 + minra_image
    if maxra_image < 0:
        maxra_image = 360 + maxra_image
    if mindec_image < 0:
        mindec_image = 360 + mindec_image
    if maxdec_image < 0:
        maxdec_image = 360 + maxdec_image

    # TODO add beam sensitivity condition: e.g. when sensitivity drops to 50%
    if ramin > minra_image:
        if ramax < maxra_image:
            if decmin > mindec_image:
                if decmax < maxdec_image:
                    checked = pointing

                    return checked

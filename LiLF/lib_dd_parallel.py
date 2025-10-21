import os, sys
import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs as pywcs
import pyregion
import mocpy
import astropy.units as u
from pyregion.parser_helper import Shape
from matplotlib.path import Path
from scipy.ndimage import binary_dilation, generate_binary_structure
from scipy.ndimage.measurements import label, center_of_mass
try:
    from scipy.spatial import Voronoi, voronoi_plot_2d
except:
    logger.error("Load latest scipy with 'use Pythonlibs'")
    sys.exit(1)

from LiLF.lib_log import logger
from LiLF import lib_img


def check_lotss_coverage(center, size):
    """ check if is in LoTSS DR3, this is mostly borrowed from RAPTHOR / D. Rafferty

    Parameters
    ----------
    center: [ra,deg] in degrees
    size: float, square size in degrees

    Returns
    -------
    is_covered: bool,
    """
    ra, dec = center
    logger.debug('Checking LoTSS coverage for the requested centre and radius.')

    moc = mocpy.MOC.from_fits(os.path.dirname(__file__) + '/../models/lotss_dr3_moc.fits')
    covers_centre = moc.contains(ra * u.deg, dec * u.deg)

    # Checking single coordinates, so get rid of the array
    covers_left = moc.contains(ra * u.deg - size * u.deg, dec * u.deg)[0]
    covers_right = moc.contains(ra * u.deg + size * u.deg, dec * u.deg)[0]
    # Ensure dec-size does not exceed -90 deg (south celestial pole)
    dec_bottom = max(dec - size, -90.0)
    covers_bottom = moc.contains(ra * u.deg, dec_bottom * u.deg)[0]
    # Ensure dec+size does not exceed 90 deg (north celestial pole)
    dec_top = min(dec + size, 90.0)
    covers_top = moc.contains(ra * u.deg, dec_top * u.deg)[0]

    fully_covered = False
    if covers_left and covers_right and covers_bottom and covers_top and covers_centre:
        fully_covered = True
    return fully_covered


def closest_distance_between_patches(skymodel):
    """
    Return the name of the two patches which are closes and their distance

    Parameters
    ----------
    skymodel

    Returns
    -------
    closest_name: (str, str) - names of the two closest patches
    closest_distance: float - distance
    """

    closest_patch = np.zeros(len(skymodel.getPatchNames()))
    closest_name = []
    names = skymodel.getPatchNames()
    for i, (name, pos) in enumerate(skymodel.getPatchPositions().items()):
        distances = skymodel.getDistance(*pos, byPatch=True)
        closest_patch[i] = np.sort(distances)[1] # [0] is the patch itself if there is only ONE patch with that distance
        nearby_name = names[distances <= closest_patch[i]] # select all patches with this distance in case multiple have SAME distance
        # Case multiple patches at same distance
        if len(nearby_name) > 1:
            if nearby_name[0] == name:
                nearby_name = nearby_name[1]
            else:
                nearby_name = nearby_name[0]
            if name == nearby_name:  # sanity check, it should not be identical.
                raise ValueError(f'A: patch {name} is identical to closest patch {nearby_name} at distance {closest_patch[i]}!')
        else:
            nearby_name = nearby_name[0]
            if name == nearby_name:  # sanity check, it should not be identical.
                raise ValueError(f'B: patch {name} is identical to closest patch {nearby_name} at distance {closest_patch[i]}!')
        closest_name.append([name,nearby_name])

    name_closest, dist_closest = closest_name[np.argmin(closest_patch)], np.min(closest_patch)
    return name_closest, dist_closest


def merge_nearby_bright_facets(skymodel, max_distance, min_flux, applyBeam=False):
    """
    Merge all bright patches of a skymodel that are within min_distance of another patch

    Parameters
    ----------
    skymodel
    max_distance: max distance to merge
    min_flux: min flux of facets to be considered bright
    applyBeam

    Returns
    -------
    merged skymodel
    """
    skymodel = skymodel.copy()
    skymodel_bright = skymodel.copy()
    skymodel_bright.select(f'I>{min_flux}', aggregate='sum', applyBeam=applyBeam)
    if len(skymodel_bright) > 1:
        # loop over the bright patches as long as the distance between the two closest patches is less than max_distance
        while closest_distance_between_patches(skymodel_bright)[1] < max_distance:
            closest_patches, closest_distance = closest_distance_between_patches(skymodel_bright)
            logger.info(f'Merging nearby bright patches {closest_patches[0]} {closest_patches[1]} (distance={closest_distance*3600:.2f}arcsec')
            skymodel_bright.merge(closest_patches)
            skymodel.merge(closest_patches)
    else:
        logger.warning(f'Only one bright source - nothing to merge.')
    skymodel.setPatchPositions(method='wmean', applyBeam=applyBeam)
    return skymodel


def rename_skymodel_patches(skymodel, applyBeam=False):
    """
    Rename the patches in the input sky model according to flux

    Parameters
    ----------
    skymodel : LSMTool skymodel.SkyModel object
        Input sky model
    applyBeam : bool, intrinsic/apparent
    """
    if not skymodel.hasPatches:
        raise ValueError('Cannot rename patches since the input skymodel is not grouped '
                         'into patches.')
    patch_names = skymodel.getPatchNames()
    patch_fluxes = skymodel.getColValues('I', aggregate='sum', applyBeam=applyBeam)
    patch_pos = skymodel.getPatchPositions()

    old_new_dict = {}
    for i, id in enumerate(np.argsort(patch_fluxes)[::-1]):
        old_new_dict[patch_names[id]] = f'patch_{i:02.0f}'

    patch_col = skymodel.getColValues('Patch')
    for old_name, new_name in old_new_dict.items():
        patch_col[patch_col == old_name] = new_name
        patch_pos[new_name] = patch_pos.pop(old_name)

    skymodel.setColValues('Patch', patch_col)
    skymodel.setPatchPositions(patch_pos)

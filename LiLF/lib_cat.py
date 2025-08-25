# Copyright (C) 2022 - Henrik Edler
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

import os
import numpy as np
import pandas as pd
import copy
import logging
import lsmtool
from astropy.table import Table, Column
import astropy.units as u
from astropy.units import UnitConversionError
from astropy.coordinates import SkyCoord

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return SkyCoord(c_ra,c_dec).separation(SkyCoord(ra,dec)).value
    #return np.sqrt((np.cos(c_dec.to_value('rad'))*(ra-c_ra))**2.0+(dec-c_dec)**2.0)

class Cat():
    """ Wrapper class to include filtering and cross-matching utilities for astropy tables"""
    def __init__(self, cat, catname, log=None, ra=None, dec=None, wcs=None):
        if log is None:
            log = logging.getLogger('lib_catalog.logger')
            log.setLevel('INFO')

        self.log = log
        if isinstance(cat, str):
            cat = Table.read(cat)
        elif isinstance(cat, Table):
            pass
        else:
            raise TypeError(f"cat must be a filename or astropy Table, not {type(cat)}")
        # fix column names and units
        if ra:
            try:
                cat['RA'] = cat[ra].to('degree')
            except UnitConversionError:
                log.warning('No unit in RA input - assume deg.')
                cat['RA'] = cat[ra]*u.deg
        else:
            for racol in ['RA', 'RAJ2000', 'RA_J2000', '_RAJ2000']:
                try:
                    try:
                        cat['RA'] = cat[racol].to('degree')
                    except UnitConversionError:
                        cat['RA'] = cat[racol]*u.deg
                        log.warning('No unit in RA input - assume deg.')
                        break
                except KeyError:
                    pass
        # be sure it's between 0 and 360
        cat['RA'] = cat['RA'] % 360
        
        if dec:
            try:
                cat['DEC'] = cat[dec].to('degree')
            except UnitConversionError:
                print('No unit in DEC input - assume deg.')
                cat['DEC'] = cat[dec]*u.deg
        else:
            for decol in ['DEC', 'DECJ2000', 'DEJ2000', 'DE_J2000', 'DEC_J2000', '_DEJ2000']:
                try:
                    try:
                        cat['DEC'] = cat[decol].to('degree')
                    except UnitConversionError:
                        print('No unit in DEC input - assume deg.')
                        cat['DEC'] = cat[decol]*u.deg
                    break
                except KeyError:
                    pass
        self.cat = cat
        self.wcs = wcs
        self.catname = catname

    def __len__(self):
        return len(self.cat)

    def __str__(self):
        return self.cat.__str__()

    def __iter__(self):
        self.index = 0
        return self.cat

    def __next__(self):
        try:
            nextrow = self.cat[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return nextrow

    def __getitem__(self, key):
        if not isinstance(key, str) and not isinstance(key, int) and hasattr(key, '__len__'): # case1: index with boolean mask - return a copy of object with boolean mask applied to internal cat
            try:
                mask = np.array(key).astype(bool)
                return Cat(self.cat[mask], self.catname)
            except ValueError: pass
        else:# case2: string mask, return column of internal cat
            return self.cat[key]

    def __setitem__(self, key, value):
        self.cat[key] = value

    def copy(self):
        return copy.deepcopy(self)

    def keys(self):
        return self.cat.keys()

    def filter(self, isolation=None, rectangle=None, circle=None, ellipse=None):
        """
        Filter the catalogue on a variety of properties.

        Parameters
        ----------
        isolation: float, optional.
            Source isolation in arcsec.
        rectangle: (4,) array-like.
            Filter for sources in a rectangle defined by [minra, maxra, mindec, maxdec]
        circle: (3,) array-like.
            [RA, DEC, radius] in deg for circle filtering.
        ellipse: (5,) array-like
            [RA, DEC, height, width, angle] of filter ellipse in degree.
        """
        cat = self.cat
        logstr = f'Filter {self.catname}. Initial:{len(cat)}'
        if rectangle: # filter in rectangle
            i = (cat['RA'] > rectangle[0]) & (cat['RA'] < rectangle[1]) & (cat['DEC'] > rectangle[2]) & (cat['DEC'] < rectangle[3])
            cat =  cat[i]
            logstr += f'->rectangle:{len(cat)}'
        if circle: # filter in circle
            r = separation(circle[0], circle[1], cat['RA'], cat['DEC'])
            cat['center_dist'] = r
            cat = cat[cat['center_dist'] < circle[2]]
            logstr += f'->circle:{len(cat)}'
        if ellipse: # filter in ellipse (e.g. for primary beam)
            from regions import EllipseSkyRegion
            if not self.wcs:
                raise ValueError('Need catalog with wcs to apply ellipse filter.')
            sc = SkyCoord(ellipse[0], ellipse[1], unit='degree')
            region_sky = EllipseSkyRegion(center=sc, height=2 * ellipse[2] * u.deg, width=2 * ellipse[3] * u.deg, angle= ellipse[4] * u.deg)
            in_beam = region_sky.contains(SkyCoord(cat['RA'], cat['DEC'], unit='degree'), self.wcs)
            cat = cat[in_beam]
            logstr += f'->ellipse:{len(cat)}'
        if isolation:
            cat['NN_dist']=np.nan
            for row in cat:
                dist=3600.0*separation(row['RA']*u.deg,row['DEC']*u.deg,cat['RA'],cat['DEC'])
                dist.sort()
                row['NN_dist']=dist[1]
            cat=cat[cat['NN_dist']> isolation]
            logstr += f'->isolation:{len(cat)}'
        # log.info(logstr)
        self.log.info(logstr)
        self.cat = cat

    def match(self, cat2, radius, cols=[], unique=True):
        """
        Match this catalogue with another.
        Parameters
        ----------
        cat2: RadioCat object, catalogue to match with.
        radius: float, match radius in arcsec
        cols: list, non-standard columns of matchcat to transfer to new cat.

        Returns
        -------
        n_match: int, number of matches
        """
        cat = self.cat
        label = cat2.catname
        oldv = ['RA', 'DEC'] + cols
        for o in oldv:
            cat.add_column(Column(0, name=f'{label}_{o}', dtype=cat2[o].dtype, unit=cat2[o].unit))
        mcols = [f'{label}_{suff}' for suff in ['separation','dRA','dDEC']]
        for mc in mcols:
            cat.add_column(Column(0., name=mc, unit=u.arcsec))
        cat[f'{label}_match'] = False
        rdeg = (radius / 3600.0)* u.deg
        minra = np.min(cat['RA'] - rdeg)
        maxra = np.max(cat['RA'] + rdeg)
        mindec = np.min(cat['DEC'] - rdeg)
        maxdec = np.max(cat['DEC'] + rdeg)
        # pre-filter tab, which may be all-sky
        cat2 = cat2[(cat2['RA'] > minra) & (cat2['RA'] < maxra) & (cat2['DEC'] > mindec) & (cat2['DEC'] < maxdec)]
        matches = 0
        for r in cat:
            dist = separation(r['RA']*cat['RA'].unit, r['DEC']*cat['DEC'].unit, cat2['RA'], cat2['DEC'])*u.deg
            stab = cat2[dist < radius*u.arcsec]
            if len(stab) > 0:
                df = dist[dist < radius*u.arcsec]
                # got at least one match
                if not unique and len(stab) > 1:
                    print(f"Non-unique match {len(stab)} - using closest match.")
                    stab = cat2[dist == np.min(dist)]
                    df = dist[dist == np.min(dist)]
                elif len(stab) > 1:
                    continue
                matches += 1
                for i in range(len(oldv)):
                    r[f'{label}_{oldv[i]}'] = stab[0][oldv[i]]
                r[label + '_separation'] = df[0].to_value('arcsec')
                r[label + '_match'] = True
                r[label + '_dRA'] =  (np.cos(r['DEC']*np.pi/180) * (r['RA']*u.deg - stab[0]['RA']*u.deg)).to_value('arcsec')
                r[label + '_dDEC'] =  (r['DEC']*u.deg - stab[0]['DEC']*u.deg).to_value('arcsec')
        self.log.info(f'{label}: matched {matches} sources')
        return self.get_matches(label)

    def get_matches(self, labels):
        """
        Return the catalogue where we have common matches between the provided labels
        Parameters
        ----------
        labels: list, labels of the catalogues where we want the matches

        Returns
        -------
        Cat object, catalogue where we have the matches
        """
        if isinstance(labels,str): labels = [labels]
        keys = [label+'_match' for label in labels]
        for i, key in enumerate(keys):
            if i == 0:
                matched = self[key].copy()
            else:
                matched &= self[key]
        self.log.info(f'Found {np.sum(matched)} common matches')
        return Cat(self.cat[matched], self.catname)

    def write(self, path, overwrite=False, format=None):
        """
        Write to a catalogue file.
        Parameters
        ----------
        path: string or path-like
        overwrite: boolean
        """
        self.cat.write(path, overwrite=overwrite, format=format)


class RadioCat(Cat):
    """ Wrapper class to include filtering and cross-matching utilities for astropy tables"""
    def __init__(self, cat, catname, log=None, col_tflux='Total_flux', col_pflux='Peak_flux', col_maj='Maj', ra=None, dec=None, wcs=None, resolution=None, scale=1.0):
        # super class
        Cat.__init__(self, cat, catname, log, ra, dec, wcs)
        cat = self.cat
        # fix column names and units
        if col_tflux:
            cat['Total_flux'] = scale*cat[col_tflux].to('Jy')
            try:
                cat['E_Total_flux'] = scale*cat['e_' + col_tflux].to('Jy')
            except KeyError:
                cat['E_Total_flux'] = scale*cat['E_' + col_tflux].to('Jy')
        if col_pflux:
            cat['Peak_flux'] = scale*cat[col_pflux].to('Jy/beam')
            try:
                cat['E_Peak_flux'] = scale*cat['e_' + col_pflux].to('Jy/beam')
            except KeyError:
                cat['E_Peak_flux'] = scale*cat['E_' + col_pflux].to('Jy/beam')
        if col_maj:
            cat['Maj'] = cat[col_maj].to('degree')
        self.resolution = resolution

    def filter(self, isolation=None, rectangle=None, circle=None, ellipse=None, sigma=5, size=None, peak_to_total=None, minflux=None):
        """
        Filter the catalogue on a variety of properties.

        Parameters
        ----------
        isolation: float, optional.
            Source isolation in arcsec.
        rectangle: (4,) array-like.
            Filter for sources in a rectangle defined by [minra, maxra, mindec, maxdec]
        circle: (3,) array-like.
            [RA, DEC, radius] in deg for circle filtering.
        ellipse: (5,) array-like
            [RA, DEC, height, width, angle] of filter ellipse in degree.
        sigma: float.
            Signal to noise.
        size: float.
            major axis filtering in arcsec.
        peak_to_total: float.
            Minimal peak to total flux factor
        """
        cat = self.cat
        logstr = f'Filter {self.catname}. Initial:{len(cat)}'
        if rectangle: # filter in rectangle
            i = (cat['RA'] > rectangle[0]) & (cat['RA'] < rectangle[1]) & (cat['DEC'] > rectangle[2]) & (cat['DEC'] < rectangle[3])
            cat =  cat[i]
            logstr += f'->rectangle:{len(cat)}'
        if circle: # filter in circle
            r = separation(circle[0]*u.deg, circle[1]*u.deg, cat['RA'], cat['DEC'])
            print(circle)
            cat['center_dist'] = r
            cat = cat[cat['center_dist'] < circle[2]]
            logstr += f'->circle:{len(cat)}'
        if ellipse: # filter in ellipse (e.g. for primary beam)
            from regions import EllipseSkyRegion
            if not self.wcs:
                raise ValueError('Need catalog with wcs to apply ellipse filter.')
            sc = SkyCoord(ellipse[0], ellipse[1], unit='degree')
            region_sky = EllipseSkyRegion(center=sc, height=2 * ellipse[2] * u.deg, width=2 * ellipse[3] * u.deg, angle= ellipse[4] * u.deg)
            in_beam = region_sky.contains(SkyCoord(cat['RA'], cat['DEC'], unit='degree'), self.wcs)
            cat = cat[in_beam]
            logstr += f'->ellipse:{len(cat)}'
        if isolation:
            cat['NN_dist']=np.nan*u.deg
            for row in cat:
                dist=separation(row['RA']*u.deg,row['DEC']*u.deg,cat['RA'],cat['DEC'])
                dist.sort()
                row['NN_dist']=dist[1]
            cat=cat[cat['NN_dist'] > isolation/3600.] # deg and arcsec->deg
            logstr += f'->isolation:{len(cat)}'
        if size:
            cat = cat[cat['Maj'] < size*u.arcsec]
            logstr += f'->size:{len(cat)}'
        if minflux:
            cat = cat[cat['Total_flux'] > minflux*u.Jy]
            logstr += f'->fluxcut:{len(cat)}'
        if sigma:
            cat = cat[cat['Total_flux'] > sigma*cat['E_Total_flux']]
            logstr += f'->sigma:{len(cat)}'
        if peak_to_total:
            cat = cat[cat['Peak_flux'] / sigma*cat['Total_flux'] < peak_to_total]
            logstr += f'->peak_to_total:{len(cat)}'
        self.log.info(logstr)
        self.cat = cat

    def match(self, cat2, radius, cols=[], **kwargs):
        """
        Match this catalogue with another.
        Parameters
        ----------
        cat2: RadioCat object, catalogue to match with.
        radius: float, match radius in arcsec
        cols: list, non-standard columns of matchcat to transfer to new cat.

        Returns
        -------
        n_match: int, number of matches
        """
        if type(cat2) == type(self):
            cols = list(np.unique(['Total_flux', 'E_Total_flux'] + cols))
        return Cat.match(self, cat2, radius, cols=cols, **kwargs)

    def flux_ratio(self, label):
        t = self.get_matches(label)
        ratios = t['Total_flux'] / t[label + '_Total_flux']
        self.log.info(f"{self.catname}/{label}")
        median = np.median(ratios)
        self.log.info(f"median: {median}")
        self.log.info(f"percentiles/median: 16%:{np.percentile(ratios, 16)/median}, 84%:{np.percentile(ratios, 84)/median}")
        return median


def get_LOTSS_DR3_cone_as_skymodel(centre, radius, filename, beamMS=None):
    """Do a cone search on LoTSS DR-3 and return the result as skymodel
    center: [ra,dec] in degree
    radius: radius in degree
    filename: skymodel path
    beamMS: path to beamMS
    """
    with open(os.path.dirname(__file__) + '/../models/lotss_dr3_gaus_110325.skymodel', 'r') as f:
        header = f.readline()
        original_colnames = header.replace('\n', '').split(",")
        colnames = header.replace('\n', '').split(" = ")[-1].split(", ")

    table = pd.read_csv(
        os.path.dirname(__file__) + '/../models/lotss_dr3_gaus_110325.skymodel',
        names=colnames, skiprows=1)

    table = table[table['Dec'] >= centre[1] - radius ]
    table = table[table['Dec'] <= centre[1] + radius]

    phasecentre_sky = SkyCoord(centre[0], centre[1], unit=(u.deg, u.deg), frame='icrs')
    skycoords = SkyCoord(table['Ra'], table['Dec'], unit=(u.deg, u.deg), frame='icrs')
    seperations = skycoords.separation(phasecentre_sky).degree
    table = table[seperations < radius]

    table.to_csv(filename, index=False, header=original_colnames)
    sm = lsmtool.load(filename, beamMS=beamMS)
    sm.setColValues('SpectralIndex', [[-0.7]] * len(sm.getColValues('I')))  # add standard spidx
    sm.setColValues('LogarithmicSI', ['True']*len(sm.getColValues('I'))) # add standard spidx
    return sm
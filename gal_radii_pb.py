from astropy.coordinates import ICRS, Galactic, Distance, Angle
from astropy import units as u
import math as mt
import numpy as np
from astropy.table import Table, Column

# inspired by
#http://idl-moustakas.googlecode.com/svn-history/r560/trunk/impro/hiiregions/im_hiiregion_deproject.pro
def compute_rgc(coord, glx_ctr=ICRS('00h42m44.33s +41d16m07.5s'), glx_PA=Angle('37d42m54s'), glx_incl=Angle('77.5d'), glx_dist=Distance(783,unit=u.kpc),deproject=True):
    '''computes deprojected galactocentric distance for an object at coord,
       relative to galaxy at glx_ctr with PA glx_PA, inclination glx_incl, distance glx_distance
       (defaults are for M31)

       deproject parameter controls whether to correct for galaxy inclination or not

       example usage: [NB - may need to be updated for changes in coord package]
       from astropy.coordinates import ICRS
       region = ICRS('00h42m44.31s +41d16m09.4s')
       gal_radii_pb.compute_rgc(region)
         --> <Distance 0.0240779583299347 kpc>
    '''

    # distance from coord to glx centre
    sky_radius = glx_ctr.separation(coord)
    avg_dec = 0.5*(glx_ctr.dec + coord.dec).radian
    x = (glx_ctr.ra - coord.ra)*mt.cos(avg_dec)
    y = glx_ctr.dec - coord.dec
    # azimuthal angle from coord to glx  -- not completely happy with this
    phi = glx_PA - Angle('90d') + Angle(mt.atan(y.arcsec/x.arcsec),unit=u.rad) 

    # convert to coordinates in rotated frame, where y-axis is galaxy major axis
    # have to convert to arcmin b/c can't do sqrt(x^2+y^2) when x and y are angles
    xp = (sky_radius * mt.cos(phi.radian)).arcmin
    yp = (sky_radius * mt.sin(phi.radian)).arcmin

    # de-project, or not
    if deproject:
        ypp = yp/mt.cos(glx_incl.radian)
        obj_radius = mt.sqrt(xp**2+ypp**2) # in arcmin
    else:
        obj_radius = mt.sqrt(xp**2+yp**2) # in arcmin

    # TODO: options for units of output (might want angle rather than distance)
    obj_dist = Distance(Angle(obj_radius, unit=u.arcmin).radian*glx_dist, unit=glx_dist.unit)

    # don't really need this but compute it just for fun
    if abs(Angle(xp, unit = u.arcmin)) < Angle(1e-5, unit = u.rad): # why 1e-5? no idea.
        obj_phi = Angle(0.0, unit = u.rad)
    else:
        if deproject:
            obj_phi = Angle(mt.atan(ypp/xp), unit = u.rad)
        else:
            obj_phi = phi
    return(obj_dist)


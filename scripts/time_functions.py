# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:42:38 2016

@author: rstreet
"""
from sys import argv
from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime
from pyslalib import slalib
import math

def calc_hjd(date_obs,ra,dec,debug=False):
    """Function to convert a timestamp in %Y-%m-%dT%H:%M:%S UTC format
    to Heliocentric Julian Date for a given RA, Dec of target"""
    
    # Convert RA, Dec to radians:
    c = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    if debug==True: print 'RA '+ra+' -> decimal radians '+str(c.ra.radian)
    if debug==True: print 'Dec '+dec+' -> decimal radians '+str(c.dec.radian)
    
    # Convert the timestamp into a DateTime object:
    if 'T' in date_obs: 
        try:        
            dt = datetime.strptime(date_obs,"%Y-%m-%dT%H:%M:%S")
        except ValueError:
            dt = datetime.strptime(date_obs,"%Y-%m-%dT%H:%M:%S.%f")
    else: 
        try:        
            dt = datetime.strptime(date_obs,"%Y-%m-%d %H:%M:%S")
        except ValueError:
            dt = datetime.strptime(date_obs,"%Y-%m-%d %H:%M:%S.%f")
    
    # Calculate the MJD (UTC) timestamp:
    mjd_utc = datetime_to_mjd_utc(dt)
    if debug==True: print 'MJD_UTC = '+str(mjd_utc)
    
    # Correct the MJD to TT:
    mjd_tt = mjd_utc_to_mjd_tt(mjd_utc)
    if debug==True: print 'MJD_TT = '+str(mjd_tt)
    
    # Calculating MJD of 1st January that year:
    (mjd_jan1,iexec) = slalib.sla_cldj(dt.year,1,1)
    if debug==True: print 'MJD of Jan 1, '+str(dt.year)+' = '+str(mjd_jan1)
    
    # Calculating the MJD difference between the DateObs and Jan 1 of the same year:
    tdiff = mjd_tt - mjd_jan1
    if debug==True: print 'Time difference from Jan 1 - dateobs, '+str(dt.year)+' = '+str(tdiff)
    
    # Calculating the RV and time corrections to the Sun:
    (rv,tcorr) = slalib.sla_ecor(c.ra.radian,c.dec.radian,\
                                    dt.year,int(tdiff),(tdiff-int(tdiff)))
    if debug==True: print 'Time correction to the Sun = '+str(tcorr)
    
    # Calculating the HJD:
    hjd = mjd_tt + tcorr/86400.0 + 2400000.5
    if debug==True: print 'HJD = '+str(hjd)
    
    return hjd

def datetime_to_mjd_utc(d):
    """Function to calculate MJD for a given UTC"""

    (mjd, status) = slalib.sla_cldj(d.year, d.month, d.day)
    if status != 0:
        return None
    (fday, status ) = slalib.sla_dtf2d(d.hour, d.minute, d.second+(d.microsecond/1e6))
    if status != 0:
        return None
    mjd_utc = mjd + fday
    
    return mjd_utc

def mjd_utc_to_mjd_tt(mjd_utc, dbg=False):
    """Converts a MJD in UTC (MJD_UTC) to a MJD in TT (Terrestial Time) which is
    needed for any position/ephemeris-based calculations.UTC
    UTC->TT consists of: UTC->TAI = 10s offset + 24 leapseconds (last one 2009 Jan 1.)
    TAI->TT  = 32.184s fixed offset"""

    tt_utc = slalib.sla_dtt(mjd_utc)
    if dbg: print 'TT-UTC(s)=', tt_utc

    mjd_tt = mjd_utc + (tt_utc/86400.0)
    if dbg: print 'MJD(TT)  =  ', mjd_tt

    return mjd_tt


if __name__ == '__main__':
    
    if len(argv) == 1:
        print 'Call sequence: python calctime.py DateObsStr RA Dec'
        exit()
    else:
        dateobs = argv[1]
        ra = argv[2]
        dec = argv[3]

    hjd = calc_hjd(dateobs,ra,dec,debug=True)
    
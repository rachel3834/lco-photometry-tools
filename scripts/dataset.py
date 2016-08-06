# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:22:41 2016

@author: rstreet
"""

import numpy as np
from os import path
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import time_functions

class Star():
    def __init__(self,idstr,rastr,decstr):
        self.id = idstr
        self.coord = None
        
        if ':' in str(rastr) and ':' in str(decstr):
            self.coord = SkyCoord(rastr,decstr,unit=(u.hourangle,u.deg))
        else:
            self.coord = SkyCoord(ra=rastr*u.deg,dec=decstr*u.deg)

class DataSet():
    def __init__(self):
        self.input_dir = None
        self.output_dir = None
        self.input_frames = None
        self.data_cube = np.zeros(1)
        self.nstars = 0
        self.nframes= 0
        self.star_list = []
    
    def load_star_list(self,params):
        self.star_list = []
        entries = open(params['star_file'],'r').readlines()
        for line in entries:
            if line[0:1] != '#' and len(line.replace('\n','')) > 0:
                (idstr,rastr,decstr) = line.replace('\n','').split()
                self.star_list.append( Star(idstr,rastr,decstr) )
    
    def print_star_list(self):
        for star in self.star_list:
            print(str(star.id)+' '+star.coord.to_string('hmsdms'))
    
    def load_data(self,params):
        self.nstars = len(self.star_list)
        if self.input_dir == None:
            self.input_dir = params['input_dir']
        
        self.input_frames = glob.glob(path.join(self.input_dir, '*.fits'))
        self.nframes = len(self.input_frames)
        
        self.data_cube = np.zeros([self.nstars,self.nframes,7])
        
        for i in range(0,len(self.input_frames)):
            print('Reading data for '+path.basename(self.input_frames[i]))
            table = BanzaiTable(self.input_frames[i])
            print('Extracting data for selected stars...')
            stars_data = extract_star_data(self.star_list,table)
            self.data_cube[:,i,0] = time_functions.calc_hjd(table.date_obs,\
                                            stars_data[:,2],stars_data[:,3])
            for k in range(0,stars_data.shape[1]):
                self.data_cube[:,i,k+1] = stars_data[:,k]
            
class BanzaiTable():
    def __init__(self,file_path=None):
        self.file = None
        self.data_file = None
        self.date_obs = None
        self.exptime = None
        self.obs_lat = None
        self.obs_long = None
        self.obs_height = None
        self.x = None
        self.y = None
        self.coord = None
        self.flux = None
        self.flux_err = None
        
        if file_path != None:
            self.read_banzai_table(file_path)
    
    def read_banzai_table(self,file_path):
        """Method to read the photometry table from a single reduced BANZAI 
        FITS data product"""
        
        hdu_list = fits.open(file_path)
        self.data_file = path.basename(file_path)
        self.date_obs = hdu_list[0].header['DATE-OBS']
        self.exptime = hdu_list[0].header['EXPTIME']
        self.obs_lat = hdu_list[0].header['LATITUDE']
        self.obs_long = hdu_list[0].header['LONGITUD']
        self.obs_height = hdu_list[0].header['HEIGHT']
        
        table = hdu_list[1].data
        data = np.zeros([len(table),6])
        self.x = table.field('X')
        self.y = table.field('Y')
        self.coord = SkyCoord(ra=(table.field('RA')*u.degree), \
                            dec=(table.field('DEC')*u.degree))
        self.flux = table.field('FLUX')
        self.flux_err = table.field('FLUXERR')
        
def extract_star_data(star_list,table):
    """Function to extract the data for a given star from a table array
    star should be a SkyCoord object    
    """
    data = np.zeros([len(star_list),6])
    for j in range(0,len(star_list)):
        sep = table.coord.separation(star_list[j].coord)
        idx = np.where(sep == sep.min())[0]
        data[j,0] = table.x[idx]
        data[j,1] = table.y[idx]
        data[j,2] = table.coord[idx].ra.value
        data[j,3] = table.coord[idx].dec.value
        data[j,4] = table.flux[idx]
        data[j,5] = table.flux_err[idx]
        
    return data
    
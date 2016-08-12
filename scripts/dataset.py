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
from matplotlib import use as useBackend
useBackend('Agg')
from matplotlib import pyplot
import statistics

class Star():
    def __init__(self,idstr,rastr,decstr):
        self.id = idstr
        self.coord = None
        
        if ':' in str(rastr) and ':' in str(decstr):
            self.coord = SkyCoord(rastr,decstr,unit=(u.hourangle,u.deg))
        else:
            self.coord = SkyCoord(ra=rastr*u.deg,dec=decstr*u.deg)

class DataSet():
    def __init__(self,params=None):
        self.input_dir = None
        self.output_dir = None
        self.input_frames = None
        self.valid_frames = []
        self.data_cube = np.zeros(1)
        self.residuals_cube = np.zeros(1)
        self.ensemble_flux = np.zeros(1)
        self.ensemble_fluxerr = np.zeros(1)
        self.nstars = 0
        self.nframes= 0
        self.star_list = []
        self.star_ids = []
        if params != None:
            for key, value in params.items():
                setattr(self,key,value)
    
    def load_star_list(self,params):
        self.star_list = []
        self.star_ids = []
        entries = open(params['star_file'],'r').readlines()
        for line in entries:
            if line[0:1] != '#' and len(line.replace('\n','')) > 0:
                (idstr,rastr,decstr) = line.replace('\n','').split()
                star = Star(idstr,rastr,decstr)
                self.star_list.append( star )
                self.star_ids.append( star.id )
    
    def get_star_idx(self,star):
        j = self.star_ids.index(star.id)
        return j
    
    def print_star_list(self):
        for star in self.star_list:
            print(str(star.id)+' '+star.coord.to_string('hmsdms'))
    
    def load_data(self,params):
        self.nstars = len(self.star_list)
        if self.input_dir == None:
            self.input_dir = params['input_dir']
        
        self.input_frames = glob.glob(path.join(self.input_dir, '*.fits'))
        self.input_frames.sort()
        self.nframes = len(self.input_frames)
        
        self.data_cube = np.zeros([self.nstars,self.nframes,7])
        self.residuals_cube = np.zeros([self.nstars,self.nframes,2])
        self.data_cube.fill(np.NaN)
        self.residuals_cube.fill(np.NaN)
        self.ensemble_flux = np.zeros([self.nframes])
        self.ensemble_flux.fill(np.NaN)
        self.ensemble_fluxerr = np.zeros([self.nframes])
        self.ensemble_fluxerr.fill(np.NaN)
        
        for i in range(0,len(self.input_frames)):
            print('Reading data for '+path.basename(self.input_frames[i]))
            table = BanzaiTable(self.input_frames[i])
            if table.got_data == True:
                self.valid_frames.append(table.data_file)
                stars_data = extract_star_data(self.star_list,table)
                ts = time_functions.calc_hjd(table.date_obs,\
                                            stars_data[:,2],stars_data[:,3])
                self.data_cube[:,i,0] = ts - 2450000.0
                for k in range(0,stars_data.shape[1]):
                    self.data_cube[:,i,k+1] = stars_data[:,k]
    
    def diff_photometry(self):
        """Method to compute differential photometry for the selected stars"""
        
        self.calc_residuals()
        print('Calculated residuals')
        self.plot_lightcurves(self.star_list,suffix='norm')
        print('Plotted normalized lightcurves')
        self.calc_ensemble_lightcurve()
        print('Calculated the ensemble lightcurve')
        self.calc_differential_lightcurve()
        print('\nCompleted differential photometry')
    
    def calc_residuals(self,debug=False):
        """Method to calculate the residual lightcurve for all stars in the
        data cube by dividing by the Median Absolute Deviation from each 
        lightcurve, weighted by the photometric errors"""
        
        for j in range(0,self.nstars):
            flux = self.data_cube[j,:,5]
            fluxerr = self.data_cube[j,:,6]
            (wmean,sigma_wmean) = statistics.calc_weighted_mean(flux,fluxerr)
            idx = ~np.isnan(flux)
            self.residuals_cube[j,idx,0] = flux[idx] / wmean
            self.residuals_cube[j,idx,1] = fluxerr[idx] / wmean
            if debug == True:
                print('Mean flux of star '+self.star_ids[j]+' = '+str(wmean))
                
    def calc_ensemble_lightcurve(self):
        """Method to calculate the ensemble lightcurve for differential 
        photometry"""
        
        residuals = self.residuals_cube[1:,:,0]
        sigmas = self.residuals_cube[1:,:,1]
        for i in range(0,self.nframes):
            (self.ensemble_flux[i],sig) = statistics.calc_weighted_mean(residuals[:,i],sigmas[:,i])
            self.ensemble_fluxerr[i] = statistics.calc_weighted_sigma(residuals[:,i],sigmas[:,i],self.ensemble_flux[i])
        
        
    def calc_differential_lightcurve(self):
        """Method to calculate the differential lightcurves of all stars 
        in the data cube"""
        
        ensemble_sigma_sq = 1.0 / (self.ensemble_fluxerr*self.ensemble_fluxerr)
        for j in range(0,self.nstars):
            self.residuals_cube[j,:,0] = self.residuals_cube[j,:,0] / self.ensemble_flux
            sigmas = self.residuals_cube[j,:,1]
            weights = (1.0 / (sigmas*sigmas)) + ensemble_sigma_sq
            self.residuals_cube[j,:,1] = np.sqrt(1.0/weights)
            
            
    def plot_lightcurves(self,selected_stars,suffix=None):
        """Method to plot lightcurve files of a selected range of stars
        from the star list.  Expects a list of Star objects"""
        
        for star in selected_stars:
            j = self.get_star_idx(star)
            fig = pyplot.figure(1)
            ts = self.data_cube[j,:,0]
            f = self.residuals_cube[j,:,0]
            ferr = self.residuals_cube[j,:,1]
            idx = ~np.isnan(f)
            pyplot.errorbar(ts[idx],f[idx],yerr=ferr[idx],fmt='k.', mfc='k', mec='k',ms=2, capsize=1)
            #pyplot.plot(ts[idx],f[idx],'k.',ms=2)
            #print f[idx],ferr[idx]
            (wmean,sig) = statistics.calc_weighted_mean(f[idx],ferr[idx])
            wsig = statistics.calc_weighted_sigma(f[idx],ferr[idx],wmean)
            pyplot.plot([ts[0],ts[-1]], [wmean,wmean],'r-')
            pyplot.plot([ts[0],ts[-1]], [wmean+wsig,wmean+wsig],'r-.')
            pyplot.plot([ts[0],ts[-1]], [wmean-wsig,wmean-wsig],'r-.')
            pyplot.xlabel('HJD-2450000.0')
            pyplot.ylabel('Residual flux')
            pyplot.title('Lightcurve of '+star.id)
            if suffix != None:
                plt_file = path.join(self.output_dir, 'lightcurve_'+suffix+'_'+str(j)+'.png')
            else:
                plt_file = path.join(self.output_dir, 'lightcurve_'+str(j)+'.png')
            pyplot.savefig(plt_file)
            pyplot.close(1)
    
    def plot_ensemble_lightcurve(self):
        """Method to plot the ensemble lightcurve"""
        
        fig = pyplot.figure(1)
        ts = self.data_cube[0,:,0]
        f = self.ensemble_flux
        ferr = self.ensemble_fluxerr
        idx = ~np.isnan(f)
        pyplot.errorbar(ts[idx],f[idx],yerr=ferr[idx],fmt='k.', mfc='k', mec='k',ms=2, capsize=1)
        pyplot.xlabel('HJD-2450000.0')
        pyplot.ylabel('Residual flux')
        pyplot.title('Ensemble lightcurve')
        plt_file = path.join(self.output_dir, 'ensemble_lightcurve.png')
        pyplot.savefig(plt_file)
        pyplot.close(1)
        
    def output_lightcurves(self,selected_stars):
        
        def flux_to_mag(flux):
            mag = 2.5*np.log(flux)
            return mag
        
        def calc_mag(flux,flux_err):
            mag = flux_to_mag(flux)
            df1 = flux - flux_err
            m1 = flux_to_mag(df1)
            df2 = flux + flux_err
            m2 = flux_to_mag(df2)
            mag_err = (m2-m1)/2.0
            return mag, mag_err
        
        for star in selected_stars:
            j = self.get_star_idx(star)
            out_file = path.join(self.output_dir, 'lightcurve_'+str(j)+'.txt')
            fileobj = open(out_file,'w')
            fileobj.write('# HJD-2450000.0    X[pix]   Y[pix]    Flux    Flux_err  Residual_flux Residual_flux_err  Residual_mag Residual_mag_err\n')
            frames = self.input_frames
            ts = self.data_cube[j,:,0]
            x = self.data_cube[j,:,1]
            y = self.data_cube[j,:,2]
            f = self.data_cube[j,:,5]
            ferr = self.data_cube[j,:,6]
            rf = self.residuals_cube[j,:,0]
            rferr = self.residuals_cube[j,:,1]
            (mag,mag_err) = calc_mag(rf,rferr)
            idx = np.where(f > -99.0 )
            for i in idx[0]:
                fileobj.write(path.basename(frames[i])+'  '+str(ts[i])+'  '+str(x[i])+'  '+\
                str(y[i])+'  '+str(f[i])+'  '+str(ferr[i])+'  '+\
                str(rf[i])+' '+str(rferr[i])+'  '+\
                str(mag[i])+'  '+str(mag_err[i])+'\n')
            fileobj.close()
    
class BanzaiTable():
    def __init__(self,file_path=None):
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
        self.got_data = False
        
        if file_path != None:
            self.read_banzai_table(file_path)
    
    def read_banzai_table(self,file_path):
        """Method to read the photometry table from a single reduced BANZAI 
        FITS data product"""
        
        hdu_list = fits.open(file_path)
        if len(hdu_list) == 3:
            self.data_file = path.basename(file_path)
            self.date_obs = hdu_list[0].header['DATE-OBS']
            self.exptime = float(hdu_list[0].header['EXPTIME'])
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
            self.got_data= True
        else:
            print('ERROR: '+path.basename(file_path)+' has too few FITS extensions')
            self.got_data = False
            
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
    
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:18:38 2016

@author: rstreet
"""

from sys import argv, exit
from os import path
import dataset
import glob

def run_diff_phot():
    """Function to perform difference photometry for selected stars within 
    a given dataset"""

    params = get_params()
    check_sanity(params)

    phot_data = dataset.DataSet(params=params)
    phot_data.load_star_list(params)
    phot_data.print_star_list()
    phot_data.load_data(params)
    
    phot_data.plot_lightcurves(phot_data.star_list)
    phot_data.output_lightcurves(phot_data.star_list)
    
    phot_data.diff_photometry()

# Output data products

def get_params():
    """Function to gather the required commandline arguments for different
    photometry"""
    
    params = {}
    if len(argv)!= 4:
        params['input_dir'] = raw_input('Please enter the input directory path: ')
        params['output_dir'] = raw_input('Please enter the output directory path: ')
        params['star_file'] = raw_input('Please enter the path to the starlist file: ') 
    else:
        params['input_dir'] = argv[1]
        params['output_dir'] = argv[2]
        params['star_file'] = argv[3]

    return params

def check_sanity(params):
    """Function to perform basic checks to see if the reduction can continue"""
    
    for dpath in ['input_dir','output_dir']:
        if path.isdir(params[dpath]) == False:
            print('ERROR: Cannot find directory '+params[dpath])
            exit()
    
    if path.isfile(params['star_file']) == False:
        print('ERROR: Cannot find star file '+params['star_file'])
        exit()


if __name__ == '__main__':
    run_diff_phot()
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:11:59 2016

@author: rstreet
"""

import numpy as np

def calc_mad(data,sigma):
    """Function to calculate the MAD and RMS of a data vector and uncertainties"""
    
    weights = 1.0 / (sigma*sigma)
    wdata = data * weights
    
    m = np.median(wdata) / weights.sum()
    
    mad = np.median( abs(data - m) )    
    
    return mad

if __name__ == '__main__':
    
    data = np.random.normal(12.0, 0.2, 100)
    sigma = np.random.normal(0.2,0.05,100)
    mad = calc_mad(data,sigma)
    print('MAD = '+str(mad))
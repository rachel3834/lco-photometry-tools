# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:11:59 2016

@author: rstreet
"""

import numpy as np

def calc_weighted_mean(data,sigma):
    """Function to calculate the mean of a data vector, weighted by the 
    given uncertainties"""
    
    idx = ~np.isnan(data)
    if True in idx:
        sigma_sq = 1.0/(sigma*sigma)
        wmean = (data[idx]*sigma_sq[idx]).sum() / (sigma_sq[idx].sum())
        sigma_wmean = 1.0 / (sigma_sq[idx].sum())
    else:
        wmean = np.NaN
        sigma_wmean = np.NaN
    return wmean, sigma_wmean

def calc_weighted_sigma(data, sigma, wmean):
    """Function to calculate the standard deviation of a data vector, weighted by the
    uncertainties"""

    idx = ~np.isnan(data)
    if True in idx:
        weights = 1.0/(sigma[idx]*sigma[idx])
        ddata = data[idx] - wmean
        wsigma = np.sqrt((ddata*ddata*weights).sum()/weights.sum())
    else:
        wsigma = np.NaN
    return wsigma
    
def calc_mad(data,sigma):
    """Function to calculate the MAD of a data vector and uncertainties"""
    
    idx = ~np.isnan(data)
    weights = 1.0 / (sigma*sigma)
    wdata = data * weights
    m = np.median(wdata[idx]) / weights[idx].sum()
    mad = np.median( abs(data[idx] - m) )
    return mad

if __name__ == '__main__':
    
    data = np.random.normal(12.0, 0.2, 100)
    sigma = np.random.normal(0.2,0.05,100)
    mad = calc_mad(data,sigma)
    print('MAD = '+str(mad))
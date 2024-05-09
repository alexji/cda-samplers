import numpy as np
from scipy import optimize, stats

def lnprob_emcee(theta, rv, rverr, feh, feherr):
    """ Likelihood and Prior """
    pgal, pbg1, \
    vhel, lsigv, feh0, lsigfeh, \
    vbg1, lsigvbg1, fehbg1, lsigfeh1, \
    vbg2, lsigvbg2, fehbg2, lsigfeh2 = theta
    
    ## The prior is just a bunch of hard cutoffs
    if (pgal > 1) or (pgal < 0) or (pbg1 > 1) or (pbg1 < 0) or \
        (lsigv > 3) or (lsigvbg1 > 3) or (lsigvbg2 > 3) or \
        (lsigv < -1) or (lsigvbg1 < -1) or (lsigvbg2 < -1) or \
        (lsigfeh > 1) or (lsigfeh1 > 1) or (lsigfeh1 > 1) or \
        (lsigfeh < -3) or (lsigfeh1 < -3) or (lsigfeh1 < -3) or \
        (vhel > 150) or (vhel < 50) or (vbg1 > 500) or (vbg1 < 50) or \
        (vbg2 > 50) or (vbg2 < -50):
        return -1e10 # outside of prior, return a tiny number
    
    ## Compute log likelihood in rv
    lgal_vhel = stats.norm.logpdf(rv, loc=vhel, scale=np.sqrt(rverr**2 + (10**lsigv)**2))
    lbg1_vhel = stats.norm.logpdf(rv, loc=vbg1, scale=np.sqrt(rverr**2 + (10**lsigvbg1)**2))
    lbg2_vhel = stats.norm.logpdf(rv, loc=vbg2, scale=np.sqrt(rverr**2 + (10**lsigvbg2)**2))
    
    ## Compute log likelihood in feh
    lgal_feh = stats.norm.logpdf(feh, loc=feh0,    scale=np.sqrt(feherr**2 + (10**lsigfeh)**2))
    lbg1_feh = stats.norm.logpdf(feh, loc=fehbg1, scale=np.sqrt(feherr**2 + (10**lsigfeh1)**2))
    lbg2_feh = stats.norm.logpdf(feh, loc=fehbg2, scale=np.sqrt(feherr**2 + (10**lsigfeh2)**2))
    
    ## Note: If for some reason you have covariances, e.g. for Gaia proper motions, 
    ## you can use stats.multivariate_normal.logpdf
    
    ## Combine the components
    lgal = np.log(pgal) + lgal_vhel + lgal_feh
    lbg1 = np.log(pbg1) + lbg1_vhel + lbg1_feh
    lbg2 = np.log(1-pbg1) + lbg2_vhel + lbg2_feh
    ## np.logaddexp takes the exp, adds them, and re-lns them in a numerically stable way
    lbgtot = np.logaddexp(lbg1, lbg2)
    ltot = np.logaddexp(lgal, np.log(1-pgal) + lbgtot)
    
    return np.sum(ltot)

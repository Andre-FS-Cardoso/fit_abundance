import numpy as np

## Calculation of the distances of each HII region (Boczko 1998, Scarano 2008, Scarano et al. 2008, Giovanelli et al. 1994)

def distances(ra, ra0, dec, dec0, pa, ba, d, re):

    # --- ra = array of right ascension of each HII region in degrees
    # --- ra0 = right ascension of the galaxy center in degrees
    # --- dec = array of declination of each HII region in degrees
    # --- dec0 = declination of the galaxy center in degrees
    # --- pa = position angle in degrees
    # --- ba = ratio between the minor and major axes of the galaxy
    # --- d = distance to the galaxy in Mpc
    # --- re = effective radius of the galaxy in kpc
    
    ra = np.array(ra)*np.pi/180
    dec = np.array(dec)*np.pi/180
    pa = pa*np.pi/180
    dec0 = dec0*np.pi/180
    ra0 = ra0*np.pi/180

    cos_i = np.sqrt((ba**2-0.13**2)/(1-0.13**2))

    r1 = -(ra-ra0)*np.sin(pa)*np.cos(dec) + (dec-dec0)*np.cos(pa)

    r2 = (-(ra-ra0)*np.cos(pa)*np.cos(dec) - (dec-dec0)*np.sin(pa))/cos_i

    d_kpc = d * 1e3 # Mpc p/ kpc

    r = np.sqrt(r1**2 + r2**2) * d_kpc / re
    
    return r

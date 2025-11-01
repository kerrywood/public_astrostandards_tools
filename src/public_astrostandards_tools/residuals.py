from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import scipy.optimize
from . import astro_time
from . import observations
from . import sgp4
from . import sensor
from . import orbit_utils
from . import residuals

# -----------------------------------------------------------------------------------------------------
def shortestAngle( angles : np.array ):
    # compute the RMS based on the residuals (https://stackoverflow.com/questions/1878907/how-can-i-find-the-smallest-difference-between-two-angles-around-a-point#7869457)
    return (angles + 180) % 360 - 180

# -----------------------------------------------------------------------------------------------------
def getUVW( sv_df : pd.DataFrame ):
    '''
    return the UVW vectors from a frame with teme_p and teme_v vectors
    we'll use this when projecting an observation into the computed orbit frame
    (see Figure 3 of the ROTASDll documentation (pg 19 in version 9.5)

    Note: this will annotate the input frame
    '''
    # we'll do this in numpy arrays / vectors for speed (push lists to DataFrame)
    tP = np.vstack( sv_df['teme_p'].values )
    tV = np.vstack( sv_df['teme_v'].values )
    # get the U coordinate
    U_c = tP / np.linalg.norm( tP, axis=1 )[:,np.newaxis]
    sv_df['U_c'] = U_c.tolist()
    # get the W coordinate
    W_c = np.cross( tP, tV, axis=1 )
    W_c = W_c / np.linalg.norm( W_c, axis=1 )[:,np.newaxis]
    sv_df['W_c'] = W_c.tolist()
    # get the V coordinate
    V_c = np.cross( W_c, U_c, axis=1 )
    sv_df['V_c'] = V_c.tolist()
    # return the data 
    return sv_df 

# =====================================================================================================
# ROTAS-like functions
# 
# ROTAS tries to calculate important metrics like TOES and BETA from a single ob, including EO / angles
# only.  We'll try to mimic / improve upon those functions here.  Reference the RotasDll library 
# documentation or the Ford Pilco Astrodynamics reference for background.
# 
# BLUF:
# we need to guess the object position ( UDL EO are angles-only; that obviously doesn't tell us where 
# the object is.  The ROTAS documentation has a rudimentary algorithm to guess the range.  This process
# is iterative: make your first guess for orbit radius based on computed orbit, then iterate until you 
# get close to what you think the right value is.  There are other notes about ASW algorithms as well,
# but they are poorly defined.  We'll offer up some options.
#
# Key steps
#   - identify the most likely range on the EO ob (via plane intersection or other)
#   - identify the point that represents
#   - calculate orbit residuals (TOES, BETA) based on that guess
#   - PROFIT!
# =====================================================================================================

# -----------------------------------------------------------------------------------------------------
# plane intersection : closed-form calculation to find when the ray of the observation (ECI) intersects
# the orbit plane.  Will not work for co-planar (or close to co-planar obs)
# -----------------------------------------------------------------------------------------------------
def plane_intersection( 
                       eph_df : pd.DataFrame,   # hypothesis / TLE data
                       obs_df : pd.DataFrame,   # observations (calculated from the sensor to the hypothesis)
                       sen_df : pd.DataFrame ):

    '''
    All must be time / row aligned:
        eph_df : computed / test ephemeris
        obs_df : observations 
        sen_df : sensor location

    Assumes some geometric diversity; if co-planar, this will fail
    '''
    # test ephemeris / computed location
    temep = np.vstack( eph_df['teme_p'] )
    temev = np.vstack( eph_df['teme_v'] )
    # sensor eci vector
    senp  = np.vstack( sen_df['teme_p'] )
    # observation look vectors
    obslv = np.vstack( obs_df['teme_lv'] )

    # orbit momentum vector
    omom  = np.cross( temep, temev, axis=1 )
    omom  = omom / np.linalg.norm(omom,axis=1)[:,np.newaxis]
   
    # numerator / denominator
    num   = np.sum( -senp * omom, axis=1 )
    den   = np.sum( obslv * omom, axis=1 )

    # TODO: if the denominator is zero.. we are coplanar('ish); need another method
    # badidx = np.abs( den ) < 1e-3  # what is the right limit here?
    ranges =num / den
    return ranges, obslv*ranges[:,np.newaxis]

# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
def UDL_residuals( udl_obs : pd.DataFrame, hypothesis_obs : pd.DataFrame ):
    '''
    udl_obs :   UDL obs that have gone through prepUDLObs and have been rotated.  Should contain
                fields like teme_ra, teme_dec

    hypothesis_obs :    returned from `sensor.compute_looks`, should have fields like `XA_TOPO_DEC` and `XA_TOPO_RA`

    TODO : ROTAS-like fields
    '''
    rv = pd.DataFrame()
    if 'teme_ra' in udl_obs     : rv['ra']    = shortestAngle( udl_obs['teme_ra'] - hypothesis_obs['XA_TOPO_RA'] )
    if 'teme_dec' in udl_obs    : rv['dec']   = shortestAngle( udl_obs['teme_dec'] - hypothesis_obs['XA_TOPO_DEC'] )
    if 'azimuth' in udl_obs     : rv['az']    = shortestAngle( udl_obs['azimuth'] - hypothesis_obs['XA_TOPO_AZ'] )
    if 'elevation' in udl_obs   : rv['el']    = shortestAngle( udl_obs['elevation'] - hypothesis_obs['XA_TOPO_EL'] )
    if 'range' in udl_obs       : rv['range'] =  udl_obs['range'] - hypothesis_obs['XA_TOPO_RANGE'] 
    return rv

# =====================================================================================================
if __name__ == '__main__':
    import os
    import public_astrostandards as PA
    from . import test_helpers
    PA.init_all()

    # load the stored time constants file 
    astro_time.load_time_constants( test_helpers.get_test_time_constants(), PA )

    # TLE data 
    tleF     = os.path.join( test_helpers.get_test_dir(), '19548.tle' )
    #tleF     = os.path.join( test_helpers.get_test_dir(), '40294.tle' )
    with open( tleF ) as F: lines = F.readlines()
    L1 = lines[0].strip()
    L2 = lines[1].strip()

    # load obs 
    obsF      = os.path.join( test_helpers.get_test_dir(), '19548.json.gz' )
    print(obsF)
    #obsF      = os.path.join( test_helpers.get_test_dir(), '40294.json.gz' )
    obs_df    = pd.read_json( obsF ).sort_values(by='obTime')
    # reformat UDL obs (these are our actual TEST obs ; from a sensor)
    obs_df    = observations.prepUDLObs( obs_df, PA )

    # now that those have the correct date fields, peel the dates out of the obs.. we'll need those later
    date_df   = obs_df[ astro_time.DATE_FIELDS ].copy()
    # get a sensor frame (from the OBS again; we're using those ground sensors)
    sensor_df = sensor.prepUDLSensor( obs_df.copy(), PA )

    # --------------------- HYPOTHESIS (from a TLE)
    # get hypothesis obs.. (from a TLE)
    eph_df    = sgp4.propTLE_df( date_df, L1, L2, PA )
    # turn each P,V into osculating elements (for ROTAS comparison)
    eph_df    = orbit_utils.sv_to_osc_df( eph_df, PA )
    # find the U,V,W frame
    eph_df    = getUVW(eph_df)

    # -------------------- start the residual calc
    # now, generate hypothesis looks ( from sensor to eph frame )
    looks_df  = sensor.compute_looks( sensor_df, eph_df, PA )
    ranges, intersect_points = plane_intersection( eph_df, obs_df, sensor_df )
    looks_df['guess_ranges'] = ranges
    looks_df['guess_points'] = intersect_points.tolist()
    residuals_df = residuals.UDL_residuals( obs_df, looks_df )

    # delta_nu -> equation (6) from ROTAS : atan( U_o \dot V_c / U_o \dot U_c ) 
    # (notation is wrong in document; replace U_c with O_i (observed unit vector)
    V_c = np.vstack( eph_df['V_c'] )
    U_c = np.vstack( eph_df['U_c'] )
    O_i = np.vstack( obs_df['teme_lv'] )
    num = np.sum( O_i * V_c , axis=1 )
    den = np.sum( O_i * U_c , axis=1 )
    del_nu = np.arctan2( num, den )
    residuals_df['del_nu'] = del_nu.tolist()

    # approximate true anomaly from hypothesis (computed) and the delta-nu value
    nu_c      = eph_df['XA_KEP_TA']
    nu_o      = nu_c + del_nu

    # now calculate TOES via the SPADOC 4 method (8.3.1.2 in ROTAS documentation, version 9.5)
    ecc_term =  np.sqrt( (1-eph_df['XA_KEP_E']) / (1 + eph_df['XA_KEP_E']) )
    tan_V_c  = np.arctan2( nu_c, 2 )
    tan_V_o  = np.arctan2( nu_o, 2 )
    Ec       = 2 * np.arctan2( ecc_term * tan_V_c , 1 )
    Eo       = 2 * np.arctan2( ecc_term * tan_V_o , 1 )
    Mc       = Ec - eph_df['XA_KEP_E'] * Ec 
    Mo       = Eo - eph_df['XA_KEP_E'] * Eo
    # calculate the actual TOES (eq. 11)
    n        = np.sqrt( 398600.5 / eph_df['XA_KEP_A'] ** 3)
    del_t    = (Mc - Mo) / n
    residuals_df['del_t'] = del_t

#    print('{:20} {:20} '.format( 'Calc/plane range','Obs range'))
#    for A,B in zip(ranges,obs_df['range']):
#        print('{:<20.3f} {:<20.3f} '.format( A,B ) )

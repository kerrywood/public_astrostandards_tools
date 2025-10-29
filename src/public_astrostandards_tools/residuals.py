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

# -----------------------------------------------------------------------------------------------------
# we need to guess the object position ( UDL EO are angles-only; that obviously doesn't tell us where 
# the object is.  The ROTAS documentation has a rudimentary algorithm to guess the range.  This process
# is iterative: make your first guess for orbit radius based on computed orbit, then iterate until you 
# get close to what you think the right value is.
# 
# we'll try a plane intersection approach to find the direction
# -----------------------------------------------------------------------------------------------------
def plane_intersection( 
                       eph_df : pd.DataFrame, 
                       obs_df : pd.DataFrame,
                       sen_df : pd.DataFrame ):

    '''
    eph_df : computed / test ephemeris
    obs_df : observations 
    '''
    # test ephemeris / computed location
    temep = np.vstack( eph_df['teme_p'] )
    temev = np.vstack( eph_df['teme_v'] )
    # look vector from the obs (computed in TEME)
    lookv = np.vstack( obs_df['teme_lv'] )
    # pull the sensor location in TEME
    senp  = np.vstack( sen_df['teme_p'] )
    # compute normals and 
    orbit_normal = np.cross( temep, temev, axis=1 )
    orbit_normal = orbit_normal / np.linalg.norm( orbit_normal, axis=1)[:,np.newaxis]
    t     = - np.sum( orbit_normal * temep) 
    t    /= np.sum( orbit_normal * lookv )
    return t

# -----------------------------------------------------------------------------------------------------
def UDL_residuals( udl_obs : pd.DataFrame, hypothesis_obs : pd.DataFrame ):
    '''
    udl_obs :   UDL obs that have gone through prepUDLObs and have been rotated.  Should contain
                fields like teme_ra, teme_dec

    hypothesis_obs :    returned from `sensor.compute_looks`, should have fields like `XA_TOPO_DEC` and `XA_TOPO_RA`
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

    # TLE
    tleF     = os.path.join( test_helpers.get_test_dir(), '19548.tle' )
    with open( tleF ) as F: lines = F.readlines()
    L1 = lines[0].strip()
    L2 = lines[1].strip()

    # load obs 
    obsF      = os.path.join( test_helpers.get_test_dir(), '19548.json.gz' )
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

    # -------------------- RESIDUAL
    # now, generate hypothesis looks ( from sensor to eph frame )
    looks_df  = sensor.compute_looks( sensor_df, eph_df, PA )
    residuals_df = residuals.UDL_residuals( obs_df, looks_df )

    plane_intersection( eph_df, obs_df, sensor_df )

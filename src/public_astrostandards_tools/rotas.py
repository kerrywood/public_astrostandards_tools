from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import scipy.optimize

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
    print(t)
    return t



# =====================================================================================================
if __name__ == '__main__':
    import astro_time
    import observations
    import sgp4
    import sensor
    import utils
    import time

    import public_astrostandards as PA
    PA.init_all()

    # TLE
    with open('../../tests/27566.tle') as F: lines = F.readlines()
    L1 = lines[0].strip()
    L2 = lines[1].strip()

    # load obs 
    obs_df    = pd.read_json('../../tests/27566.json.gz').sort_values(by='obTime')
    # reformat UDL obs (these are our actual TEST obs ; from a sensor)
    obs_df    = observations.prepUDLObs( obs_df, PA )

    st = time.time()
    # now that those have the correct date fields, peel the dates out of the obs.. we'll need those later
    date_df   = obs_df[ astro_time.DATE_FIELDS ].copy()
    # get a sensor frame (from the OBS again; we're using those ground sensors)
    sensor_df = sensor.prepUDLSensor( obs_df.copy(), PA )

    # --------------------- HYPOTHESIS (from a TLE)
    # get hypothesis obs.. (from a TLE)
    eph_df    = sgp4.propTLE_df( date_df, L1, L2, PA )
    # turn each P,V into osculating elements (for ROTAS comparison)
    eph_df    = utils.sv_to_osc_df( eph_df, PA )
    # find the U,V,W frame
    eph_df    = getUVW(eph_df)

    # -------------------- RESIDUAL
    # now, generate hypothesis looks ( from sensor to eph frame )
    looks_df  = sensor.compute_looks( sensor_df, eph_df, PA )
    residuals_df = observations.UDL_residuals( obs_df, looks_df )

    plane_intersection( eph_df, obs_df, sensor_df )

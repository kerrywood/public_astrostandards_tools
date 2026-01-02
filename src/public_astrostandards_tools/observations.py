from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import scipy.optimize

from . import astro_time

# -----------------------------------------------------------------------------------------------------
def UDL_rotate_TEME_ob( O , harness ):
    '''
    given a single ob (O) with ra / declination fields (J2K), convert to TEME

    assume that date has already been converted to ds50_utc with prepObs or astro_time.convert_time
    '''
    newRA  = (harness.ctypes.c_double)()
    newDec = (harness.ctypes.c_double)()
    harness.AstroFuncDll.RotRADec_EqnxToDate( 
                                        106,
                                        2,
                                        O['ds50_utc'],
                                        O['ra'],
                                        O['declination'],
                                        newRA,
                                        newDec )
    return ( np.float64(newRA),np.float64(newDec) )

# -----------------------------------------------------------------------------------------------------
def ra_dec_to_lv( ra, dec ):
    x = np.cos( np.radians(dec) )* np.cos( np.radians( ra ) )
    y = np.cos( np.radians(dec) ) * np.sin( np.radians( ra ) )
    z = np.sin( np.radians(dec ) )
    lv = np.hstack( ( x.values[:,np.newaxis], y.values[:,np.newaxis], z.values[:,np.newaxis] )  )
    return lv

# -----------------------------------------------------------------------------------------------------
# rotate a dataframe of obs into TEME and then also get a TEME look vector (for solving)
def UDL_rotate_TEME_df( df, harness ):
    ''''
    given a set of UDL obs in a dataframe that have been annotated with astro_time.convert_time,
    rotate all from J2K into TEME
    '''
    tv = df.apply( lambda X : UDL_rotate_TEME_ob( X, harness ) , axis=1 )
    df['teme_ra']  = [ X[0] for X in tv ]
    df['teme_dec'] = [ X[1] for X in tv ]
    # x = np.cos( np.radians(df['teme_dec']) ) * np.cos( np.radians( df['teme_ra'] ) )
    # y = np.cos( np.radians(df['teme_dec']) ) * np.sin( np.radians( df['teme_ra'] ) )
    # z = np.sin( np.radians(df['teme_dec'] ) )
    # lv = np.hstack( ( x.values[:,np.newaxis], y.values[:,np.newaxis], z.values[:,np.newaxis] )  )
    # df['teme_lv'] = lv.tolist()
    df['teme_lv'] = ra_dec_to_lv( df['teme_ra'], df['teme_dec'] ).tolist()
    return df

# -----------------------------------------------------------------------------------------------------
# grab some raw UDL obs; convert the time; sort
def prepUDLObs( o_df, harness ):
    '''
    given raw UDL obs, sort them and convert dates
    then run it through astro_time.convert_times to get a standard frame
    '''
    o_df['obTime_dt'] = pd.to_datetime( o_df['obTime'] )
    o_df = o_df.sort_values(by='obTime_dt')
    t_df = astro_time.convert_times( o_df['obTime_dt'], harness )
    o_df = pd.concat( (t_df.reset_index(drop=True),o_df.reset_index(drop=True) ), axis=1 )
    o_df = UDL_rotate_TEME_df( o_df, harness )
    return o_df 
    
# -----------------------------------------------------------------------------------------------------
def synthetic_to_UDL_like(  
                          sensor_df : pd.DataFrame,
                          obs_df    : pd.DataFrame
                          ) : 
    '''
    this maps outputs from our `compute_looks` function to a UDL-like output
    (useful for testing)
    '''
    sensor_df['teme_ra']  = obs_df['XA_TOPO_RA']   # these are already in TEME, good to go.. no need to convert from J2K
    sensor_df['teme_dec'] = obs_df['XA_TOPO_DEC']
    sensor_df['range']    = obs_df['XA_TOPO_RANGE']
    sensor_df['senlat']   = sensor_df['lat']
    sensor_df['senlon']   = sensor_df['lon']
    sensor_df['senalt']   = sensor_df['height']
    sensor_df['teme_lv']  = ra_dec_to_lv( sensor_df['teme_ra'], sensor_df['teme_dec'] ).tolist()
    return sensor_df
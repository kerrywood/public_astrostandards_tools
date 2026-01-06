import ctypes
import numpy as np
import pandas as pd
from . import coordinates

# -----------------------------------------------------------------------------------------------------
def sun_at_time(  df : pd.DataFrame, # must have the times set
                  INTERFACE ):
    '''
    given a set of dates in the format output by time_helpers.convert_times, output the 
    sun position at those times
    '''
    sun_v  = (ctypes.c_double * 3)()
    sun_m  = ctypes.c_double()
    # the routine gives us a look vector and magnitude
    # compute the sun location at times 
    def getSun( X ):
        INTERFACE.AstroFuncDll.CompSunPos( X, sun_v, sun_m ) 
        return list( (ctypes.c_double * 3)( * (np.array( sun_v ) * sun_m ) ) )
    return [ getSun(X)  for X in df['ds50_et'] ]

# -----------------------------------------------------------------------------------------------------
def moon_at_time(  df : pd.DataFrame, # must have the times set
                  INTERFACE ):
    '''
    given a set of dates in the format output by time_helpers.convert_times, output the 
    moon position at those times

    NOTE: this does not annotate or return a DataFrame; it just returns an array of positions
    '''
    sun_v  = (ctypes.c_double * 3)()
    sun_m  = ctypes.c_double()
    # the routine gives us a look vector and magnitude
    # compute the sun location at times 
    def getMoon( X ):
        INTERFACE.AstroFuncDll.CompMoonPos( X, sun_v, sun_m ) 
        return list( (ctypes.c_double * 3)( * (np.array( sun_v ) * sun_m ) ) )
    return [ getMoon(X)  for X in df['ds50_et'] ]

# -----------------------------------------------------------------------------------------------------
def is_sunlit(  df : pd.DataFrame,
                INTERFACE ):
    '''
    given a DataFrame with `teme_p` and `ds50_et` set, determine if the point is sunlit

    NOTE: this does not annotate or return a DataFrame; it just returns an array of positions
    '''
    def closure( ds50_et, teme ):
        tt = (ctypes.c_double * 3)(* teme )
        return INTERFACE.AstroFuncDll.IsPointSunlit( ds50_et, tt )
    return [ closure(A,B) for A,B in zip( df['ds50_et'], df['teme_p'] ) ]  


# -----------------------------------------------------------------------------------------------------
def compute_looks(     
                   df_sensor : pd.DataFrame,
                   df_target : pd.DataFrame,
                   INTERFACE ):
    '''
    those frames must be time-aligned; that's up to you
    
    each row must have 
        ds50_utc 
        lon
        theta
        teme_p
        lat
        lon
        height

    if you're a ground-based sensor, make sure to use LLH_to_TEME to annotate your frame with eci data
    if you're a space-based sensor, make sure you use TEME_to_LLH to get your lon field

    once you do that, you can run this function and it'll compute looks (of course, some fields might not 
    make sense.. if you compute az/el for a space-based sensor.. that means.. something)

    NOTE: this will return a concat'd version of both DataFrames (make copies if you're worried)
        '''
    # we need a data holder for the output of ECIToTopoComps
    TOPO = INTERFACE.helpers.astrostd_named_fields( INTERFACE.AstroFuncDll, prefix='XA_TOPO_' )

    # check that the dates are aligned
    del_t = np.abs( df_target['ds50_utc'].values - df_sensor['ds50_utc'].values )
    assert np.max( np.abs(del_t) ) < 0.00001
    
    # the ECIToTopoComps call relies on astronomical latitude; if it isn't in the frame.. calculate it
    if 'astrolat' not in df_sensor:
        df_sensor['astrolat'] = coordinates.lat_to_astronomical_lat( df_sensor['lat'] )
    
    # concat the sensor and target dataframes and append suffixes
    tdf = pd.concat( (df_sensor.reset_index(drop=True).add_suffix('_sensor'), 
                      df_target.reset_index(drop=True).add_suffix('_target')), 
                      axis=1 )
    

    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # MANUAL VERSION
    # tar_pos = np.vstack( df_target['teme_p'] )
    # sen_pos = np.vstack( df_sensor['teme_p'] )
    
    # # compute the TEME look vector
    # teme_lv      = tar_pos - sen_pos
    # teme_range   = np.linalg.norm( teme_lv, axis=1 )
    # teme_lv_n    = teme_lv / teme_range[:,np.newaxis]

    # # turn those into RA/DEC
    # ra  = np.arctan( teme_lv_n[:,1], teme_lv_n[:,0] )
    # dec = np.arctan( teme_lv_n[:,2] / np.sqrt( teme_lv_n[:,0] ** 2 + teme_lv_n[:,1] ** 2) )
    # ra  = np.rad2deg( ra )
    # dec = np.rad2deg( dec )
    
    # az = INTERFACE.ctypes.c_double()
    # el = INTERFACE.ctypes.c_double()
    
    # # pre-allocate the storage for az / el calcs.. we'll call in a loop
    # azel = np.zeros( (len(ra), 2) )
    # for i, args in enumerate( zip( df_sensor['ds50_utc'], df_sensor['lat'], df_sensor['lon'], ra, dec) ):
    #     INTERFACE.AstroFuncDll.RaDecToAzElTime( *args, az, el)
    #     azel[i,0] = az.value
    #     azel[i,1] = el.value
    
    # tdf['XA_TOPO_RA']    = ra
    # tdf['XA_TOPO_DEC']   = dec
    # tdf['XA_TOPO_AZ']    = azel[:,0]
    # tdf['XA_TOPO_EL']    = azel[:,1]
    # tdf['XA_TOPO_RANGE'] = teme_range
    
    # return tdf
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # ECITopoComps VERSION
    def calcLooks( R ):
        lst = np.radians( R['lon_sensor'] ) + R['theta_sensor']
        if 'eci_v_target' in R: 
            eci_v_target = (ctypes.c_double * 3)( *R['eci_v_target'] )
        elif 'teme_v_target' in R:
            eci_v_target = (ctypes.c_double * 3)( *R['teme_v_target'] )
        else: 
            eci_v_target = (ctypes.c_double * 3)(0,0,0)

        INTERFACE.AstroFuncDll.ECIToTopoComps( lst,
                                               R['astrolat_sensor'],
                                               (ctypes.c_double * 3) (*R['teme_p_sensor']),
                                               (ctypes.c_double * 3) (*R['teme_p_target']),
                                               eci_v_target,
                                               TOPO.data )

        # INTERFACE.AstroFuncDll.ECIToTopoComps( lst, lat, sen_eci, sun_p, fake_v, SUN_TOPO.data )
        return TOPO.toDict()    

    ans =  pd.DataFrame.from_records( tdf.apply( calcLooks, axis=1 ).values )
    rv = pd.concat( (tdf.reset_index(drop=True),ans.reset_index(drop=True)), axis=1 ) 

    return rv


# -----------------------------------------------------------------------------------------------------
def prepUDLSensor( obs_df : pd.DataFrame, INTERFACE ):
    '''
    UDL observations have some consistent fields; set those up for processing by our standard functions

    NOTE: this is mostly for *ground-based* optics
    EX:
        senlat : sensor latitude (deg)
        senlon : sensor longitude (deg)
        senalt : sensor height (km)
    '''
    sensor_df  = obs_df[['ds50_utc','senlat','senlon','senalt','theta']].copy()
    sensor_df  = sensor_df.rename( columns = {'senlat' : 'lat','senlon' : 'lon', 'senalt' : 'height' } )
    sensor_df  = coordinates.LLH_to_TEME( sensor_df, INTERFACE )
    return sensor_df

# -----------------------------------------------------------------------------------------------------
def setup_ground_site( dates_df : pd.DataFrame,
                       lat : float,
                       lon : float,
                       alt : float,
                       INTERFACE ):
    # for now, modify the passed in dataframe
    sensor_f = dates_df
    # test the LLH_to_TEME function
    sensor_f['lat']     = lat
    sensor_f['lon']     = lon
    sensor_f['height']  = alt
    # we need ECI to feed later routines
    sensor_f = coordinates.LLH_to_TEME( sensor_f, INTERFACE )
    return sensor_f
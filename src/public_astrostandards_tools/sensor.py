import sys
import ctypes
import numpy as np
import pandas as pd
# from astrostandards.utils import helpers

# -----------------------------------------------------------------------------------------------------
def llh_to_eci( df : list[ float ],
                INTERFACE ) :
    '''
    given a lat / lon / alt tuple and a set of astrostandard epoch'd dates,
    give back the ECI position (TEME)
    '''
    sen_eci = (ctypes.c_double * 3)()
    def getECI( R ):
        llh = (ctypes.c_double * 3)( R['lat'], R['lon'], R['height'] )
        INTERFACE.AstroFuncDll.LLHToXYZTime( R['ds50_utc'], llh, sen_eci )
        return list( sen_eci )
    df['teme_p'] =  df.apply( getECI, axis=1 )
    return df

# -----------------------------------------------------------------------------------------------------
def llh_to_efg( df : list[ float ],
                INTERFACE) :
    '''
    given a lat / lon / height tuple ,
    give back the EFG position (ECEF)
    '''
    sen_efg = (ctypes.c_double * 3)()
    def getEFG( R ):
        llh = (ctypes.c_double * 3)(R['lat'], R['lon'], R['height'])
        INTERFACE.AstroFuncDll.LLHToEFGPos(llh, sen_efg)
        return list( sen_efg )
    df['efg_p'] = df.apply( getEFG, axis=1 )
    return df

# -----------------------------------------------------------------------------------------------------
def eci_to_llh( df : list[ float ],
                INTERFACE ) :
    '''
    given a dataframe with columns `teme_p` and `ds50_utc`, covert the 
    eci coordinates to llh
    '''
    llh  = (ctypes.c_double * 3)()
    
    def getLLH( R ):
        eci = (ctypes.c_double * 3)( *R['teme_p'] )

        INTERFACE.AstroFuncDll.XYZToLLHTime( R['ds50_utc'], eci, llh )

        return list( llh )
    
    tv = df.apply( getLLH, axis=1 )
    df['lat']    = [ T[0] for T in tv ]
    df['lon']    = [ T[1] for T in tv ]
    df['height'] = [ T[2] for T in tv ]
    return df


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
    sun_p = (ctypes.c_double * 3)( * (np.array( sun_v ) * sun_m ) )
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
    '''
    sun_v  = (ctypes.c_double * 3)()
    sun_m  = ctypes.c_double()
    # the routine gives us a look vector and magnitude
    sun_p = (ctypes.c_double * 3)( * (np.array( sun_v ) * sun_m ) )
    # compute the sun location at times 
    def getMoon( X ):
        INTERFACE.AstroFuncDll.CompMoonPos( X, sun_v, sun_m ) 
        return list( (ctypes.c_double * 3)( * (np.array( sun_v ) * sun_m ) ) )
    return [ getMoon(X)  for X in df['ds50_et'] ]

# -----------------------------------------------------------------------------------------------------
def is_sunlit(  df : pd.DataFrame,
                INTERFACE ):
    teme_p = (ctypes.c_double * 3)()
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
        '''
    # we need a data holder for the output of ECIToTopoComps
    TOPO = INTERFACE.helpers.astrostd_named_fields( INTERFACE.AstroFuncDll, prefix='XA_TOPO_' )

    # check that the dates are aligned
    del_t = np.abs( df_target['ds50_utc'].values - df_sensor['ds50_utc'].values )
    assert np.max(del_t) < 0.00001
    
    # concat the sensor and target dataframes and append suffixes
    tdf = pd.concat( (df_sensor.reset_index(drop=True).add_suffix('_sensor'), 
                      df_target.reset_index(drop=True).add_suffix('_target')), 
                    axis=1 )
    
    def calcLooks( R ):
        lst = np.radians( R['lon_sensor'] ) + R['theta_sensor']
        if 'eci_v_target' in R: 
            eci_v_target = (ctypes.c_double * 3)( *R['eci_v_target'] )
        elif 'teme_v_target' in R:
            eci_v_target = (ctypes.c_double * 3)( *R['teme_v_target'] )
        else: 
            eci_v_target = (ctypes.c_double * 3)(0,0,0)

        INTERFACE.AstroFuncDll.ECIToTopoComps( lst,
                                               R['lat_sensor'],
                                               (ctypes.c_double * 3) (*R['teme_p_sensor']),
                                               (ctypes.c_double * 3) (*R['teme_p_target']),
                                               eci_v_target,
                                               TOPO.data )

        # INTERFACE.AstroFuncDll.ECIToTopoComps( lst, lat, sen_eci, sun_p, fake_v, SUN_TOPO.data )
        return TOPO.toDict()    

    ans =  pd.DataFrame.from_records( tdf.apply( calcLooks, axis=1 ).values )
    rv = pd.concat( (tdf.reset_index(drop=True),ans.reset_index(drop=True)), axis=1 ) 
    return rv


    

# =====================================================================================================
if __name__ == "__main__":
    from datetime import datetime,timedelta,timezone
    import public_astrostandards as harness

    import astro_time 
    import sgp4 

    # init all the Dll's
    harness.init_all()

    # use the TimeFunc to load the time parameters file (need to upate this periodically)
    astro_time.load_time_constants( '/opt/astrostandards/reduced_time_constants.dat', harness )

    # generate some test data
    dates = pd.date_range( '2025-10-01', '2025-10-15',  freq='1min' )
    #dates = pd.date_range( datetime.now( timezone.utc ), datetime.now( timezone.utc ) + timedelta(days=14), freq='1min' )

    # use the astro_time to initialize the dataframe with times
    # note that the sensor and target dataframes must be time aligned
    dates_f = astro_time.convert_times( dates, harness )

    # -----------------------------------------------------------------------------------------------------
    # TEST case 1 : look vectors to sun; find out when sun is down
    # make a sensor frame
    sensor_f = dates_f.copy()
    
    # make a target frame
    target_f = dates_f.copy()

    # test the llh_to_eci function
    sensor_f['lat']     = 38.83 
    sensor_f['lon']     = -104.82
    sensor_f['height']  = 1.832
    # we need ECI to feed later routines
    sensor_f = llh_to_eci( sensor_f, harness )

    # annotate that dataframe wih the sun position
    target_f['teme_p']  = sun_at_time( dates_f, harness ) 
    
    # compute looks to the sun
    looks_f = compute_looks( sensor_f, target_f, harness )

    # find those times when the sun is down; NOTE we're indexing subsets of the DATE and SENSOR frame
    dates_sundown_f  = dates_f[ looks_f['XA_TOPO_EL'] < -4 ].copy()
    sensor_sundown_f = sensor_f[ looks_f['XA_TOPO_EL'] < -4 ].copy()

    # -----------------------------------------------------------------------------------------------------
    # TEST case 2 : now that we know when sun is down, can we find those times when ISS is visible?
    # for now, ignore if it is solar illuminated
    # re-use the dates  
    L1 = '1 25544U 98067A   25274.49975208  .00018288  00000-0  33242-3 0  9998'
    L2 = '2 25544  51.6325 140.1428 0001055 183.8834 176.2147 15.49589290531650'
    sensor_f['lat']     = 38.83 
    sensor_f['lon']     = -104.82
    sensor_f['height']  = 1.832

    # note that we're only checking at those times that the sun is down
    # this limits the propagation calls (more efficient)
    target_f = sgp4.propTLE_df( dates_sundown_f, L1, L2, harness ) 
    # check if the target is sunlit...
    target_f['is_sunlit'] = is_sunlit( target_f, harness )
    looks_f = compute_looks( sensor_sundown_f, target_f, harness )

    # note that here, _target is appended to a field we pushed into compute_looks
    # so we look for _is_sunlit_target
    good = looks_f[ (looks_f['XA_TOPO_EL'] > 5) * (looks_f['is_sunlit_target'] == 1 ) ] 
    print( good[['datetime_sensor','XA_TOPO_EL','XA_TOPO_AZ','is_sunlit_target']] )

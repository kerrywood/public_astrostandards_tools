from datetime import datetime
import numpy as np
import pandas as pd
from . import astro_time
from . import coordinates

# -----------------------------------------------------------------------------------------------------
def UDL_rotate_TEME_ob( udlob , harness ):
    '''
    given a single ob (O) with ra / declination fields (J2K), convert to TEME

    assume that date has already been converted to ds50_utc with prepObs or astro_time.convert_time
    '''
    newRA  = (harness.ctypes.c_double)()
    newDec = (harness.ctypes.c_double)()
    harness.AstroFuncDll.RotRADec_EqnxToDate( 
                                        106,
                                        2,
                                        udlob['ds50_utc'],
                                        udlob['ra'],
                                        udlob['declination'],
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

# -----------------------------------------------------------------------------------------------------
def val_mapper( val, digits=5 ):
    modu = 10 ** digits
    try: 
        return int( val ) % modu
    except:
        pass
    try:  
        return abs(hash(s)) % mod   
    except: pass
    return 99999

# -----------------------------------------------------------------------------------------------------
def satNo( ob : dict ):
    if 'satNo' in ob :
        return val_mapper( ob['satNo'] )
    if 'origObjectId' in ob:
        return val_mapper( ob['origObjectId'] )
    if 'idOnOrbit' in ob:
        return val_mapper( ob['idOnOrbit'] )
    return 999

# -----------------------------------------------------------------------------------------------------
def idSensor( ob : dict ):
    if 'idSensor' in ob:
        return val_mapper( ob['idSensor'] )
    return 999

# -----------------------------------------------------------------------------------------------------
# convert UDL obs to B3 using the astrostandards; assume you pass in a helper to avoid rebuilding it
def UDLEOObtoB3Type9( ob : dict, OBSHELPER, harness ):
    # these obs should already have been converted to EFG
    assert 'efg_p' in ob
    efgx,efgy,efgz = ob['efg_p']
    # get the satNo we should use (modify the function above to change)
    satno = satNo( ob )
    OBSHELPER.clear()
    OBSHELPER['XA_OBS_SECCLASS']  = 1 
    OBSHELPER['XA_OBS_SATNUM']    = satno
    OBSHELPER['XA_OBS_SITETAG']   = satno
    OBSHELPER['XA_OBS_SPADOCTAG'] = satno
    OBSHELPER['XA_OBS_SENNUM']    = ob['fake_sensor_number'] if 'fake_sensor_number' in ob else idSensor( ob )
    OBSHELPER['XA_OBS_DS50UTC']   = ob['ds50_utc']
    OBSHELPER['XA_OBS_ELORDEC']   = ob['declination']
    OBSHELPER['XA_OBS_AZORRA']    = ob['ra']
    OBSHELPER['XA_OBS_POSX']      = efgx
    OBSHELPER['XA_OBS_POSY']      = efgy
    OBSHELPER['XA_OBS_POSZ']      = efgz
    OBSHELPER['XA_OBS_OBSTYPE']   = 9 
    OBSHELPER['XA_OBS_TRACKIND']  = ob['track_indicator'] if 'track_indicator' in ob else 3
    OBSHELPER['XA_OBS_YROFEQNX']  = 2 # J2K equinox
    # --------------- inject into the astrostandards
    ob['asObId']                  = harness.ObsDll.ObsAddFrArray( OBSHELPER.getData() )
    
    # --------------- check for failure, print
    if ob['asObId'] < 0 : 
        print('----------------------------------------- ERROR -----------------------------------------')
        print( json.dumps( ob, indent=4, default=str ) ) 
        print( OBSHELPER.toDict() )
    
    # --------------- extract and return the B3
    b3str = harness.Cstr('',512)
    if harness.ObsDll.ObsGetB3Card( ob['asObId'], b3str ) != 0: 
        return 'ERR'
    return b3str.value.decode('utf-8').rstrip()

# -----------------------------------------------------------------------------------------------------
# convert UDL obs to B3 using the astrostandards; assume you pass in a helper to avoid rebuilding it
def UDLEOObstoB3Type9( obs_df, harness ):
    # setup times and convert coordinates (though we won't use that)
    obs_df = prepUDLObs( obs_df, harness )
    # copy fields (UDL calls them senlat, senlon, senalt.. we need lat, lon, height)
    # we could rename them, but this preseves data all the way through at the expense of duplication
    obs_df['lat']    = obs_df['senlat']
    obs_df['lon']    = obs_df['senlon']
    obs_df['height'] = obs_df['senalt']
    # convert to EFG (for type9)
    obs_df = coordinates.LLH_to_EFG( obs_df, harness )
    OBSHELPER = harness.helpers.astrostd_named_fields( harness.ObsDll, prefix='XA_OBS_' )
    b3 = obs_df.apply( lambda udlob: UDLEOObtoB3Type9(udlob,OBSHELPER,harness), axis=1 )
    obs_df['B3'] = b3.tolist()
    return obs_df
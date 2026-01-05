import ctypes
from datetime import datetime, timezone
import numpy as np
import pandas as pd

'''
Key dataframe names (semi-canonical):
    teme_p      : TEME position vector
    teme_v      : TEME velocity vector
    j2k_p       : J2K position vector
    j2k_v       : J2K velocity vector
    efg_p       : EFG position vector
    efg_v       : EFG velocity vector
'''

# -----------------------------------------------------------------------------------------------------
def TEME_to_J2K( teme: pd.DataFrame , harness ):
    '''
    assume j2k has the following columns:
        ds50_utc   (from astro_time)
        ds50_tai   (from astro_time)
        teme_p     (J2K position vector)
        teme_v     (J2K velocity vector)

    and return...
        ds50_utc
        ds50_tai
        j2k_p
        j2k_v
    '''
    j2k_p = (harness.ctypes.c_double * 3)()
    j2k_v = (harness.ctypes.c_double * 3)()

    def processLine( L ):
        teme_p      = (harness.ctypes.c_double * 3)( *L['teme_p'] )
        teme_v      = (harness.ctypes.c_double * 3)( *L['teme_v'] )
        ds50_tai    = harness.TimeFuncDll.UTCToTAI( L['ds50_utc'] )
        harness.AstroFuncDll.RotDateToJ2K( 1, 106, ds50_tai, teme_p, teme_v, j2k_p, j2k_v )
        return [ list(j2k_p), list(j2k_v) ]

    tv = teme.apply( processLine, axis=1 ).values
    teme['j2k_p' ] = [ X[0] for X in tv ]
    teme['j2k_v' ] = [ X[1] for X in tv ]
    return teme 
        

# -----------------------------------------------------------------------------------------------------
def J2K_to_TEME( j2k : pd.DataFrame , harness ):
    '''
    assume j2k has the following columns:
        ds50_utc   (from astro_time)
        ds50_tai   (from astro_time)
        j2k_p      (J2K position vector)
        j2k_v      (J2K velocity vector)

    and return...
        ds50_utc
        ds50_tai
        teme_p
        teme_v
    '''
    teme_p = (harness.ctypes.c_double * 3)()
    teme_v = (harness.ctypes.c_double * 3)()

    def processLine( L ):
        j2k_p = (harness.ctypes.c_double * 3)( *L['j2k_p'] )
        j2k_v = (harness.ctypes.c_double * 3)( *L['j2k_v'] )
        ds50_tai    = harness.TimeFuncDll.UTCToTAI( L['ds50_utc'] )
        harness.AstroFuncDll.RotJ2KToDate( 1, 106, ds50_tai, j2k_p, j2k_v, teme_p, teme_v )
        return [ list(teme_p), list(teme_v) ]

    tv = j2k.apply( processLine, axis=1 ).values
    j2k['teme_p' ] = [ X[0] for X in tv ]
    j2k['teme_v' ] = [ X[1] for X in tv ]
    return j2k

# -----------------------------------------------------------------------------------------------------
def LLH_to_TEME( df : list[ float ],
                INTERFACE ) :
    '''
    given a lat / lon / alt tuple and a set of astrostandard epoch'd dates,
    give back the ECI position (TEME)

    df must have columns 'lat', 'lon', 'height', and 'ds50_utc'
    '''
    sen_eci = (ctypes.c_double * 3)()
    def getECI( R ):
        llh = (ctypes.c_double * 3)( R['lat'], R['lon'], R['height'] )
        INTERFACE.AstroFuncDll.LLHToXYZTime( R['ds50_utc'], llh, sen_eci )
        return list( sen_eci )
    df['teme_p'] =  df.apply( getECI, axis=1 )
    return df

# -----------------------------------------------------------------------------------------------------
def LLH_to_EFG( df : list[ float ],
                INTERFACE) :
    '''
    given a lat / lon / height tuple ,
    give back the EFG position (ECEF)

    df must have columns 'lat', 'lon', 'height'
    '''
    sen_efg = (ctypes.c_double * 3)()
    def getEFG( R ):
        llh = (ctypes.c_double * 3)(R['lat'], R['lon'], R['height'])
        INTERFACE.AstroFuncDll.LLHToEFGPos(llh, sen_efg)
        return list( sen_efg )
    df['efg_p'] = df.apply( getEFG, axis=1 )
    return df

# -----------------------------------------------------------------------------------------------------
def TEME_to_LLH( df : list[ float ],
                 INTERFACE ) :
    '''
    given a dataframe with columns `teme_p` and `ds50_utc`, covert the 
    eci coordinates to llh

    df must have 'teme_p' and 'ds50_utc'
    '''
    llh  = (ctypes.c_double * 3)()
    
    def convert( R ):
        eci = (ctypes.c_double * 3)( *R['teme_p'] )
        INTERFACE.AstroFuncDll.XYZToLLHTime( R['ds50_utc'], eci, llh )
        return list( llh )
    
    tv = df.apply( convert, axis=1 )
    df['lat']    = [ T[0] for T in tv ]
    df['lon']    = [ T[1] for T in tv ]
    df['height'] = [ T[2] for T in tv ]
    return df

# -----------------------------------------------------------------------------------------------------
def TEME_to_EFG( df : list[ float ],
                 INTERFACE ) :
    '''
    given a dataframe with columns `teme_p` and `ds50_utc`, covert the 
    coordinates to EFG
    '''
    efg_p = (ctypes.c_double * 3)()
    efg_v = (ctypes.c_double * 3)()
    
    def convert( R ):
        teme_p = (ctypes.c_double * 3)( *R['teme_p'] )
        teme_v = (ctypes.c_double * 3)( *R['teme_v'] )
        INTERFACE.AstroFuncDll.ECIToEFGTime( R['ds50_utc'], teme_p, teme_v, efg_p, efg_v )
        return ( list(efg_p), list(efg_v) )
    
    tv = df.apply( convert, axis=1 )
    df['efg_p'] = [ T[0] for T in tv ]
    df['efg_v'] = [ T[1] for T in tv ]
    return df

# -----------------------------------------------------------------------------------------------------
def EFG_to_TEME( df : list[ float ],
                 INTERFACE ) :
    '''
    given a dataframe with columns `efg_p` and `efg_v` and `ds50_utc` convert the
    coordinates to teme
    '''
    teme_p = (ctypes.c_double * 3)()
    teme_v = (ctypes.c_double * 3)()
    
    def convert( R ):
        efg_p = (ctypes.c_double * 3)( *R['efg_p'] )
        efg_v = (ctypes.c_double * 3)( *R['efg_v'] )
        INTERFACE.AstroFuncDll.EFGToECITime( R['ds50_utc'], efg_p, efg_v, teme_p, teme_v )
        return ( list(teme_p), list(teme_v) )
    
    tv = df.apply( convert, axis=1 )
    df['teme_p'] = [ T[0] for T in tv ]
    df['teme_v'] = [ T[1] for T in tv ]
    return df

# -----------------------------------------------------------------------------------------------------
def test( ) :
    import public_astrostandards as PA
    from . import astro_time
    from . import sgp4
    from . import utils

    print('*'*100)
    print('Will convert random coordinates from one frame to another and back.. then quantify error')
    print('*'*100)

    # init the astrostandards (this starts logging)
    PA.init_all()
    # load the time constants
    astro_time.load_time_constants(  utils.get_test_time_constants(), PA )
    # test TLE lines
    L1 = '1 25544U 98067A   25301.52216109  .00016210  00000-0  29595-3 0  9990'
    L2 = '2 25544  51.6346   6.3420 0004740 349.4592  10.6297 15.49553329535849'
    # generate some test data
    DATES      = pd.date_range('2025-11-1', '2025-12-1', freq='5min')
    dates_f    = astro_time.convert_times( DATES, PA )
    original_f = sgp4.propTLE_df( dates_f, L1, L2, PA )

    def findErr( A, B ):
        err = np.vstack(A) - np.vstack(B)
        return np.sum( np.abs(err) )

    # convert to J2K
    to_j2k     = TEME_to_J2K( original_f.copy(), PA ) 
    # now convert back; they should be the same 
    from_j2k   = J2K_to_TEME( to_j2k, PA )

    N = original_f.shape[0]

    err = findErr( original_f['teme_p'].values, from_j2k['teme_p'].values ) 
    print('Error from position TEME->J2K->TEME conversion: {} over {} points'.format( err, N ) )
    err = findErr( original_f['teme_v'].values, from_j2k['teme_v'].values ) 
    print('Error from velocity TEME->J2K->TEME conversion: {} over {} points'.format( err, N ) )

    # convert to J2K
    to_efg   = TEME_to_EFG( original_f.copy(), PA ) 
    # now convert back; they should be the same 
    from_efg = EFG_to_TEME( to_efg, PA )

    err = findErr( original_f['teme_p'].values, from_efg['teme_p'].values ) 
    print('Error from position TEME->EFG->TEME conversion: {} over {} points'.format( err, N ) )
    err = findErr( original_f['teme_v'].values, from_efg['teme_v'].values ) 
    print('Error from velocity TEME->EFG->TEME conversion: {} over {} points'.format( err, N ) )



# =====================================================================================================
if __name__ == '__main__':
    test()

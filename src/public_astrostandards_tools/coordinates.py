from datetime import datetime, timedelta, timezone
import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------------------------------
def TEME_to_J2K( teme: pd.DataFrame , harness):
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
        teme_p = (harness.ctypes.c_double * 3)( *L['teme_p'] )
        teme_v = (harness.ctypes.c_double * 3)( *L['teme_v'] )
        ds50_tai    = PA.TimeFuncDll.UTCToTAI( L['ds50_utc'] )
        harness.AstroFuncDll.RotDateToJ2K( 1, 106, ds50_tai, teme_p, teme_v, j2k_p, j2k_v )
        return [ list(j2k_p), list(j2k_v) ]

    tv = teme.apply( processLine, axis=1 ).values
    teme['j2k_p' ] = [ X[0] for X in tv ]
    teme['j2k_v' ] = [ X[1] for X in tv ]
    return teme 
        

# -----------------------------------------------------------------------------------------------------
def J2K_to_TEME( j2k : pd.DataFrame , harness):
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
    j2k_v  = (harness.ctypes.c_double * 3)()
    teme_p = (harness.ctypes.c_double * 3)()
    teme_v = (harness.ctypes.c_double * 3)()

    def processLine( L ):
        j2k_p = (harness.ctypes.c_double * 3)( *L['j2k_p'] )
        j2k_v = (harness.ctypes.c_double * 3)( *L['j2k_v'] )
        ds50_tai    = PA.TimeFuncDll.UTCToTAI( L['ds50_utc'] )
        harness.AstroFuncDll.RotJ2KToDate( 1, 106, ds50_tai, j2k_p, j2k_v, teme_p, teme_v )
        return [ list(teme_p), list(teme_v) ]

    tv = j2k.apply( processLine, axis=1 ).values
    j2k['teme_p' ] = [ X[0] for X in tv ]
    j2k['teme_v' ] = [ X[1] for X in tv ]
    return j2k
        
# -----------------------------------------------------------------------------------------------------
def toCCSDS( df : pd.DataFrame ):
    '''
    assume that the input dataframe has
        datetime :   (datetime or datetime-like)
        j2k_p    :   J2K position vector
        j2k_v    :   J2K velocity vector

        generate an output that is CCSDS-like
    '''
    odf = pd.DataFrame()
    odf['datetime']  = df['datetime'].apply( lambda X: X.strftime('%Y-%m-%dT%H:%M:%S.%f') )
    odf['j2k_x']     = df['j2k_p'].apply( lambda X: X[0] )
    odf['j2k_y']     = df['j2k_p'].apply( lambda X: X[1] )
    odf['j2k_z']     = df['j2k_p'].apply( lambda X: X[2] )
    odf['j2k_dx']    = df['j2k_v'].apply( lambda X: X[0] )
    odf['j2k_dy']    = df['j2k_v'].apply( lambda X: X[1] )
    odf['j2k_dz']    = df['j2k_v'].apply( lambda X: X[2] )
    return odf.to_csv( index=None )

# -----------------------------------------------------------------------------------------------------
def fromCCSDS( lines ): 
    '''
    lines : file data
    '''
    def parseLine(L):
        try:
            flds = L.strip().split(',')
            assert len(flds) == 7 
            dt   = datetime.strptime(flds[0], '%Y-%m-%dT%H:%M:%S.%f' )
            return (dt, [ float(X) for X in flds[1:4] ], [float(X) for X in flds[4:] ] )
        except Exception as e: 
            return None
    prsd  = [ parseLine(X) for X in lines ]
    lines = list( filter( lambda X: X is not None, prsd ) )
    tv    = pd.DataFrame( lines, columns=['datetime','j2k_p','j2k_v'] )
    return tv

# =====================================================================================================
if __name__ == '__main__':
    import public_astrostandards as PA
    from . import astro_time
    from . import sgp4
    from . import test_helpers

    PA.init_all()

    astro_time.load_time_constants(  test_helpers.get_test_time_constants(), PA )

    L1 = '1 25544U 98067A   25301.52216109  .00016210  00000-0  29595-3 0  9990'
    L2 = '2 25544  51.6346   6.3420 0004740 349.4592  10.6297 15.49553329535849'

    DATES      = pd.date_range('2025-11-1', '2025-12-1', freq='5min')
    dates_f    = astro_time.convert_times( DATES, PA )
    original_f = sgp4.propTLE_df( dates_f, L1, L2, PA )
    to_j2k     = TEME_to_J2K( original_f.copy(), PA ) 

    from_j2k   = J2K_to_TEME( to_j2k, PA )

    X  = toCCSDS( from_j2k ) 
    Y  = fromCCSDS( X.split('\n') )
    print(Y)

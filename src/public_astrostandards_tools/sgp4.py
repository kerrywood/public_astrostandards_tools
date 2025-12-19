import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------------------------------
def addTLE( L1 : str, L2 : str, INTERFACE ):
    # load the TLE
    tleid = INTERFACE.TleDll.TleAddSatFrLines( 
        INTERFACE.Cstr( L1, 512 ),
        INTERFACE.Cstr( L2, 512 )
        )
    return tleid

# -----------------------------------------------------------------------------------------------------
def initTLE( tleid : float, INTERFACE ):
    return INTERFACE.Sgp4PropDll.Sgp4InitSat( tleid ) == 0

# -----------------------------------------------------------------------------------------------------
def propTLEToDS50s( tleid, ds50_l : list[ float ], INTERFACE ):
    '''
    take a tleID; it must be initialized and ready to go

    take a list of ds50 UTC values (dates in AstroStandard epoch) 
    and return <ds50><teme_pos (3)><teme_vel (3)>
    '''
    # data holders for OUTPUT
    ds50 = INTERFACE.ctypes.c_double() 
    pos  = (INTERFACE.ctypes.c_double * 3)()
    vel  = (INTERFACE.ctypes.c_double * 3)()
    llh  = (INTERFACE.ctypes.c_double * 3)()

    def propToDS50( dsutc ):
        INTERFACE.Sgp4PropDll.Sgp4PropDs50UtcPosVel( tleid, dsutc, pos, vel )
        return np.hstack( (dsutc, np.array(pos), np.array(vel) ) )
    return np.vstack( [propToDS50(X) for X in ds50_l ] ) 

# -----------------------------------------------------------------------------------------------------
def propTLE_byID_df( tleid, 
                tle_df,
                INTERFACE,
                clear_all = True ):
    '''
    this function assumes that your TLE has already been loaded into the AstroStandards
    pass in the ID.  This is useful when you're modifying a TLE via the array, or you
    don't want to parse and re-parse
    '''
    assert initTLE( tleid, INTERFACE ) 
    eph   = propTLEToDS50s( tleid, tle_df['ds50_utc'] , INTERFACE )
    tle_df['teme_p'] = eph[:,1:4].tolist()
    tle_df['teme_v'] = eph[:,4:7].tolist()
    return tle_df 


# -----------------------------------------------------------------------------------------------------
def propTLE_df( dates : pd.DataFrame, 
                line1 : str,
                line2 : str,
                INTERFACE,
                clear_all = True ):

    if clear_all :
        INTERFACE.TleDll.TleRemoveAllSats()
        INTERFACE.Sgp4PropDll.Sgp4RemoveAllSats()

    # add the TLE and init it
    tleid = addTLE( line1, line2, INTERFACE )    
    assert initTLE( tleid , INTERFACE )

    eph   = propTLEToDS50s( tleid, dates['ds50_utc'] , INTERFACE )
    rv    = dates
    rv['teme_p'] = eph[:,1:4].tolist()
    rv['teme_v'] = eph[:,4:7].tolist()
    return rv

# -----------------------------------------------------------------------------------------------------
def test():
    from . import astro_time
    from datetime import datetime, timedelta, timezone 

    import public_astrostandards as PA
    from public_astrostandards import helpers
    from . import test_helpers 

    PA.init_all()

    astro_time.load_time_constants(  test_helpers.get_test_time_constants(), PA )
    
    # test TLE
    L1 = '1 25544U 98067A   24365.67842578  .00026430  00000-0  46140-3 0  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    
    # set up your dates
    dates    = pd.date_range(   datetime.now( timezone.utc), 
                                datetime.now( timezone.utc ) + timedelta(days=1), 
                                freq='5 min' )
    time_df  = astro_time.convert_times( dates, PA )

    # now propagate into that frame
    testout =  propTLE_df( time_df.copy(), L1, L2, PA ) 
    print('-'*100)
    print('Type 0 test :\n\t{}\n\t{}'.format(L1,L2))
    print(testout)

    nL1 = '1 25545U 98067A   24365.67842578  .00026430  00000-0  46140-3 4  9990'
    nL2 = '2 25545  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    testID = addTLE( nL1, L2, PA )
    assert initTLE( testID, PA )
    print('-'*100)
    print('Fake type 4 test :\n\t{}\n\t{}'.format(nL1,nL2))
    print( propTLE_byID_df( testID, testout, PA ) )


# =====================================================================================================
if __name__ == '__main__':
    test()
    

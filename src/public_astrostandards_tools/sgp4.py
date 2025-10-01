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

# =====================================================================================================
if __name__ == '__main__':
    import astro_time
    from datetime import datetime, timedelta 

    import public_astrostandards as PA
    from public_astrostandards import helpers
    
    astro_time.load_time_constants(  '/opt/astrostandards/reduced_time_constants.dat' , PA )

    # test TLE
    L1 = '1 25544U 98067A   24365.67842578  .00026430  00000-0  46140-3 0  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    
    # set up your dates
    dates    = pd.date_range( datetime.utcnow(), datetime.utcnow() + timedelta(days=1), freq='5 min' )
    time_df  = astro_time.convert_times( dates, PA )

    # now propagate into that frame
    print( propTLE_df( time_df.copy(), L1, L2, PA ) )

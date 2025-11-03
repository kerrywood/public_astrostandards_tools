import ctypes 
from datetime import datetime, timedelta, timezone
import numpy as np
import pandas as pd
import copy

# -----------------------------------------------------------------------------------------------------
_DSEPOCH = datetime(year=1950,month=1,day=1)
def datetime_from_ds50( ds50 : float ):
    # minus 1 because astrostandards epoch is 1-indexed.. not 0
    return _DSEPOCH + timedelta( days=ds50-1 )

# -----------------------------------------------------------------------------------------------------
def getRIC( sv ):
    P = np.array( sv['teme_p'] )
    V = np.array(sv['teme_v'] )
    R = np.array(P) / np.linalg.norm(P)
    I = np.array(V) / np.linalg.norm(V)
    C = np.cross( I, R )
    return R, I, C

# -----------------------------------------------------------------------------------------------------
def XA_TLE_to_str( XA_TLE, harness, satno=None ):
    TSTR = harness.Cstr('',512)
    if satno : XA_TLE['XA_TLE_SATNUM'] = satno
    nL1, nL2 = harness.Cstr('',512),harness.Cstr('',512)
    harness.TleDll.TleGPArrayToLines( XA_TLE.data, TSTR, nL1, nL2 )
    return nL1.value.decode('utf-8').strip(), nL2.value.decode('utf-8').strip()

# -----------------------------------------------------------------------------------------------------
def sv_to_osc( sv, PA ):
    '''
    sv : <teme_pos><teme_vel>
    return XA_KEP
    '''
    # this is an expensive operation
    XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    # we'll use the conversion in the astrostandards
    PA.AstroFuncDll.PosVelToKep( 
        (ctypes.c_double*3)(*sv['teme_p']), 
        (ctypes.c_double*3)(*sv['teme_v']), 
        XA_KEP.data )
    return XA_KEP

#  -----------------------------------------------------------------------------------------------------
def osc_to_true_anomaly( xkep , PA ):
    return PA.AstroFuncDll.CompTrueAnomaly( xkep.data )

#  -----------------------------------------------------------------------------------------------------
def sv_to_osc_df( sv_df : pd.DataFrame, PA ) :
    ''' 
    given a dataframe with 'teme_p' and 'teme_v' on each row, annotate each row with XA_KEP data
    '''
    # getting a helpers object is expensive, so rather than re-use the sv_to_osc routine, we'll do 
    # something more efficient here and re-use the holder and a closure
    XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    def _kep_closure( D ):
        PA.AstroFuncDll.PosVelToKep( 
            (ctypes.c_double*3)(*D['teme_p']), 
            (ctypes.c_double*3)(*D['teme_v']), 
            XA_KEP.data )
        true_anom = PA.AstroFuncDll.CompTrueAnomaly( XA_KEP.data )
        rv = XA_KEP.toDict()
        rv.update( {'XA_KEP_TA' : true_anom } )
        return rv
    tv = sv_df.apply( lambda X: _kep_closure(X), axis=1 )
    tv = pd.DataFrame.from_dict( tv.values.tolist() )
    rv = pd.concat( (sv_df.reset_index(drop=True), tv.reset_index(drop=True) ), axis=1 )
    # add in the true anomaly data 
    return rv

# -----------------------------------------------------------------------------------------------------
def osc_to_mean( XA_KEP, PA ):
    '''
    take a XA_KEP structure and return a XA_KEP with *mean* fields
    '''
    XA_KEP_MEAN = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.KepOscToMean( XA_KEP.data, XA_KEP_MEAN.data )
    return XA_KEP_MEAN

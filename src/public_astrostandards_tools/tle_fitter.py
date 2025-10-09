import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import scipy.optimize
import sgp4

# what fields will we optimize over?  This doubles as a field accessor list for the optimizer..
FIT_TYPE4 = [
        'XA_TLE_BTERM',
        'XA_TLE_NDOT',
        'XA_TLE_AGOM_GP',
        'XA_TLE_INCLI',
        'XA_TLE_NODE',
        'XA_TLE_ECCEN',
        'XA_TLE_OMEGA',
        'XA_TLE_MNANOM',
        'XA_TLE_MNMOTN',
    ]

FIT_CORE = [
    'XA_TLE_INCLI',
    'XA_TLE_NODE',
    'XA_TLE_ECCEN',
    'XA_TLE_OMEGA',
    'XA_TLE_MNANOM',
    'XA_TLE_MNMOTN',
    ]

FIT_TYPE0 = [
    'XA_TLE_INCLI',
    'XA_TLE_NODE',
    'XA_TLE_ECCEN',
    'XA_TLE_OMEGA',
    'XA_TLE_MNANOM',
    'XA_TLE_MNMOTN',
    'XA_TLE_BSTAR',
    'XA_TLE_BSTAR',
    ]
#FIELDS = _FULL_FIELDS
FIELDS = FIT_TYPE0

# -----------------------------------------------------------------------------------------------------
def XA_TLE_to_str( XA_TLE, PA ):
    PA.TleDll.TleRemoveAllSats()
    XS_TLE =  PA.Cstr('',512)
    tleid = PA.TleDll.TleAddSatFrArray( XA_TLE.data, XS_TLE )
    assert tleid > 0
    outL1, outL2 = PA.Cstr('',512), PA.Cstr('',512)
    assert PA.TleDll.TleGetLines( tleid, outL1, outL2 ) == 0
    return outL1.value.decode('utf-8').strip(), outL2.value.decode('utf-8').strip()

# -----------------------------------------------------------------------------------------------------
def TLE_str_to_XA_TLE( L1 : str, L2 : str , PA ):
    # load the TLE
    PA.TleDll.TleRemoveAllSats()
    tleid = PA.TleDll.TleAddSatFrLines( PA.Cstr(L1,512), PA.Cstr(L2,512) )
    if tleid <=0 : return None
    XA_TLE = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
    XS_TLE = PA.Cstr('',512)
    PA.TleDll.TleDataToArray( tleid, XA_TLE.data, XS_TLE )  # <--- note that you pass the "data" holder in
    return XA_TLE

# -----------------------------------------------------------------------------------------------------
def sv_to_osc( sv, PA ):
    '''
    sv : <teme_pos><teme_vel>
    return XA_KEP
    '''
    # we'll use the conversion in the astrostandards
    XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.PosVelToKep( 
        (ctypes.c_double*3)(*sv['teme_p']), 
        (ctypes.c_double*3)(*sv['teme_v']), 
        XA_KEP.data )
    return XA_KEP

# -----------------------------------------------------------------------------------------------------
def osc_to_mean( XA_KEP, PA ):
    XA_KEP_MEAN = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.KepOscToMean( XA_KEP.data, XA_KEP_MEAN.data )
    return XA_KEP_MEAN

# -----------------------------------------------------------------------------------------------------
def insert_kep_to_TLE( TLE_DATA, KEP_DATA, PA ):
    # now set the values from KEP_DATA into TLE_DATA
    TLE_DATA['XA_TLE_INCLI']  = KEP_DATA['XA_KEP_INCLI']
    TLE_DATA['XA_TLE_NODE']   = KEP_DATA['XA_KEP_NODE']
    TLE_DATA['XA_TLE_ECCEN']  = KEP_DATA['XA_KEP_E']
    TLE_DATA['XA_TLE_MNANOM'] = KEP_DATA['XA_KEP_MA']
    TLE_DATA['XA_TLE_OMEGA']  = KEP_DATA['XA_KEP_OMEGA']
    TLE_DATA['XA_TLE_MNMOTN'] = PA.AstroFuncDll.AToN( KEP_DATA['XA_KEP_A'] )
    return TLE_DATA

# -----------------------------------------------------------------------------------------------------
class tle_fitter:
    def __init__( self, PA ):
        self.PA         = PA        # this is the harness for public_astrostandards
        self.FIELDS     = FIELDS    # what fields are we optimizing over (from XA_TLE)
        self.init_tle   = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
        self.new_tle    = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
        self.satno      = None
        self.epoch_idx  = None
    
    def set_type0( self ):
        self.new_tle['XA_TLE_EPHTYPE'] = 0
        return self

    def set_type2( self ):
        self.new_tle['XA_TLE_EPHTYPE'] = 2
        return self

    def set_type4( self ):
        self.new_tle['XA_TLE_EPHTYPE'] = 2
        return self

    def set_satno( self, satno=99999 ):
        self.satno = satno
        self.new_tle['XA_TLE_SATNUM'] = self.satno
        return self

    def get_init_fields( self ):
        return [ self.init_tle[X] for X in self.FIELDS ]

    def reset_tle( self ):
        for k in self.FIELDS:
            self.new_tle[k] = self.init_tle[k]
        return self


# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import sgp4


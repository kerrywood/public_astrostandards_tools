import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import scipy.optimize
from . import sgp4
from . import orbit_utils

# what fields will we optimize over?  This doubles as a field accessor list for the optimizer..
FIT_TYPE4 = [
        'XA_TLE_BTERM',
        'XA_TLE_NDOT',
        'XA_TLE_AGOMGP',
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
    ]

NON_CONSERVATIVES = [
    'XA_TLE_AGOMGP',
    'XA_TLE_NDOT',
    'XA_TLE_NDOTDOT',
    'XA_TLE_BSTAR',
    'XA_TLE_BTERM',
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
    return XA_TLE, XS_TLE 

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
        self.init_str   = PA.Cstr('',512)
        self.new_tle    = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
        self.satno      = None
        self.epoch_idx  = None
        self.tle_type   = 0
        self.mm_kozai   = None      # we might switch between these; store both "out of band" (not in the data array)
        self.mm_brouwer = None

    # when we set TLE data from some external source (like state vector), compute the Kozai and 
    # Brouwer elements for those data; we can then set the fields when we switch types
    def _brouwer_to_kozai( self ):
        self.mm_kozai = self.PA.AstroFuncDll.BrouwerToKozai( 
                                self.init_tle['XA_TLE_ECCEN'], 
                                self.init_tle['XA_TLE_INCLI'],
                                self.init_tle['XA_TLE_MNMOTN'] )

    def _kozai_to_brouwer( self ):
        self.mm_brouwer = self.PA.AstroFuncDll.KozaiToBrouwer( 
                                self.init_tle['XA_TLE_ECCEN'], 
                                self.init_tle['XA_TLE_INCLI'],
                                self.init_tle['XA_TLE_MNMOTN'] )
    
    def set_from_lines( self, L1 : str, L2 : str ):
        self.init_tle, self.init_str = TLE_str_to_XA_TLE( L1, L2, self.PA )
        self.tle_type  = self.init_tle['XA_TLE_EPHTYPE']
        if self.tle_type:
            self.mm_kozai   = self.init_tle['XA_TLE_MNMOTN']
            self._kozai_to_brouwer()
        else:
            self.mm_brouwer = self.init_tle['XA_TLE_MNMOTN']
            self._brouwer_to_kozai()
            
    def set_type0( self ):
        self.tle_type = 0
        self.FIELDS   = FIT_TYPE0
        self.init_tle['XA_TLE_EPHTYPE'] = self.tle_type
        self.init_tle['XA_TLE_MNMOTN']  = self.mm_kozai
        return self

    def set_type2( self ):
        self.tle_type = 2
        self.FIELDS   = FIT_TYPE0  # <--- fields are the same as TYPE0
        self.init_tle['XA_TLE_EPHTYPE'] = self.tle_type
        self.init_tle['XA_TLE_MNMOTN']  = self.mm_brouwer
        return self

    def set_type4( self ):
        self.tle_type = 4
        self.FIELDS   = FIT_TYPE4
        self.new_tle['XA_TLE_EPHTYPE'] = self.tle_type
        self.new_tle['XA_TLE_MNMOTN']  = self.mm_kozai
        return self

    def set_AGOM( self, agom=0.01 ):
        self.init_tle['XA_TLE_AGOMGP'] = agom
        return self

    def set_satno( self, satno=99999 ):
        self.satno = satno
        self.new_tle['XA_TLE_SATNUM'] = self.satno
        return self

    def get_init_fields( self ):
        # if we are changing type-to-type, we might not get all the field names right
        return [ self.init_tle[X] for X in self.FIELDS ]

    def reset_tle( self ):
        for k in self.FIELDS:
            self.new_tle[k] = self.init_tle[k]
        return self

    def clear_nonconservatives( self ):
        for F in NON_CONSERVATIVES: self.new_tle[F] = 0.
        return self

    def set_from_sv( self, sv ):
        # set the date
        self.new_tle['XA_TLE_EPOCH'] = sv['ds50_utc']
        # set the elements according to the state vector
        osc  = orbit_utils.sv_to_osc( sv, self.PA )
        mean = orbit_utils.osc_to_mean( osc, self.PA )
        insert_kep_to_TLE( self.new_tle, mean, self.PA )
        self._brouwer_to_kozai()

    def getLines( self ):
        return XA_TLE_to_str( self.new_tle, self.PA )
    
    def initial_simplex( self, delta=0.15):
        '''
        take our initial fields and perturb each entry delta% in either direction (up and down
        this should give us a good search space
        '''
        X     = self.get_init_fields()
        # +/5 15% on the above diagonal
        smplx = np.ones( shape=( len(X), len(X) ) )
        smplx += np.diag( np.ones( len(X)-1 ), -1 ) * -delta
        smplx += np.diag( np.ones( len(X)-1 ), 1 ) * delta
        smplx = np.vstack( (np.ones(len(X)), smplx ) )

        # random normal ...
        #smplx = np.random.normal( 1, 0.2, size=( len(X)+1, len(X) ) )

        smplx *= X
        return smplx


# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import sgp4

    A = tle_fitter()


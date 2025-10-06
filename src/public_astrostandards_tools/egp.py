import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import scipy.optimize
import sgp4

# what fields will we optimize over?  This doubles as a field accessor list for the optimizer..
_FULL_FIELDS = [
    'XA_TLE_BTERM',
    'XA_TLE_NDOT',
    'XA_TLE_SP_BTERM',
    'XA_TLE_INCLI',
    'XA_TLE_NODE',
    'XA_TLE_ECCEN',
    'XA_TLE_OMEGA',
    'XA_TLE_MNANOM',
    'XA_TLE_MNMOTN',
]

_REDUCED_FIELDS = [
    'XA_TLE_INCLI',
    'XA_TLE_NODE',
    'XA_TLE_ECCEN',
    'XA_TLE_OMEGA',
    'XA_TLE_MNANOM',
    'XA_TLE_MNMOTN',
]

#FIELDS = _FULL_FIELDS
FIELDS = _REDUCED_FIELDS

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
def optFunction( X, EH, return_scalar=True ):
    PA      = EH.PA
    XS_TLE  = PA.Cstr('',512)
    # take the function parameters (X) and overwrite the "new_tle" values based on FIELDS 
    for k,v in zip(FIELDS,X) :  EH.new_tle[ k ] = v
    # --------------------- clear state
    PA.TleDll.TleRemoveAllSats()
    PA.Sgp4PropDll.Sgp4RemoveAllSats()
    # --------------------- init our test TLE from the modified data
    tleid = PA.TleDll.TleAddSatFrArray( EH.new_tle.data, XS_TLE )
    if tleid <= 0: return np.inf
    if PA.Sgp4PropDll.Sgp4InitSat( tleid ) != 0: return np.inf
    # --------------------- generate our test ephemeris
    test_eph = sgp4.propTLEToDS50s( tleid, EH.truth_date, PA )
    # use numpy to return the distance between our hypothesis and truth
    resids = test_eph[:,1:4] - EH.truth_eph
    rms    = np.sqrt( np.sum( np.linalg.norm( resids, axis=1 ) ) / resids.shape[0] )
    print( 'RMS : {:08.3f}                    '.format(rms) , end='\r')
    if return_scalar:
        return rms
        # return np.sum( np.linalg.norm( resids[:,:3], axis=1 ) ) 
    else:
        np.linalg.norm( resids[:,:3], axis=1 ) 

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
class egp_helper:
    def __init__( self, PA ):
        self.PA         = PA        # this is the harness for public_astrostandards
        self.FIELDS     = FIELDS    # what fields are we optimizing over (from XA_TLE)
        self.init_tle   = None      # holder for the initial TLE data 
        self.new_tle    = None      # holder for the TLE data we'll modify
        self.satno      = None
        self.epoch_idx  = None
    
    def init( self ):
        # pick one of those points as the new epoch (put it in the middle)   # <----- epoch and sv choice
        if self.epoch_idx == None:
            self.epoch_idx = self.truth_df.shape[0]//2 
        epoch_sv        = self.truth_df.iloc[ self.epoch_idx ]
        osc_data        = sv_to_osc( epoch_sv , self.PA )
        mean_data       = osc_to_mean( osc_data, self.PA )
        # update our "fit" TLE with the new osculating data
        self.new_tle    = insert_kep_to_TLE( self.new_tle, mean_data, self.PA )
        # if this is a type-0, we need Kozai mean motion   
        if self.new_tle['XA_TLE_EPHTYPE'] == 0:
            self.new_tle['XA_TLE_MNMOTN'] = self.PA.AstroFuncDll.BrouwerToKozai( 
                    self.new_tle['XA_TLE_ECCEN'], 
                    self.new_tle['XA_TLE_INCLI'],
                    self.new_tle['XA_TLE_MNMOTN'] )
        self.new_tle['XA_TLE_EPOCH'] = epoch_sv['ds50_utc']
        return self

    
    def set_data( self, L1 : str, L2 : str, dates : list[ datetime ] ):
        ''' 
        take a TLE (from lines) and a set of dates we'll optimize over,
        setup everything we need for a fit
        '''
        dates_f = astro_time.convert_times( dates, self.PA )
        # crack open this TLE
        self.init_tle   = TLE_str_to_XA_TLE( L1, L2, self.PA )
        #self.new_tle    = TLE_str_to_XA_TLE( L1, L2, self.PA )
        self.new_tle    = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
        self.set_satno()
        # clear all sats
        self.PA.TleDll.TleRemoveAllSats()
        # propagate this TLE to the dates.. this is our truth s
        self.truth_df   = sgp4.propTLE_df( dates_f, L1, L2, self.PA )
        self.truth_dt   = dates_f['datetime'].values
        self.truth_date = dates_f['ds50_utc'].values
        self.truth_eph  = np.vstack( (self.truth_df['teme_p']) )
        self.init()
        return self

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

    def egp_tle( self ):
        # -----------------------------  nelder mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        #X = EH.get_init_fields()
        #smplx = np.random.normal( 1, 0.15, size=( len(X) + 1, len(X) ) )
        #smplx *= X
        ans = scipy.optimize.minimize( optFunction, 
                                       self.get_init_fields(),
                                       args    = (EH,True),
                                       method  = 'Nelder-Mead' )
                                       #options = {'initial_simplex' : smplx } )
        if ans.success:
            self.reset_tle()
            # now update with perturbed values (not sure if this is necessary... last step should be there)
            for k,v in zip(FIELDS,ans.x) : EH.new_tle[k] = v
            return  {   'input_tle'  : (L1,L2),
                        'output_tle' : XA_TLE_to_str( self.new_tle, PA ),
                        'fields'     : self.FIELDS,
                        'dates'      : self.truth_dt,
                        'points'     : len( self.truth_date),
                        'rms'        : ans.fun,
                        'rms_per_pt' : ans.fun / len( self.truth_date ) }

    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import sgp4

    PA.init_all()
    # example TLE 
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    # this is your fit range
    #DATES = pd.date_range( '2025-1-7', '2025-1-9', freq='5min' )
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )

    # setup the job
    EH = egp_helper( PA ).set_data(L1, L2, DATES ).set_satno(77777).set_type0()

    print()
    print('-'*100)
    print('Fitting :\n\t{}\n\t{}\n\t{}'.format( L1, L2, DATES ) )
    print('-'*100)
    print()
    # holder

    output = EH.egp_tle( )

    print()
    print('Your original TLE was')
    print('\n'.join( [L1,L2] ) ) 


    print('Your new TLE is :')
    print('\n'.join(output['output_tle'] ) )
    
    #print()
    #print('-'*100)
    #for k,v in output.items():
    #    print('{:13} {}'.format(k,v))
    #print('-'*100)

import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import scipy.optimize
import astro_time
import sgp4
from enum import Enum

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
    for k,v in zip(EH.FIELDS,X) :  EH.new_tle[ k ] = v
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
    def __init__( self, PA, FIELDS=FIT_CORE ):
        self.PA         = PA        # this is the harness for public_astrostandards
        self.FIELDS     = FIELDS    # what fields are we optimizing over (from XA_TLE)
        self.init_tle   = None      # holder for the initial TLE data 
        self.new_tle    = None      # holder for the TLE data we'll modify
        self.satno      = None
        self.epoch_idx  = None
        self.init_l1    = None
        self.init_l2    = None
    
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
        self.init_l1    = L1
        self.init_l2    = L2
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
        self.init()
        return self

    def set_type2( self ):
        self.new_tle['XA_TLE_EPHTYPE'] = 2
        self.init()
        return self

    def set_type4( self ):
        self.new_tle['XA_TLE_EPHTYPE'] = 4
        self.init()
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

    def get_init_tle( self ):
        return self.init_tle, XA_TLE_to_str( self.init_tle, self.PA )
    
    def get_new_tle( self ):
        return self.new_tle, XA_TLE_to_str( self.new_tle, self.PA )

    def summarize_results( self ):
        ntle, text = self.get_new_tle()
        rv = ntle.toDict()
        rv.update( {'new_line1' : text[0], 'new_line2' : text[1] })
        if self.init_l1 and self.init_l2 : 
            rv.update({'ini_line1' : self.init_l1, 'ini_line2' : self.init_l2 })
        return rv

    def egp_tle( self ):
        # -----------------------------  nelder mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        X = EH.get_init_fields()
        smplx = np.random.normal( 1, 0.15, size=( len(X) + 1, len(X) ) )
        smplx *= X
        ans = scipy.optimize.minimize( optFunction, 
                                       self.get_init_fields(),
                                       args    = (self,True),
                                       method  = 'Nelder-Mead' ,
                                       options = {'fatol' : 0.2, 'xatol' : 0.01, 'initial_simplex' : smplx} )
        if ans.success:
            self.reset_tle()
            # now update with perturbed values (not sure if this is necessary... last step should be there)
            for k,v in zip(self.FIELDS,ans.x) : self.new_tle[k] = v
            return self
        print(ans)


# ----------------------------------------------------------------------------------------------------- 
def test():
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import sgp4

    print()
    print('-'*100)
    print('Running tests')
    print('-'*100)
    PA.init_all()
    # example TLE 
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    # this is your fit range
    #DATES = pd.date_range( '2025-1-7', '2025-1-9', freq='5min' )
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )

    # type 4 to type 0
    # setup the job
    print()
    print('Your original TLE was')
    print('\n'.join( [L1,L2] ) ) 
    print('Fitting {} points from {} -- {}'.format( len(DATES), DATES[0], DATES[-1]) )

    EH = egp_helper( PA ).set_data(L1, L2, DATES ).set_satno(77777).set_type0()
    output = EH.egp_tle( )
    print()
    print('\nType 0 fit:')
    print(EH.summarize_results())

    print('\nType 2 fit:')
    EH.set_type2()
    output = EH.egp_tle( )
    print(EH.summarize_results())

    print('\nType 4 fit:')
    EH.set_type4()
    output = EH.egp_tle( )
    print(EH.summarize_results())

    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import sys
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import argparse
    import json

    PA.init_all()

    parser = argparse.ArgumentParser(
        prog='egp converter',
        description='''take a type-4 TLE (or really any type) and convert it to a fit region''',
        epilog = ''
    )

    parser.add_argument('--line1',"-l1",
            required=False,
            help='line 1 of a TLE', 
            default='1 25544U 98067A   24365.67842578    .00000000  00000-0  00000-0 4  9990')

    parser.add_argument('--line2','-l2',
            required=False,
            help = 'line2 of TLE',
            default = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028' )

    parser.add_argument('--startdate',"-sd", required=False,                    help='first date (in a format that pd.date_range can understand)')
    parser.add_argument('--enddate',  "-ed", required=False,                    help='last date (in a format that pd.date_range can understand)')
    parser.add_argument('--spacing',  "-sp", required=False, default='5min',    help='ephemeris spacing (in format pd.date_range can understand)')
    parser.add_argument('--type',     "-ty", required=False, default='0',       help='TLE type (0,2,4)', type=int)
    parser.add_argument('--satno',    "-sno",required=False, default=99999,     help='new TLE satno',type=int)
    parser.add_argument('--test',     "-ts", required=False, action='store_true',help='run a test (all other arguments ignored)')

    # parse the arguments
    args = parser.parse_args()

    if args.test:
        test()
        sys.exit()

    # now start adding params
    assert  args.type  in set([0,2,4])
    dates = pd.date_range( args.startdate, args.enddate, freq=args.spacing )

    EH = egp_helper( PA ).set_data( args.line1, args.line2, dates )
    if args.type == 0: EH.set_type0()
    if args.type == 2: EH.set_type2()
    if args.type == 4: EH.set_type4()
    EH.set_satno( args.satno )

    #    test()
    print()
    print()
    print( json.dumps( EH.egp_tle().summarize_results(), indent=4) )

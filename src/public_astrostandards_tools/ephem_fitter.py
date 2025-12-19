from datetime import datetime, timedelta, timezone
import numpy as np
import pandas as pd
import scipy.optimize

from . import astro_time
from . import tle_fitter
from . import sgp4


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
class ephem_fitter( tle_fitter.tle_fitter ):
    def __init__( self, PA ):
        super().__init__( PA )
        self.line1 = None
        self.line2 = None

    def _init_fields( self ):
        # pick one of those points as the new epoch (put it in the middle)   # <----- epoch and sv choice
        if self.epoch_idx == None:
            self.epoch_idx = self.truth_df.shape[0]//2 
        epoch_sv        = self.truth_df.iloc[ self.epoch_idx ]
        self.set_from_sv( epoch_sv )
        self.new_tle['XA_TLE_EPOCH'] = epoch_sv['ds50_utc']
        return self

    def set_ephemeris( self, eph_df : pd.DataFrame ):
        self.truth_df   = eph_df
        self.truth_dt   = eph_df['datetime'].values
        self.truth_date = eph_df['ds50_utc'].values
        self.truth_eph  = np.vstack( (self.truth_df['teme_p']) )
        self._init_fields()

    def set_from_tle( self, L1 : str, L2 : str, dates : list[ datetime ] ):
        ''' 
        take a TLE (from lines) and a set of dates we'll optimize over,
        setup everything we need for a fit
        '''
        dates_f = astro_time.convert_times( dates, self.PA )
        # crack open this TLE
        self.set_from_lines( L1, L2 )
        # clear all sats
        self.PA.TleDll.TleRemoveAllSats()
        # propagate this TLE to the dates.. this is our truth s
        eph_df = sgp4.propTLE_df( dates_f, L1, L2, self.PA )
        self.set_ephemeris( eph_df )
        self.line1      = L1
        self.line2      = L2
        return self

    def summarize_results( self ):
        return  {   'input_tle'  : (self.line1,self.line2),
                    'output_tle' : tle_fitter.XA_TLE_to_str( self.new_tle, self.PA ),
                    'fields'     : self.FIELDS,
                    'dates'      : ( self.truth_dt[0], self.truth_dt[-1] ),
                    'points'     : len( self.truth_date),
                    'rms'        : float(self.ans.fun),
                    'rms_per_pt' : float(self.ans.fun) / len( self.truth_date ) }

    def fit_tle( self ):
        # -----------------------------  nelder mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        ans   = scipy.optimize.minimize(optFunction, 
                                        self.get_init_fields(),
                                        args    = (self,True),
                                        method  = 'Nelder-Mead' ,
                                        #options = {'xatol' : 0.01, 'fatol' : 0.9 } )
                                        options = {'xatol' : 0.01, 'fatol' : 0.1, 'initial_simplex' : self.initial_simplex() } )
                                       #options = {'initial_simplex' : smplx } )
        if ans.success:
            self.reset_tle()
            # now update with perturbed values (not sure if this is necessary... last step should be there)
            for k,v in zip(self.FIELDS,ans.x) : self.new_tle[k] = v
            self.ans = ans
            return self

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    PA.init_all()

    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'

    print('-'*100)
    print('Performing fit test')
    print()
    # example TLE 
    # this is your fit range
    # DATES = pd.date_range( '2025-1-7', '2025-1-9', freq='5min' )
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )


    # setup the job
    EH = ephem_fitter( PA ).set_from_tle(L1, L2, DATES ).set_satno(77777).set_type0()

    print('-'*100)
    print('Fitting :\n\t{}\n\t{}'.format( L1, L2 ) )
    print('\tto: {} -- {}'.format( DATES[0], DATES[-1] ) )
    print('\t{} points'.format( len(DATES) ) )
    print('-'*100)
    print()

    output = EH.fit_tle( )

    print()
    print('\nYour original TLE was')
    print('\n'.join( [L1,L2] ) ) 


    print('\nYour new TLE is :')
    print('\n'.join( output.getLines() ) )
    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    test()

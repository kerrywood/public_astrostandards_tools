from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import scipy.optimize

import astro_time
import sgp4
import sensor
import observations
import tle_fitter

# -----------------------------------------------------------------------------------------------------
def optFunction( X, EH, return_scalar=True ):
    XS_TLE  = EH.PA.Cstr('',512)
    # take the function parameters (X) and overwrite the "new_tle" values based on FIELDS 
    for k,v in zip(EH.FIELDS,X) :  EH.new_tle[ k ] = v
    # --------------------- clear state
    EH.PA.TleDll.TleRemoveAllSats()
    EH.PA.Sgp4PropDll.Sgp4RemoveAllSats()
    # --------------------- init our test TLE from the modified data
    tleid = PA.TleDll.TleAddSatFrArray( EH.new_tle.data, XS_TLE )
    if tleid <= 0: return np.inf
    if EH.PA.Sgp4PropDll.Sgp4InitSat( tleid ) != 0: return np.inf
    # --------------------- generate our test ephemeris
    test_frame    = sgp4.propTLE_byID_df( tleid, EH.date_f, EH.PA )
    looks         = sensor.compute_looks( EH.sensor_df, test_frame, EH.PA )
    RA_residuals  = EH.obs_df['teme_ra'] - looks['XA_TOPO_RA']
    DEC_residuals = EH.obs_df['teme_dec'] - looks['XA_TOPO_DEC']
    # shortest angle / path computation ( -180 and 180 are NOT 360 apart )
    RA_residuals  = (RA_residuals + 180) % 360 - 180
    DEC_residuals = (DEC_residuals + 180) % 360 - 180
    # compute the RMS based on the residuals (https://stackoverflow.com/questions/1878907/how-can-i-find-the-smallest-difference-between-two-angles-around-a-point#7869457)
    N   = len(RA_residuals)
    rms = np.sum( RA_residuals.values ** 2 + DEC_residuals.values ** 2 ) / ( 2 * N )
    rv  = np.sqrt( rms )
    print('{:10.7f}                '.format(rv), end='\r')
    if return_scalar : return rv
    return {'observations'  : EH.obs_df.to_dict(orient='records'), 
            'looks'         : looks.to_dict(orient='records'), 
            'initial_line1' : EH.line1, 
            'initial_line2' : EH.line2 }

# -----------------------------------------------------------------------------------------------------
class eo_fitter( tle_fitter.tle_fitter ):
    def __init__( self, PA ):
        super().__init__( PA )
        self.line1 = None
        self.line2 = None

    def _set_new_epoch( self, epoch ):
        ''' assume that epoch is set, and that line1, line2 are also set '''
        self.PA.TleDll.TleRemoveAllSats()
        tleid = sgp4.addTLE( self.line1, self.line2, self.PA )
        assert sgp4.initTLE( tleid, self.PA )
        rv    = sgp4.propTLEToDS50s( tleid, [ epoch ], self.PA )[0]
        rv    = { 'teme_p' : rv[1:4], 'teme_v' : rv[4:7], 'ds50_utc' : self.epoch_ds50 }
        return rv

    def set_data( self, L1 : str, L2 : str, inobs : list[ dict ] ):
        ''' 
        take an initial TLE as a guess (L1,L2) 
        take a list of JSON formatted obs (directly from UDL)
        solve for a new TLE
        '''
        self.line1      = L1
        self.line2      = L2
        self.obs        = inobs

        # everything builds off obs; set up the frame and pull off the key date fields
        self.obs_df     = observations.prepObs( inobs, self.PA )
        self.obs_df     = observations.rotateTEMEdf( self.obs_df, self.PA )
        self.date_f     = self.obs_df[ ['ds50_utc','ds50_et','theta']].copy()

        # init the TLE from the lines data
        self.init_tle    = tle_fitter.TLE_str_to_XA_TLE( L1, L2, self.PA )
        #self.new_tle     = tle_fitter.TLE_str_to_XA_TLE( L1, L2, self.PA )
        # pick the last ob as the epoch of the TLE
        self.epoch_ds50  = self.obs_df.iloc[ -1 ]['ds50_utc']
        self.epoch_dt    = self.obs_df.iloc[ -1 ]['obTime_dt']
        epoch_sv         = self._set_new_epoch( self.epoch_ds50 )
        self.set_from_sv( epoch_sv )
        # if this is a type-0, we need Kozai mean motion   
        if self.new_tle['XA_TLE_EPHTYPE'] == 0:
            self.new_tle['XA_TLE_MNMOTN'] = self.PA.AstroFuncDll.BrouwerToKozai( 
                                                self.new_tle['XA_TLE_ECCEN'], 
                                                self.new_tle['XA_TLE_INCLI'],
                                                self.new_tle['XA_TLE_MNMOTN'] )

        # setup the sensor frame (for generating looks)
        self.sensor_df        = self.obs_df[['ds50_utc','senlat','senlon','senalt','theta']]
        self.sensor_df        = self.sensor_df.rename( columns = {'senlat' : 'lat','senlon' : 'lon', 'senalt' : 'height' } )
        self.sensor_df        = sensor.llh_to_eci( self.sensor_df, self.PA )
        return self


    def fit_tle( self ):
        # -----------------------------  nelder mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        # TODO : termination conditions should be set AND we should weight according to obs calibration
        ans   = scipy.optimize.minimize(optFunction, 
                                        self.get_init_fields(),
                                        args    = (self,True),
                                        method  = 'Nelder-Mead' ,
                                        options = {'xatol' : 0.10, 'fatol' : 30/3600 * self.obs_df.shape[0], 'initial_simplex' : self.initial_simplex() } )

        self.ans = ans
        if ans.success:
            self.reset_tle()
            # now update with perturbed values (not sure if this is necessary... last step should be there)
            for k,v in zip(self.FIELDS,ans.x) : self.new_tle[k] = v
            self.final_answer = optFunction( ans.x, self, False )
            return True
        return False


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

    parser = argparse.ArgumentParser(
        prog='UDL EO obs fitter',
        description='''take a set of obs downloaded from UDL and an initial TLE, and fit it''',
        epilog = ''
    )

    parser.add_argument('--line1',"-l1",
            required = True,
            help     = 'line 1 of a TLE')

    parser.add_argument('--line2','-l2',
            required = True,
            help     = 'line2 of TLE')

    parser.add_argument('--infile',   "-F", 
            required = True, 
            default  = './19548.json.gz',
            help     ='load obs from this file')

    parser.add_argument('--outfile',  "-O", 
            required = True,
            default  = './out.json',
            help     = 'store output from loaded jobs')

    parser.add_argument('--type',  "-T",
            required = False, 
            type     = int,
            help     = 'TLE type to fit (0,2,4)')
    
    parser.add_argument('--verbose',  "-v",
            required = False, 
            default  = False,
            action   = 'store_true',
            help     = 'print debugging info')
    

    # parse the arguments
    args = parser.parse_args()

    # init the public_astrostandards harness 
    PA.init_all()
    astro_time.load_time_constants('./reduced_time_constants.dat', PA )

    # init the fitter
    FIT = eo_fitter( PA ) 
    
    # check on type
    if type in args:
        assert args.type in set([0,2,4])
        if args.type == 0 : FIT.set_type0()
        if args.type == 2 : FIT.set_type2()
        if args.type == 4 : FIT.set_type4()

    # load up the obs
    obs    = pd.read_json( args.infile ).sort_values(by='obTime').reset_index(drop=True)    

    # add in the data
    FIT = FIT.set_data( args.line1, args.line2, obs )
    FIT = FIT.set_satno(99999)

    if args.verbose:
        print('Fitting')
        print('initial TLE:\n\t{}\n\t{}'.format( FIT.line1, FIT.line2 ) )
        print('\tnew epoch : {}'.format( FIT.epoch_dt ) )
        print('\tfitting over {} obs'.format( len(obs) ) )
        print('\tspan {} -- {}'.format( FIT.obs_df.iloc[0]['obTime_dt'],
                                         FIT.obs_df.iloc[-1]['obTime_dt'] ) ) 

    # do the solve
    converged = FIT.fit_tle()
    if converged:
        if args.verbose: rv = FIT.final_answer
        else: rv = {}
        nl1, nl2 = FIT.getLines()
        print('New TLE:\n\t{}\n\t{}'.format( nl1, nl2 ) )
        rv.update({ 'new_line1' : nl1, 'new_line2' : nl2 })
        with open( args.outfile, 'wt') as F: json.dump( rv, F, default=str )
    else:
        print('Did not converge')

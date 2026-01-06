import numpy as np
import pandas as pd
import scipy.optimize

from . import astro_time
from . import sgp4
from . import sensor
from . import observations
from . import tle_fitter
from . import residuals

# -----------------------------------------------------------------------------------------------------
def optFunction( X, EH, return_scalar=True ):
    XS_TLE  = EH.PA.Cstr('',512)
    # take the function parameters (X) and overwrite the "new_tle" values based on FIELDS 
    for k,v in zip(EH.FIELDS,X) : 
        EH.new_tle[ k ] = v
    # --------------------- clear state
    EH.PA.TleDll.TleRemoveAllSats()
    EH.PA.Sgp4PropDll.Sgp4RemoveAllSats()
    # --------------------- init our test TLE from the modified data
    tleid = EH.PA.TleDll.TleAddSatFrArray( EH.new_tle.data, XS_TLE )
    if tleid <= 0: 
        return np.inf
    if EH.PA.Sgp4PropDll.Sgp4InitSat( tleid ) != 0: 
        return np.inf
    # --------------------- generate our test ephemeris
    target_frame  = sgp4.propTLE_byID_df( tleid, EH.date_f, EH.PA )
    # --------------------- generate looks from our sensor positinos
    looks         = sensor.compute_looks( EH.sensor_df, target_frame, EH.PA )
    # --------------------- get the residuals of these frames / obs
    resids        = residuals.UDL_residuals( EH.obs_df, looks )
    N   = resids.shape[0]
    rms = np.sum( resids['ra'].values ** 2 + resids['dec'].values**2 ) / (2 * N )
    rv  = np.sqrt( rms )
    print('RMS : {:10.7f}                '.format(rv), end='\r')
    if return_scalar : 
        return rv
    return {'observations'  : EH.obs_df.to_dict(orient='records'), 
            'looks'         : looks.to_dict(orient='records'), 
            'initial_line1' : EH.line1, 
            'initial_line2' : EH.line2,
            'residuals'     : resids.to_dict( orient='records') }

# -----------------------------------------------------------------------------------------------------
class eo_fitter( tle_fitter.tle_fitter ):
    def __init__( self, PA ):
        super().__init__( PA )
        self.line1 = None
        self.line2 = None

    def _move_epoch( self, epoch ):
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
        self.obs_df     = observations.prepUDLObs( inobs, self.PA )
        self.date_f     = self.obs_df[ ['ds50_utc','ds50_et','theta']].copy()

        # init the TLE from the lines data
        self.set_from_lines( L1, L2 )

        # pick the last ob as the epoch of the TLE (propagate your search TLE and set up data)
        self.epoch_ds50  = self.obs_df.iloc[ -1 ]['ds50_utc']
        self.epoch_dt    = self.obs_df.iloc[ -1 ]['obTime_dt']
        epoch_sv         = self._move_epoch( self.epoch_ds50 )
        self.set_from_sv( epoch_sv )

        # setup the sensor frame (for generating looks)
        self.sensor_df        = self.obs_df[['ds50_utc','senlat','senlon','senalt','theta']]
        self.sensor_df        = self.sensor_df.rename( columns = {'senlat' : 'lat','senlon' : 'lon', 'senalt' : 'height' } )
        self.sensor_df        = sensor.llh_to_eci( self.sensor_df, self.PA )
        return self

    def fit_tle( self ):
        # -----------------------------  nelder mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        # TODO : termination conditions should be set AND we should weight according to obs calibration
        # for now, we're using fatol as the terminating condition (with a huge xatol).  
        #FATOL = np.sqrt(30 / 3600 * self.obs_df.shape[0] )
        FATOL =  1
        ans   = scipy.optimize.minimize(optFunction, 
                                        self.get_init_fields(),
                                        args    = (self,True),
                                        method  = 'Nelder-Mead' ,
                                        options = {'xatol' : FATOL, 'fatol' : FATOL, 'initial_simplex' : self.initial_simplex() } )

        # -----------------------------  L-BFGS mead -----------------------------
        # if your seed is not near the final, nelder works great (at the expense of time)
        # TODO : termination conditions should be set AND we should weight according to obs calibration
        # for now, we're using fatol as the terminating condition (with a huge xatol).  
        #FATOL = np.sqrt(30 / 3600 * self.obs_df.shape[0] )
        #FATOL = 1
        #ans   = scipy.optimize.minimize(optFunction, 
        #                                self.get_init_fields(),
        #                                args    = (self,True),
        #                                method  = 'L-BFGS-B' ,
        #                                options = {'xatol' : FATOL, 'fatol' : FATOL, 'finite_diff_rel_step' : 1} )

        self.ans = ans
        if ans.success:
            self.reset_tle()
            # now update with perturbed values (not sure if this is necessary... last step should be there)
            for k,v in zip(self.FIELDS,ans.x) : 
                self.new_tle[k] = v
            self.final_answer = optFunction( ans.x, self, False )
            return True
        return False


# -----------------------------------------------------------------------------------------------------
def test():
    '''
    test by propagating a satellite forward, generating synthetic (perfect) obs.. and fitting those
    in essence: this is a TLE epoch mover that uses synthetic obs
    '''
    import public_astrostandards as PA
    L1       = '1 88888           25289.62065319 +.00000000  00000 0  00000 0 0 0000'
    L2       = '2 88888  56.8172  40.3738 0061389  69.9068 290.7527  2.0056030900000'
    dts      = pd.date_range( '2025-12-01', '2025-12-14', freq='5min')
    dates_df = astro_time.convert_times( dts, PA )
    print('Generating synthetic obs for:\n{}\n{}'.format( L1, L2 ) )
    print('Over {} -- {} with 5 minute spacing\n'.format( dts[0], dts[-1] ))
    #dates_df['obTime'] = dts.values
    # create a sensor frame
    sens_df  = dates_df.copy()
    # fake the sensor at the center of the earth
    sens_df['senlat'], sens_df['senlon'], sens_df['senalt'] = 0., 0., 0.
    sens_df            = sensor.prepUDLSensor( sens_df, PA )
    # propagate the TLE to the dates
    eph_df   = sgp4.propTLE_df( dates_df, L1, L2, PA )
    # generate some fake looks of that target from our fake sensor
    looks    = sensor.compute_looks( sens_df, eph_df, PA )
    # now fake a UDL set of obs
    udlobs   = pd.DataFrame()
    udlobs['obTime']        = dts
    udlobs['ra']            = looks['XA_TOPO_RA']
    udlobs['declination']   = looks['XA_TOPO_DEC']
    udlobs['senlat'], udlobs['senlon'], udlobs['senalt'] = 0., 0., 0.
    # set up your fitter 
    FIT      = eo_fitter( PA )
    FIT = FIT.set_data( L1, L2, udlobs ).set_satno(88880).set_type0()
    converged = FIT.fit_tle()
    if converged:
        rv = FIT.final_answer
        nl1, nl2 = FIT.getLines()
        print('Type 0 fit-test')
        print('\nOld TLE:\n{}\n{}'.format( L1, L2 )) 
        print('New TLE:\n{}\n{}'.format( nl1, nl2 ) )
        rv.update({ 'new_line1' : nl1, 'new_line2' : nl2 })
    else:
        print('Did not converge')

    FIT = FIT.set_data( L1, L2, udlobs ).set_satno(88881).set_type4()
    converged = FIT.fit_tle()
    if converged:
        rv = FIT.final_answer
        nl1, nl2 = FIT.getLines()
        print('Type 4 fit-test')
        print('\nOld TLE:\n{}\n{}'.format( L1, L2 )) 
        print('New TLE:\n{}\n{}'.format( nl1, nl2 ) )
        rv.update({ 'new_line1' : nl1, 'new_line2' : nl2 })
    else:
        print('Did not converge')



# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import sys
    import pandas as pd
    import public_astrostandards as PA
    import argparse
    import json
    from . import astro_time
    from . import utils

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

    parser.add_argument('--satno',  "-N",
            required = False, 
            default  = 99999,
            type     = int,
            help     = 'new TLE satno')
    
    

    # parse the arguments
    args = parser.parse_args()

    # init the public_astrostandards harness 
    PA.init_all()
    #astro_time.load_time_constants('./reduced_time_constants.dat', PA )
    astro_time.load_time_constants( utils.get_test_time_constants(), PA )

    # init the fitter
    FIT = eo_fitter( PA ) 
    
    # check on type (this should be called before set_data to get Kozai/Brouwer right)
    if args.type not in set([0,2,4]) :
        print('Valid TLE types are 0,2,4.  You input {}'.format( args.type ) )
        sys.exit(1)
    
    # load up the obs
    obs    = pd.read_json( args.infile ).sort_values(by='obTime').reset_index(drop=True)    

    # add in the data
    FIT = FIT.set_data( args.line1, args.line2, obs ).set_satno(args.satno)

    if args.type == 0 : 
        FIT.set_type0()
    if args.type == 2 :
        FIT.set_type2()
    if args.type == 4 : 
        FIT.set_type4()

    # always set your non-conservatives last
    if args.type == 4:
        FIT = FIT.set_AGOM(0.01)

    if args.verbose:
        print('Fitting')
        print('\ninitial TLE:\n\n{}\n{}'.format( FIT.line1, FIT.line2 ) )
        print('\tnew epoch : {}'.format( FIT.epoch_dt ) )
        print('\tfitting over {} obs'.format( len(obs) ) )
        print('\tspan {} -- {}'.format( FIT.obs_df.iloc[0]['obTime_dt'],
                                         FIT.obs_df.iloc[-1]['obTime_dt'] ) ) 

    # do the solve
    converged = FIT.fit_tle()
    if converged:
        if args.verbose: 
            rv = FIT.final_answer
        else: 
            rv = {}
        nl1, nl2 = FIT.getLines()
        print('\n\nNew TLE:\n{}\n{}'.format( nl1, nl2 ) )
        rv.update({ 'new_line1' : nl1, 'new_line2' : nl2 })
        with open( args.outfile, 'wt') as F: 
            json.dump( rv, F, default=str )
    else:
        print('Did not converge')

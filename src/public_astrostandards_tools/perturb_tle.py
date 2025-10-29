import numpy as np
import pandas as pd
from . import astro_time
from . import sgp4
from . import orbit_utils
from . import tle_fitter
from . import ephem_fitter

# =====================================================================================================
# Suppose we want to perturb a TLE at epoch.  How can we do that?  It is very nonlinear to encode a 
# perturbation into a TLE.  We'll do this local linearization:
#   - get the osculating elements from a TLE at epoch
#   - apply the perturbation(s) to the epoch state vector
#   - find the deltas in the osculating elements
#   - apply those deltas to the original osculating elements
#
# KN Wood (2025/10/29)
# =====================================================================================================

# -----------------------------------------------------------------------------------------------------
def subtract_dicts( A, B ):
    # routine to subtract one dict from another (by key) 
    # returns a dict of A - B
    rv = {}
    for k in A:  rv[k] = A[k] - B[k]
    return rv

# -----------------------------------------------------------------------------------------------------
def get_perturbed_XA_TLE( matrix, harness ):
    '''
    matrix is a Mx6 matrix of M state vectors (to find the perturbation from)
    '''
    XA_KEP  = [orbit_utils.sv_to_osc( X, harness).toDict() for X in matrix ]
    XA_delt = [subtract_dicts( X, epochKep.toDict() ) for X in XA_KEP ]
    return [ self.perturb_XA_TLE( self.XA_TLE, X ) for X in XA_delt ]

# -----------------------------------------------------------------------------------------------------
def perturb_XA_TLE( XA_TLE_original, XA_KEP_perturb, harness, satno=99999 ):
    TSTR = harness.Cstr('',512)
    XA_TLE_new = harness.helpers.astrostd_named_fields( harness.TleDll, prefix='XA_TLE_' )
    for k,v in XA_TLE_original.toDict().items(): XA_TLE_new[k] = v
    XA_TLE_new[ 'XA_TLE_ECCEN' ]  += XA_KEP_perturb['XA_KEP_E']
    XA_TLE_new[ 'XA_TLE_ECCEN' ]  = np.abs( XA_TLE_new[ 'XA_TLE_ECCEN' ] )
    XA_TLE_new[ 'XA_TLE_INCLI' ]  += XA_KEP_perturb['XA_KEP_INCLI']
    XA_TLE_new[ 'XA_TLE_MNANOM' ] += XA_KEP_perturb['XA_KEP_MA']
    XA_TLE_new[ 'XA_TLE_NODE' ]   += XA_KEP_perturb['XA_KEP_NODE']
    XA_TLE_new[ 'XA_TLE_OMEGA' ]  += XA_KEP_perturb['XA_KEP_OMEGA']

    # find the new mean motion (rounadbout way)
    old_A  = harness.AstroFuncDll.NToA( XA_TLE_original['XA_TLE_MNMOTN'] )
    new_A  = old_A + XA_KEP_perturb['XA_KEP_A']
    XA_TLE_new['XA_TLE_MNMOTN'] = harness.AstroFuncDll.AToN( new_A )
    return XA_TLE_new

# -----------------------------------------------------------------------------------------------------
def perturbTLE( EF : ephem_fitter.ephem_fitter,
                RIC : list[ float, float, float],
                satnos : list[ int ] = None):
    '''
    assume that you get an ephem_fitter (or tle_fitter) as an input,
    we'll use this to give back a perturbed set of TLE
    '''
    # did we get satno's to assign
    if satnos == None:
        satnos = range(90000,90000+len(RIC))
    assert len(satnos) == len(RIC)

    # get the state vector at epoch
    date_df = astro_time.convert_times( [ EF.getOriginalEpoch() ] , EF.PA )
    L1, L2  = EF.getOriginalLines()
    sv_df   = sgp4.propTLE_df( date_df , L1, L2, EF.PA ).iloc[0]

    # get the initial keplerian elements at epoch
    iXA_KEP = orbit_utils.sv_to_osc( sv_df, EF.PA )

    # find the RIC frame there
    R,I,C   = orbit_utils.getRIC( sv_df )
    P,V     = sv_df['teme_p'], sv_df['teme_v']

    # build the perturbed state vectors
    svs     = []
    for X in RIC:
        dR = R * X[0]
        dI = I * X[1]
        dC = C * X[2]
        newV = V + dR + dI + dC
        svs.append( {'teme_p' : P, 'teme_v' : newV} )

    # turn these states into osculating; find the delta to the initial sv
    XA_KEP  = [orbit_utils.sv_to_osc( X, EF.PA ).toDict() for X in svs]
    XA_delt = [subtract_dicts( X, iXA_KEP ) for X in XA_KEP ]

    # apply those delta to the original TLE and get the new result
    new_XA = [ perturb_XA_TLE( EF.init_tle, X, EF.PA ) for X in XA_delt ]
    # and turn these back into strings
    str_XA = [orbit_utils.XA_TLE_to_str(X,EF.PA,satno=Y) for X,Y in zip(new_XA,satnos) ]
    for X in str_XA: print('\n'.join(X))
    return str_XA

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    PA.init_all()

    print('-'*100)
    print('Performing fit test')
    print()
    # example TLE 
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    # this is your fit range
    #DATES = pd.date_range( '2025-1-7', '2025-1-9', freq='5min' )
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )

    # setup the job
    EH = ephem_fitter.ephem_fitter( PA ).set_from_tle(L1, L2, DATES ).set_satno(77777).set_type0()
    perturbTLE(EH,[ [.1,.1,.1], [.1,0,0], [0,.1,0], [0,0,.1] ] )
    print(EH.getOriginalEpoch())

    print('-'*100)
    print('Fitting :\n\t{}\n\t{}'.format( L1, L2 ) )
    print('\t{} -- {}'.format( DATES[0], DATES[-1] ) )
    print('\t{} points'.format( len(DATES) ) )
    print('-'*100)
    print()

    output = EH.fit_tle( )

    print()
    print('Your original TLE was')
    print('\n'.join( [L1,L2] ) ) 


    print('Your new TLE is :')
    print('\n'.join( output.getLines() ) )
    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    test()

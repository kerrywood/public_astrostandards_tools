import numpy as np
import pandas as pd
from . import astro_time
from . import sgp4
from . import orbit_utils
from . import tle_fitter

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
    for k in A:  
        rv[k] = A[k] - B[k]
    return rv

# # -----------------------------------------------------------------------------------------------------
# def get_perturbed_XA_TLE( matrix, harness ):
#     '''
#     matrix is a Mx6 matrix of M state vectors (to find the perturbation from)
#     '''
#     XA_KEP  = [orbit_utils.sv_to_osc( X, harness).toDict() for X in matrix ]
#     XA_delt = [subtract_dicts( X, XA_KEP.toDict() ) for X in XA_KEP ]
#     return [ perturb_XA_TLE( XA_TLE, X ) for X in XA_delt ]

# -----------------------------------------------------------------------------------------------------
def perturb_XA_TLE( XA_TLE_original, XA_KEP_perturb, harness, satno=99999 ):
    # TSTR = harness.Cstr('',512)
    XA_TLE_new = harness.helpers.astrostd_named_fields( harness.TleDll, prefix='XA_TLE_' )
    for k,v in XA_TLE_original.toDict().items(): 
        XA_TLE_new[k] = v
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
def anomaly_TLE(TF : tle_fitter.tle_fitter,
                samples : int = 20 ,
                satnos : list[ int ] = None):
    '''
    break up the anomalies and develop TLE for that
    '''
    # did we get satno's to assign
    if satnos is None:
        satnos = range(90000,90000+samples)
    if isinstance(satnos, int):
        satnos = range( satnos, satnos + samples )
    assert len(satnos) == samples

    # our new holder; copy over all the original data
    XA_TLE_new = TF.PA.helpers.astrostd_named_fields( TF.PA.TleDll, prefix='XA_TLE_' )
    for i,d in enumerate(TF.init_tle.data): 
        XA_TLE_new.data[i] = TF.init_tle.data[i]

    # linspace 
    new_anomaly  = ( TF.init_tle['XA_TLE_MNANOM'] + np.linspace(0,360,samples) ) % 360

    # generate the perturbed / shifted TLE
    def genNewTLE( anom, sno ):
        XA_TLE_new['XA_TLE_MNANOM'] = anom
        return orbit_utils.XA_TLE_to_str( XA_TLE_new, TF.PA, satno = sno )

    str_XA = [ genNewTLE( X, Y ) for X,Y in zip( new_anomaly, satnos ) ]
    return str_XA

# -----------------------------------------------------------------------------------------------------
def perturbTLE( EF : tle_fitter.tle_fitter,
                RIC : list[ float, float, float],
                satnos : list[ int ] = None):
    '''
    assume that you get an test_fitter as an input,
    we'll use this to give back a perturbed set of TLE
    '''
    # did we get satno's to assign
    if satnos is None:
        satnos = range(90000,90000+len(RIC))
    if isinstance(satnos, int):
        satnos = range( satnos, satnos + len(RIC) )
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
    return str_XA
    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import sys
    import public_astrostandards as PA
    import argparse
    import json
    from . import test_helpers

    parser = argparse.ArgumentParser(
        prog='TLE perturb-er',
        description='''given a TLE, do an osculating approximation around epoch for your chosen perturbation''',
        epilog = ''
    )

    parser.add_argument('--line1',"-l1",
            required = True,
            help     = 'line 1 of a TLE')

    parser.add_argument('--line2','-l2',
            required = True,
            help     = 'line2 of TLE')

    parser.add_argument('--RIC', '-r',
            required = False,
            help     = 'RIC perturbations (list of lists; km/s ; JSON parse-able) : [[1,0,0],[.1,.1,.1]]' )

    parser.add_argument('--satno',  "-N",
            required = False, 
            type     = int,
            help     = 'new TLE satno')
    
    
    # parse the arguments
    args = parser.parse_args()

    # init the public_astrostandards harness 
    PA.init_all()
    astro_time.load_time_constants( test_helpers.get_test_time_constants(), PA )

    # init the fitter
    TF = tle_fitter.tle_fitter( PA ).set_from_lines(args.line1, args.line2 ).set_satno(77777).set_type0()

    # -----------------------------------------------------------------------------------------------------
    # if this is a RIC call
    if args.RIC:
        ric = json.loads( args.RIC )
        if args.satno : newtle = perturbTLE( TF, ric, args.satno )
        else : newtle = perturbTLE( TF, ric )
        print()
        print('{}\n{}'.format( args.line1, args.line2 ) )
        for t in newtle: print('\n'.join(t))
        sys.exit()

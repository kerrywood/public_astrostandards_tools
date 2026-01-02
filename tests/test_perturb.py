import pandas as pd
import public_astrostandards as PA
import public_astrostandards_tools as PAT

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    PA.init_all()

    # example TLE 
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    RIC_vecs = [ 
                [.1,.1,.1], # 100 m/s in all RIC directions
                [.1,0.,0.],   # 100 m/s radial
                [0.,.1,0.],   # 100 m/s in-track
                [0.,0.,.1],
                [0.,1.,0.],
                [0.,1.,0.],
                [1.,0.,0.],
                [0.,0.,1.],
                [0.,2.,2.] ]  # 2 km/s in-track and cross

    # setup the job
    EH = PAT.tle_fitter.tle_fitter( PA ).set_from_lines(L1, L2 ).set_satno(77777).set_type0()
    
    print('-'*100)
    print('Perturbing:\n{}\n{}'.format( L1,L2 ) )
    print()
    print('Epoch : {}'.format(EH.getOriginalEpoch()) )
    print('RIC perturbations (km/s)')
    for V in RIC_vecs:
        print('{:8.3f} {:8.3f} {:8.3f}'.format( *V ) )
    print('TLE')
    
    PAT.perturb_tle.anomaly_TLE( EH,20 )


    perturbed_set = PAT.perturb_tle.perturbTLE(EH, RIC_vecs )


    for tle in perturbed_set:
        print('\n'.join( tle ))

# =====================================================================================================
if __name__ == "__main__":
    test()

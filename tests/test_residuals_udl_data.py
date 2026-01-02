# ====================================================================================================
# use with output from `UDLOB_to_B3_ROTAS` output from `weekly_dnd` code base
# that code will use the ITAR astrostandards and ROTAS to get output; we can compare our tools here
# ====================================================================================================

import pandas as pd
import public_astrostandards as PA
import public_astrostandards_tools as PAT

# -----------------------------------------------------------------------------------------------------
def test():
    PA.init_all()
    # load the stored time constants file 
    PAT.astro_time.load_time_constants( PAT.test_helpers.get_test_time_constants(), PA )

    # test data
    test_data = pd.read_csv('~/Downloads/rotas_output.csv')

    # pull out the TLE (though each line has data, there should be one TLE)
    TLE1 = test_data['tle_line1'].unique()[0]
    TLE2 = test_data['tle_line2'].unique()[0]
    
    # prep that for comparison
    obs_df    = PAT.observations.prepUDLObs( test_data.copy(), PA )

    # calculate residuals using open source code
    output = PAT.residuals.UDL_ROTAS( obs_df, TLE1, TLE2, PA )

    # /////////////////////////////////// summarize output \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # field pairings
    pairs = ( 
             ( 'residual_ra', 'XA_OBSRES_RA'),
             ( 'residual_dec', 'XA_OBSRES_DEC'),
             ( 'residual_az', 'XA_OBSRES_AZ'),
             ( 'residual_el', 'XA_OBSRES_EL'),
             ( 'residual_range', 'XA_OBSRES_RANGE' )
    )
    
    compare = pd.DataFrame()
    for A,B in pairs:
        if 'range' in A: 
            compare['{}-{}'.format(A, B) ] = output[A] - test_data[B]
        else:
            compare['{}-{}'.format(A, B) ] = 3600 * (output[A] - test_data[B])
                     
    print(compare.describe())

# =====================================================================================================
if __name__ == "__main__":
    test()
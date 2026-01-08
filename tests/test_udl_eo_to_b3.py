# ====================================================================================================
# download some EO obs from UDL; test conversion
# ====================================================================================================
import pandas as pd
import public_astrostandards as PA
import public_astrostandards_tools as PAT

# -----------------------------------------------------------------------------------------------------
def test():
    PA.init_all()
    # load the stored time constants file 
    PAT.astro_time.load_time_constants( PAT.utils.get_test_time_constants(), PA )

    # load the test obs (direct downloaded from UDL)
    test_obs = pd.read_json('~/Downloads/27566.json')
    print( 'I see {} obs'.format( test_obs.shape[0] ))

    # do the conversion
    test_obs  = PAT.observations.UDLEOObstoB3Type9( test_obs, PA )

    print( test_obs['B3'] )

# =====================================================================================================
if __name__ == "__main__":
    test()
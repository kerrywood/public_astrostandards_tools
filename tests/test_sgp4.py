import pandas as pd

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    import public_astrostandards_tools as PAT

    # SETUP info
    ISS = ('1 25544U 98067A   25357.18166772  .00011641  00000-0  21351-3 0  9998','2 25544  51.6323  90.7678 0003190 289.6661  70.3984 15.49746572544475')
    dates = pd.date_range( '2025-12-23', '2026-1-15',  freq='1min' )
    sen_lla = (38.83, -104.82, 1.832 )
    
    print('*' * 100)
    print('TLE:')
    print('\t{}\n\t{}'.format( *ISS ) )
    print('Date range:')
    print('\t{}\n\t{}'.format( dates[0], dates[-1] ))
    print('*' * 100)

    # init all the Dll's
    # PA.init_all()

    # use the TimeFunc to load the time parameters file (need to upate this periodically)
    # PAT.astro_time.load_time_constants(  PAT.utils.get_test_time_constants(), PA )

    # use the astro_time to initialize the dataframe with times
    # note that the sensor and target dataframes must be time aligned
    dates_f = PAT.astro_time.convert_times( dates, PA )
    # make a target dates from that is identical
    target_f = dates_f.copy()

    # manually set the license path (oddly, setting a bad path seems to work)
    PAT.sgp4.setLicensePath( '', PA ) 
    print('License path: {}'.format( PAT.sgp4.getLicensePath( PA ) ) )

    target_f = PAT.sgp4.propTLE_df( dates_f.copy(), *ISS, PA ) 

    
# =====================================================================================================
if __name__ == "__main__":
    test()

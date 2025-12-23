import pandas as pd

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    import public_astrostandards_tools as PAT
    # init all the Dll's
    PA.init_all()

    # use the TimeFunc to load the time parameters file (need to upate this periodically)
    PAT.astro_time.load_time_constants(  PAT.test_helpers.get_test_time_constants(), PA )

    # generate some test data
    ISS = ('1 25544U 98067A   25357.18166772  .00011641  00000-0  21351-3 0  9998','2 25544  51.6323  90.7678 0003190 289.6661  70.3984 15.49746572544475')
    dates = pd.date_range( '2025-12-23', '2026-1-15',  freq='1min' )

    # use the astro_time to initialize the dataframe with times
    # note that the sensor and target dataframes must be time aligned
    dates_f = PAT.astro_time.convert_times( dates, PA )
    # make a target dates from that is identical
    target_f = dates_f.copy()

    # -----------------------------------------------------------------------------------------------------
    # TEST case 1 : look vectors to sun; find out when sun is down
    # make a sensor frame
    sensor_f = PAT.sensor.setup_ground_site( dates_f.copy(), 38.83, -104.82, 1.832, PA )
    # annotate that dataframe wih the sun position (why not make this a routine?  Because we can re-use sun position at these times.)
    target_f['teme_p']  = PAT.sensor.sun_at_time( dates_f, PA ) 
    # compute looks to the sun
    looks_f = PAT.sensor.compute_looks( sensor_f, target_f, PA )
    # find those times when the sun is down; NOTE we're indexing subsets of the DATE and SENSOR frame
    sensor_sundown_f = sensor_f[ looks_f['XA_TOPO_EL'] < -4 ].copy()
    sensor_sundown_dates = sensor_sundown_f[ PAT.astro_time.DATE_FIELDS ].copy()

    # -----------------------------------------------------------------------------------------------------
    # TEST case 2 : now that we know when sun is down, can we find those times when ISS is visible?
    # for now, ignore if it is solar illuminated
    # re-use the dates  
    # >> note that we're only checking at those times that the sun is down << 
    # >> this limits the propagation calls (more efficient) <<
    target_f = PAT.sgp4.propTLE_df( sensor_sundown_dates, *ISS, PA ) 
    # check if the target is sunlit... (annotate each row)
    target_f['is_sunlit'] = PAT.sensor.is_sunlit( target_f, PA )
    # now compute the actual looks....
    looks_f = PAT.sensor.compute_looks( sensor_sundown_f, target_f, PA )

    # note that here, _target is appended to a field we pushed into compute_looks
    # so we look for _is_sunlit_target
    good = looks_f[ (looks_f['XA_TOPO_EL'] > 5) * (looks_f['is_sunlit_target'] == 1 ) ] 
    pd.set_option('display.max_rows', None)
    print( good[['datetime_sensor','XA_TOPO_EL','XA_TOPO_AZ','XA_TOPO_RANGE','is_sunlit_target']] )

    
# =====================================================================================================
if __name__ == "__main__":
    test()

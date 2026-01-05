import pandas as pd
import public_astrostandards as PA
import public_astrostandards_tools as PAT


def test() :
    # init the astrostandards
    PA.init_all()
    # use the TimeFunc to load the time parameters file (need to upate this periodically)
    PAT.astro_time.load_time_constants(  PAT.utils.get_test_time_constants(), PA )

    # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # some input TLE ; hand-modify a TLE as our "error orbit" (we'll use an error orbit to fake obs with noise)
    # in this case, we'll just truncate mean motion and inclination (look at those fields) and pretend
    # that this is the "observed" orbit, where "actual" is our hypothesis orbit
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    ISS = ('1 25544U 98067A   25357.18166772  .00011641  00000-0  21351-3 0  9998','2 25544  51.6323  90.7678 0003190 289.6661  70.3984 15.49746572544475')
    TDRS_actual = ('1 27566U 02055A   25357.24095851  .00000061  00000-0  00000-0 0  9991','2 27566   9.7383  44.3591 0016647 235.9259 132.2209  0.98860736 84461')
    TDRS_mod    = ('1 27566U 02055A   25357.24090000  .00000061  00000-0  00000-0 0  9991','2 27566   9.7383  44.3591 0016647 235.9259 132.2209  0.98800000 84461')

    # generate some test dates
    test_dates = pd.date_range('2025-12-23','2025-12-30',freq='5 min')
    # format the dates for use with the other tools
    test_dates = PAT.astro_time.convert_times( test_dates, PA )

    # ---------------------------------------------------------------------------------------
    # we'll build some synthetic observations and convert them to UDL format
    # generate the ephemeris
    ISS_frame = PAT.sgp4.propTLE_df( test_dates.copy(), *ISS, PA )
    TDRS_actual_ephem = PAT.sgp4.propTLE_df( test_dates.copy(), *TDRS_actual, PA )
    TDRS_mod_ephem = PAT.sgp4.propTLE_df( test_dates.copy(), *TDRS_mod, PA )

    # to use the sensor.compute_looks, we need to add in LLH data... (required for looks)
    ISS_frame = PAT.coordinates.TEME_to_LLH( ISS_frame, PA )
    # generate the look vectors to the modified TLE
    looks = PAT.sensor.compute_looks( ISS_frame, TDRS_mod_ephem, PA)

    # ////////////////////////////////////////////////////////
    # we have to fake these *back* into UDL format....
    PAT.observations.synthetic_to_UDL_like( ISS_frame, looks )

    # now that we have faked UDL observations, feed these to the ROTAS calculator against the "real" orbit
    residuals = PAT.residuals.UDL_ROTAS( ISS_frame, *TDRS_actual, PA)
    residuals['ra_arcsec'] = residuals['residual_ra'] * 3600
    residuals['dec_arcsec'] = residuals['residual_dec'] * 3600
    print(residuals)

    # print(looks.columns)
    slatton = PAT.residuals.slatton_intersection( looks, ISS_frame) 
    print( slatton )

# =====================================================================================================
if __name__ == "__main__":
    test()

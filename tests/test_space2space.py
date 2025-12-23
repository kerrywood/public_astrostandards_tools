# #######################################################################################
# Generate a synthetic space-to-space obs scenario.
#   in this case, assume ISS is looking at TDRS 
#
# Kerry N. Wood (kerry.wood@asterism.ai)
# #######################################################################################

import numpy as np
import pandas as pd
import public_astrostandards as PA
import public_astrostandards_tools as PAT

# init the astrostandards
PA.init_all()
# use the TimeFunc to load the time parameters file (need to upate this periodically)
PAT.astro_time.load_time_constants(  PAT.test_helpers.get_test_time_constants(), PA )

ISS = ('1 25544U 98067A   25357.18166772  .00011641  00000-0  21351-3 0  9998','2 25544  51.6323  90.7678 0003190 289.6661  70.3984 15.49746572544475')
TDR = ('1 27566U 02055A   25357.24095851  .00000061  00000-0  00000-0 0  9991','2 27566   9.7383  44.3591 0016647 235.9259 132.2209  0.98860736 84461')

# generate some test dates
test_dates = pd.date_range('2025-12-23','2025-12-30',freq='5 min')

# format the dates for use with the other tools
test_dates = PAT.astro_time.convert_times( test_dates, PA )

# ---------------------------------------------------------------------------------------
# we'll build some synthetic observations and convert them to UDL format
# generate the ephemeris
ISS_ephem = PAT.sgp4.propTLE_df( test_dates.copy(), *ISS, PA )
TDR_ephem = PAT.sgp4.propTLE_df( test_dates.copy(), *TDR, PA )

# to use the sensor.compute_looks, we need to add in LLH data... (required for looks)
ISS_ephem = PAT.coordinates.TEME_to_LLH( ISS_ephem, PA )

# generate the look vectors
looks = PAT.sensor.compute_looks( ISS_ephem, TDR_ephem, PA)

# we have to fake these *back* into UDL format....
ISS_ephem['teme_ra']  = looks['XA_TOPO_RA']
ISS_ephem['teme_dec'] = looks['XA_TOPO_DEC']
ISS_ephem['senlat']   = ISS_ephem['lat']
ISS_ephem['senlon']   = ISS_ephem['lon']
ISS_ephem['senalt']   = ISS_ephem['height']
ISS_ephem['teme_lv']  = PAT.observations.ra_dec_to_lv( ISS_ephem['teme_ra'], ISS_ephem['teme_dec'] ).tolist()

residuals = PAT.residuals.UDL_ROTAS( ISS_ephem, *TDR, PA)
residuals['ra_arcsec'] = residuals['ra'] * 3600
residuals['dec_arcsec'] = residuals['dec'] * 3600
print(residuals)
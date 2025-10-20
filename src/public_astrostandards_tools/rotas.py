from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import scipy.optimize

# -----------------------------------------------------------------------------------------------------
def getUVW( sv_df : pd.DataFrame ):
    sv_df['U_c'] = sv_df['teme_p'].apply( lambda X : X / np.linalg.norm(X))

    def getW( R ):
        tv = np.cross( R['teme_p'], R['teme_v'] ) 
        return tv / np.linalg.norm(tv)

    sv_df['W_c'] = sv_df.apply( getW, axis=1 )

    def getV( R ):
        return np.cross( R['W_c'], R['U_c'] )

    sv_df['V_c'] = sv_df.apply( getV, axis=1 )
    return sv_df 


# =====================================================================================================
if __name__ == '__main__':
    import astro_time
    import observations
    import sgp4
    import sensor
    import utils

    import public_astrostandards as PA
    PA.init_all()

    # TLE
    with open('../../tests/27566.tle') as F: lines = F.readlines()
    L1 = lines[0].strip()
    L2 = lines[1].strip()

    # obs 
    obs_df    = pd.read_json('../../tests/27566.json.gz').sort_values(by='obTime')
    # reformat UDL obs (these are our actual TEST obs ; from a sensor)
    obs_df    = observations.prepUDLObs( obs_df, PA )
    # now that those have the correct date fields, peel the dates out of the obs
    date_df   = obs_df[ astro_time.DATE_FIELDS ].copy()
    # get a sensor frame (from the OBS again; we're using those ground sensors)
    sensor_df = sensor.prepUDLSensor( obs_df, PA )
    # get hypothesis obs.. (from a TLE)
    eph_df    = sgp4.propTLE_df( date_df, L1, L2, PA )
    # turn each P,V into osculating elements (for ROTAS comparison)
    eph_df    = utils.sv_to_osc_df( eph_df, PA )
    # find the U,V,W frame
    eph_df    = getUVW(eph_df)

    # now, generate hypothesis looks
    looks_df  = sensor.compute_looks( sensor_df, eph_df, PA )

    residuals_df = observations.residuals( obs_df, looks_df )
    print(residuals_df)


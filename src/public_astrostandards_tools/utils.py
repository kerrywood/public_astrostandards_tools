import ctypes 
import pandas as pd

#,  -----------------------------------------------------------------------------------------------------
def sv_to_osc( sv, PA ):
    '''
    sv : <teme_pos><teme_vel>
    return XA_KEP
    '''
    # we'll use the conversion in the astrostandards
    XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.PosVelToKep( 
        (ctypes.c_double*3)(*sv['teme_p']), 
        (ctypes.c_double*3)(*sv['teme_v']), 
        XA_KEP.data )
    return XA_KEP

#,  -----------------------------------------------------------------------------------------------------
def sv_to_osc_df( sv_df : pd.DataFrame, PA ) :
    tv = sv_df.apply( lambda X: sv_to_osc(X,PA).toDict(), axis=1 )
    tv = pd.DataFrame.from_dict( tv.values.tolist() )
    return pd.concat( (sv_df.reset_index(drop=True), tv.reset_index(drop=True) ), axis =1 )


# -----------------------------------------------------------------------------------------------------
def osc_to_mean( XA_KEP, PA ):
    XA_KEP_MEAN = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.KepOscToMean( XA_KEP.data, XA_KEP_MEAN.data )
    return XA_KEP_MEAN

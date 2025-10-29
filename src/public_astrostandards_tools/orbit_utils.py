import ctypes 
import pandas as pd

#,  -----------------------------------------------------------------------------------------------------
def sv_to_osc( sv, PA ):
    '''
    sv : <teme_pos><teme_vel>
    return XA_KEP
    '''
    XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    # we'll use the conversion in the astrostandards
    PA.AstroFuncDll.PosVelToKep( 
        (ctypes.c_double*3)(*sv['teme_p']), 
        (ctypes.c_double*3)(*sv['teme_v']), 
        XA_KEP.data )
    return XA_KEP

#,  -----------------------------------------------------------------------------------------------------
def sv_to_osc_df( sv_df : pd.DataFrame, PA ) :
    ''' 
    given a dataframe with 'teme_p' and 'teme_v' on each row, annotate each row with XA_KEP data
    '''
    tv = sv_df.apply( lambda X: sv_to_osc(X,PA), axis=1 )
    tv = [ X.toDict() for X in tv ]
    tv = pd.DataFrame.from_dict( tv )
    return pd.concat( (sv_df.reset_index(drop=True), tv.reset_index(drop=True) ), axis =1 )

# -----------------------------------------------------------------------------------------------------
def osc_to_mean( XA_KEP, PA ):
    '''
    take a XA_KEP structure and return a XA_KEP with *mean* fields
    '''
    XA_KEP_MEAN = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
    PA.AstroFuncDll.KepOscToMean( XA_KEP.data, XA_KEP_MEAN.data )
    return XA_KEP_MEAN

from datetime import datetime, timedelta, timezone
import numpy as np
import pandas as pd
from . import astro_time
   
# -----------------------------------------------------------------------------------------------------
def toCCSDS( df : pd.DataFrame ):
    '''
    assume that the input dataframe has
        datetime :   (datetime or datetime-like)
        j2k_p    :   J2K position vector
        j2k_v    :   J2K velocity vector

        generate an output that is CCSDS-like
    '''
    odf = pd.DataFrame()
    odf['datetime']  = df['datetime'].apply( lambda X: X.strftime('%Y-%m-%dT%H:%M:%S.%f') )
    odf['j2k_x']     = df['j2k_p'].apply( lambda X: X[0] )
    odf['j2k_y']     = df['j2k_p'].apply( lambda X: X[1] )
    odf['j2k_z']     = df['j2k_p'].apply( lambda X: X[2] )
    odf['j2k_dx']    = df['j2k_v'].apply( lambda X: X[0] )
    odf['j2k_dy']    = df['j2k_v'].apply( lambda X: X[1] )
    odf['j2k_dz']    = df['j2k_v'].apply( lambda X: X[2] )
    return odf.to_csv( index=None , sep='\t', header=None )

# -----------------------------------------------------------------------------------------------------
def fromCCSDS( lines, harness ): 
    '''
    lines : file data

    return an annotated pandas dataframe (in standard coordinate format with j2k_p, and j2k_v headers
    '''
    def parseLine(L):
        try:
            flds = L.strip().split('\t')
            assert len(flds) == 7 
            dt   = datetime.strptime(flds[0], '%Y-%m-%dT%H:%M:%S.%f' )
            return (dt, [ float(X) for X in flds[1:4] ], [float(X) for X in flds[4:] ] )
        except Exception as e: 
            print(e)
            return None
    prsd  = [ parseLine(X) for X in lines ]
    lines = list( filter( lambda X: X is not None, prsd ) )
    tv    = pd.DataFrame( lines, columns=['datetime','j2k_p','j2k_v'] )
    dt_df = astro_time.convert_times( tv['datetime'], harness ).drop( columns='datetime' )
    return pd.concat( (dt_df.reset_index(drop=True), tv.reset_index(drop=True) ), axis=1 )

# -----------------------------------------------------------------------------------------------------
def test():
    import public_astrostandards as PA
    from . import astro_time
    from . import sgp4
    from . import utils
    from . import coordinates

    # init the astrostandards 
    PA.init_all()
    astro_time.load_time_constants(  utils.get_test_time_constants(), PA )

    # test TLE
    L1 = '1 25544U 98067A   25301.52216109  .00016210  00000-0  29595-3 0  9990'
    L2 = '2 25544  51.6346   6.3420 0004740 349.4592  10.6297 15.49553329535849'
    
    # generate some test data
    print('*'*100)
    print('Propagating test TLE and converting to J2K')
    print('\t{}\n\t{}'.format( L1, L2 ) )
    print('*'*100)
    print()
    DATES      = pd.date_range('2025-11-1', '2025-11-2', freq='5min')
    dates_f    = astro_time.convert_times( DATES, PA )
    original_f = sgp4.propTLE_df( dates_f, L1, L2, PA )
    to_j2k     = coordinates.TEME_to_J2K( original_f.copy(), PA ) 

    X  = toCCSDS( to_j2k )
    print('*'*100)
    print('CCSDS output')
    print('*'*100)
    print(X)

    print('*'*100)
    print('PANDAS from CCSDS output')
    print('*'*100)
    Y  = fromCCSDS( X.split('\n'), PA )
    print(Y)
    print(Y.iloc[0])


# =====================================================================================================
if __name__ == '__main__':
    test()

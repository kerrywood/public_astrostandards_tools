# ###############################################################################
# MIT License

# Copyright (c) 2025 Kerry N. Wood (kerry.wood@asterism.ai)

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###############################################################################

import math
import os
from datetime import datetime, timedelta, timezone
import pandas as pd

# epochs to help us convert from and to Julian... (rather than a library)
_J2K_dt = datetime(year=2000, month=1, day=1, hour=12)
_J2K_jd = 2451545.0

DATE_FIELDS = ['datetime','theta','ds50_utc','ds50_et','ds50_ut1']

# --------------------------------------------------------------------------------------------------------
def dt2julian( DT : datetime ):
    del_d = ( DT - _J2K_dt ).total_seconds() / 86400
    return _J2K_jd + del_d

# --------------------------------------------------------------------------------------------------------
def julian2dt( jd : float ):
    return _J2K_dt + timedelta( days= jd - _J2K_jd )

# -----------------------------------------------------------------------------------------------------
def load_time_constants( filename : str, INTERFACE):
    assert os.path.isfile( filename )
    assert 0 == INTERFACE.TimeFuncDll.TimeFuncLoadFile( INTERFACE.Cstr(filename, 512 ) )

# -----------------------------------------------------------------------------------------------------
def datetime_to_ds50( dt : list[ datetime ] , INTERFACE ):
    return [ INTERFACE.helpers.datetime_to_ds50( X, INTERFACE.TimeFuncDll ) for X in dt ]

# -----------------------------------------------------------------------------------------------------
def convert_times( datetimes : list[ datetime ] ,
                   INTERFACE ):
    # convert the datetimes to astrostandard epochs
    def safeConvert( dt ):
        try :
            return INTERFACE.helpers.datetime_to_ds50( dt, INTERFACE.TimeFuncDll )
        except Exception as e: 
            print(e)
            return -1 
        
    ds50_utc = [ safeConvert(X) for X in datetimes ]

    # Convert ds50UTC to ds50UT1
    ds50_ut1 = [INTERFACE.TimeFuncDll.UTCToUT1(X) for X in ds50_utc]

    return pd.DataFrame( 
           [ { 'datetime' : X,
               'theta'    : INTERFACE.TimeFuncDll.ThetaGrnwchFK5( Z ),
               'ds50_utc' : Y,
               'ds50_et'  : INTERFACE.TimeFuncDll.UTCToET( Y ),
               'ds50_ut1' : Z } for X,Y,Z in zip(datetimes,ds50_utc,ds50_ut1) ])


# -----------------------------------------------------------------------------------------------------
def test():
    from . import test_helpers

    import public_astrostandards as harness
    # init all the Dll's
    harness.init_all()

    # use the TimeFunc to load the time parameters file (need to upate this periodically    )
    load_time_constants(  test_helpers.get_test_time_constants(), harness )

    # generate some test data
    dates = pd.date_range(  datetime.now( timezone.utc ), 
                            datetime.now( timezone.utc ) + timedelta(days=1), 
                            freq='5 min' )
    print(convert_times( dates, harness) )


# ==================================================================================================
if __name__ == '__main__':
    test()  

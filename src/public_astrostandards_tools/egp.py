import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import astro_time
import ephem_fitter

# ----------------------------------------------------------------------------------------------------- 
def test():
    import public_astrostandards as PA
    print()
    print('-'*100)
    print('Running tests')
    print('-'*100)
    PA.init_all()
    # example TLE 
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    # this is your fit range
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )

    EH = ephem_fitter.ephem_fitter( PA )
    EH.set_from_tle( L1, L2, DATES ).set_type0()
    out = EH.fit_tle()

    # type 4 to type 0
    # setup the job
    print()
    print('Your original TLE was')
    print('\n'.join( [L1,L2] ) ) 
    print('Fitting {} points from {} -- {}'.format( len(DATES), DATES[0], DATES[-1]) )

    EH = ephem_fitter.ephem_fitter( PA ).set_from_tle(L1, L2, DATES ).set_satno(77777).set_type0()
    output = EH.fit_tle( )
    print()
    print('\nType 0 fit:')
    print(EH.summarize_results())

    print('\nType 2 fit:')
    EH.set_type2()
    output = EH.fit_tle( )
    print(EH.summarize_results())

    print('\nType 4 fit:')
    EH.set_type4()
    output = EH.fit_tle( )
    print(EH.summarize_results())

    

# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import sys
    import pandas as pd
    import public_astrostandards as PA
    import astro_time
    import argparse
    import json

    PA.init_all()

    parser = argparse.ArgumentParser(
        prog='egp converter',
        description='''take a type-4 TLE (or really any type) and convert it to a fit region''',
        epilog = ''
    )

    parser.add_argument('--line1',"-l1",
            required=False,
            help='line 1 of a TLE', 
            default='1 25544U 98067A   24365.67842578    .00000000  00000-0  00000-0 4  9990')

    parser.add_argument('--line2','-l2',
            required=False,
            help = 'line2 of TLE',
            default = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028' )

    parser.add_argument('--startdate',"-sd", required=False,                    help='first date (in a format that pd.date_range can understand)')
    parser.add_argument('--enddate',  "-ed", required=False,                    help='last date (in a format that pd.date_range can understand)')
    parser.add_argument('--spacing',  "-sp", required=False, default='5min',    help='ephemeris spacing (in format pd.date_range can understand)')
    parser.add_argument('--type',     "-ty", required=False, default='0',       help='TLE type (0,2,4)', type=int)
    parser.add_argument('--satno',    "-sno",required=False, default=99999,     help='new TLE satno',type=int)
    parser.add_argument('--test',     "-ts", required=False, action='store_true',help='run a test (all other arguments ignored)')

    # parse the arguments
    args = parser.parse_args()

    if args.test:
        test()
        sys.exit()

    if not args.startdate:
        print('uhoh')
    # now start adding params
    assert  args.type  in set([0,2,4])
    dates = pd.date_range( args.startdate, args.enddate, freq=args.spacing )

    EH = egp_helper( PA ).set_data( args.line1, args.line2, dates )
    if args.type == 0: EH.set_type0()
    if args.type == 2: EH.set_type2()
    if args.type == 4: EH.set_type4()
    EH.set_satno( args.satno )

    #    test()
    print()
    print()
    print( json.dumps( EH.egp_tle().summarize_results(), indent=4) )

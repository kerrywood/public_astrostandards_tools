import os
import time
from datetime import datetime, timedelta
import ctypes
import numpy as np
import astro_time
import ephem_fitter

# ----------------------------------------------------------------------------------------------------- 
def doJob( L1, L2, DATES, PA):
    EH = ephem_fitter.ephem_fitter( PA )
    EH.set_from_tle( L1, L2, DATES ).set_type0()
    out = EH.fit_tle()
    return out

# ----------------------------------------------------------------------------------------------------- 
def fileJob( fn , PA ):
    assert os.path.isfile( fn )
    jobs = pd.read_json( fn )
    jobs = list( jobs.to_dict('records') )
    for job in jobs:
        freq = job.get('spacing','5min') 
        dates = pd.date_range( job['min_date'], job['max_date'], freq=freq )
        output = doJob( job['line1'], job['line2'], dates, PA )
        if output : job['egp'] = output.summarize_results()
        else : job['egp'] = None
    return jobs
# ----------------------------------------------------------------------------------------------------- 
def test():
    import public_astrostandards as PA
    print()
    print('-'*100)
    print('Running tests')
    print('-'*100)
    PA.init_all()
    # -----------------------------------------------------------------------------------------------------
    # this is a type-4 faked by a modified from a space-track TLE
    L1 = '1 25544U 98067A   24365.67842578  .00000000  00000-0  00000-0 4  9990'
    L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'
    # this is your fit range
    DATES = pd.date_range( '2025-6-1', '2025-6-2', freq='5min' )
    # -----------------------------------------------------------------------------------------------------
    print()
    print('Your original TLE was')
    print('\n'.join( [L1,L2] ) ) 
    print('Fitting {} points from {} -- {}'.format( len(DATES), DATES[0], DATES[-1]) )
    output = doJob( L1, L2, DATES, PA )

    print()
    print('\nType 0 fit:')
    print(output.summarize_results())

    # re-use that fitter and solve for a type-2
    print('\nType 2 fit:')
    output.set_type2()
    output = output.fit_tle( )
    print(output.summarize_results())

    print('\nType 4 fit:')
    output.set_type4()
    output = output.fit_tle( )
    print(output.summarize_results())

    

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
    parser.add_argument('--infile',   "-F", required=False, help='load jobs from a file (other arguments ignored)')
    parser.add_argument('--outfile',  "-O", required=False, help='store output from loaded jobs')

    # parse the arguments
    args = parser.parse_args()

    # =============================================================
    # CASE 1 : test
    # if test is specified, run that
    if args.test:
        test()
        sys.exit()

    # =============================================================
    # CASE 2 : file
    # if file is specified, run from that
    if args.infile:
        if not args.outfile:
            print('with infile, you must specify outfile')
            sys.exit(1)
        results = fileJob( args.infile, PA )
        with open( args.outfile, 'w') as F: 
            json.dump( results, F, default=str, indent=4 )
        sys.exit()
    
    # =============================================================
    # CASE 3 : command line
    # now start adding params
    assert  args.type  in set([0,2,4])
    dates = pd.date_range( args.startdate, args.enddate, freq=args.spacing )
    EH = ephem_fitter.ephem_fitter( PA ).set_data( args.line1, args.line2, dates )
    if args.type == 0: EH.set_type0()
    if args.type == 2: EH.set_type2()
    if args.type == 4: EH.set_type4()
    EH.set_satno( args.satno )

    print()
    print()
    print( json.dumps( EH.egp_tle().summarize_results(), indent=4) )

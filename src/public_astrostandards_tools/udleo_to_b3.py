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


# =================================================================================================
if __name__ == '__main__':
    import sys
    import pandas as pd
    import argparse
    import json
    import public_astrostandards as PA
    import public_astrostandards_tools as PAT

    parser = argparse.ArgumentParser(
        prog='UDL EO obs to B3',
        description='''take a set of obs from UDL (directly) and add a B3 column''',
        epilog = ''
    )

    parser.add_argument('--infile',   "-F", 
            required = True, 
            default  = './19548.json.gz',
            help     ='load obs from this file')

    parser.add_argument('--outfile',  "-O", 
            required = True,
            default  = './out.json',
            help     = 'store output from loaded jobs')
    
    parser.add_argument('--verbose',  "-v",
            required = False, 
            default  = False,
            action   = 'store_true',
            help     = 'print debugging info')


    # parse the arguments
    args = parser.parse_args()

    # init the public_astrostandards harness 
    PA.init_all()
    
    # config time constants from the loaded file
    PAT.astro_time.load_time_constants( PAT.utils.get_test_time_constants(), PA )
    
    if args.verbose: 
        print('Loaded time constants from : {}'.format(  PAT.utils.get_test_time_constants() ))
    
    # load up the obs
    obs = pd.read_json( args.infile ).sort_values(by='obTime').reset_index(drop=True)   
    
    if args.verbose:
        print('Loaded {} obs from {}'.format( obs.shape[0], args.infile ))
    
    if args.verbose:
        print('Converting obs...')

    # do the conversion
    obs = PAT.observations.UDLEOObstoB3Type9( obs, PA ) 
    
    bad  = obs[ obs['B3'].str.contains('ERR') ].shape[0]
    good = obs.shape[0] - bad 
    
    if args.verbose:
        print('See {} good and {} bad (B3 field contains "ERR") conversions'.format( good, bad ) )
        
    if args.outfile == '-':
        sys.stdout.write( obs.to_csv(index=None) )
        sys.exit(0)
    obs.to_csv( args.outfile, index=None )
    
import os
import requests

from . import test_helpers

ASTERISM_URL = 'https://raw.githubusercontent.com/AsterismAI/final2000_time_constants/refs/heads/main/daily/reduced_time_constants.dat' 

# -----------------------------------------------------------------------------------------------------
def update_time_constants( text ):
    PATH = test_helpers.get_test_time_constants()
    with open( PATH, 'w' ) as F: F.write( text )
    print('Saw {} lines in data'.format( len( text.split('\n') ) ) )
    print('Wrote to : {}'.format( PATH ) )
    
# -----------------------------------------------------------------------------------------------------
def update_from_github( ):
    data = requests.get( ASTERISM_URL, verify=False )
    update_time_constants( data.text )


# =====================================================================================================
# main
# =====================================================================================================
if __name__ == '__main__':
    import sys
    import os 
    import argparse

    parser = argparse.ArgumentParser(
        prog='utils',
        description='''utilities to manage the install''',
        epilog = ''
    )

    parser.add_argument('--updategithub',
            '-ugh',
            required=False,
            action='store_true',
            help='update the current reduced time constants from GitHub')

    parser.add_argument('--updatefile',
            '-ufl',
            required=False,
            help='update the current reduced time constants from a file',
            default=None)

    
    args = parser.parse_args()
    
    if args.updategithub:
        update_from_github()
        
    if args.updatefile:
        if os.path.isfile( args.updatefile):
            with open( args.updatefile, 'r' ) as F: data = F.read()
            update_time_constants( data )

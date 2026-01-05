import os
import requests

# ASTERISM_URL = 'https://raw.githubusercontent.com/AsterismAI/final2000_time_constants/refs/heads/main/daily/reduced_time_constants.dat' 
ASTERISM_URL = 'https://asterism.ai/timeconstants/reduced_time_constants.dat' 

# -----------------------------------------------------------------------------------------------------
def get_test_dir():
    cwd   = os.path.abspath( os.path.dirname(__file__ ) )
    return os.path.abspath( os.path.join( cwd, '..', '..','tests' ) )

# -----------------------------------------------------------------------------------------------------
def get_data_dir():
    cwd   = os.path.abspath( os.path.dirname(__file__ ) )
    return os.path.abspath( os.path.join( cwd, '..', '..','data' ) )

# -----------------------------------------------------------------------------------------------------
def get_lib_dir():
    return os.path.abspath( os.path.dirname(__file__ ) )

# -----------------------------------------------------------------------------------------------------
def get_test_time_constants():
    ''' 
    if writeable, we can store the time constants in with the library
    '''
    return os.path.join( get_lib_dir(), 'reduced_time_constants.dat' )

# -----------------------------------------------------------------------------------------------------
def update_time_constants( text ):
    PATH = get_test_time_constants()
    with open( PATH, 'w' ) as F: F.write( text )
    print('Saw {} lines in data'.format( len( text.split('\n') ) ) )
    print('Wrote to : {}'.format( PATH ) )

# -----------------------------------------------------------------------------------------------------
def pull_reduced_time_constants_github( ):
    data = requests.get( ASTERISM_URL, verify=False )
    assert data.status_code == 200
    return data.text

# -----------------------------------------------------------------------------------------------------
def update_from_github( ):
    update_time_constants( pull_reduced_time_constants_github() )

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

    parser.add_argument('--pullgithub',
            '-pgh',
            required=False,
            action ='store_true',
            help='pull the reduced time constants file from GitHub and echo the text',
            default=None)
    
    args = parser.parse_args()
    
    if args.updategithub:
        update_from_github()
        
    if args.updatefile:
        if os.path.isfile( args.updatefile):
            with open( args.updatefile, 'r' ) as F: data = F.read()
            update_time_constants( data )

    if args.pullgithub:
        print( pull_reduced_time_constants_github() )
import os
import sys
from datetime import timedelta, timezone

def get_test_dir():
    cwd   = os.path.abspath( os.path.dirname(__file__ ) )
    return os.path.abspath( os.path.join( cwd, '..', '..','tests' ) )

def get_test_time_constants():
    return os.path.join( get_test_dir(), 'reduced_time_constants.dat' )

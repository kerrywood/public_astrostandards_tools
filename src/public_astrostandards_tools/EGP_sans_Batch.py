if __name__ == '__main__':
    import public_astrostandards as PA
    PA.init_all()
    PA.get_versions()

    import time
    from datetime import datetime, timedelta
    import ctypes
    import numpy as np
    import scipy.optimize


# example TLE 
# this is a type-4 faked by a modified from a space-track TLE
L1 = '1 25544U 98067A   24365.67842578  .00026430  00000-0  46140-3 4  9990'
L2 = '2 25544  51.6404  61.8250 0005853  25.4579 117.0387 15.50482079489028'

# what dates should we convert to?
WINDOW_START = datetime( year=2025, month=1, day=10 )
WINDOW_END   = datetime( year=2025, month=1, day=11 )


# clear all sats
PA.TleDll.TleRemoveAllSats()
# load the TLE
tleid = PA.TleDll.TleAddSatFrLines( PA.Cstr(L1,512), PA.Cstr(L2,512) )

# ----------------------------------------------------------------------------------------------
# now, pull the data back out of AstroStandards and into an "introspector" helper class
# per the code, these are the XA_TLE lookups... init a helper to parse those data
XA_TLE = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
# it also passes back some string values.. we'll ignore those
XS_TLE = PA.Cstr('',512)
# now pull..
PA.TleDll.TleDataToArray( tleid, XA_TLE.data, XS_TLE )  # <--- note that you pass the "data" holder in
origvals = XA_TLE.toDict()
# print the data
for k,v in origvals.items():
    print('{:20} {}'.format( k,v ).replace(' ','.') )


# ## Fitting
# 
# So, how will we do this?
# 
# - copy over the data from the original TLE (as our seed orbit), but change the type to 0
# - specify which fields we will optimize over
# - let Python twiddle the bits until we have an answer that best matches our output ephemeris

# In[ ]:


# truth ephemeris; pick one minute spacing because reasons
# convert dates to the astrostandards formats using our helpers
WINDOW_START_DS50 = PA.helpers.datetime_to_ds50( WINDOW_START, PA.TimeFuncDll )
WINDOW_END_DS50   = PA.helpers.datetime_to_ds50( WINDOW_END, PA.TimeFuncDll )
DS50_DATES        = np.arange( WINDOW_START_DS50, WINDOW_END_DS50, 1/1440 )  # astrostandards does days since 1950.. this is one minute steps


# In[ ]:


# helper routine to turn data arrays into TLE lines (using astrostandards)
# data arrays are organized by the helper auto-parser
def arrayToTle( HELPER : PA.helpers.astrostd_named_fields ):
    PA.TleDll.TleRemoveAllSats()
    tleid = PA.TleDll.TleAddSatFrArray( HELPER.data, XS_TLE )
    assert tleid > 0
    outL1, outL2 = PA.Cstr('',512), PA.Cstr('',512)
    assert PA.TleDll.TleGetLines( tleid, outL1, outL2 ) == 0
    return outL1.value.decode('utf-8'), outL2.value.decode('utf-8')


# In[ ]:


# -------------------------------------------------------------------------------------
# propagate routine...
# start with our C-handles to the variables
mse = ctypes.c_double()
pos = (ctypes.c_double * 3)()
vel = (ctypes.c_double * 3)()
llh = (ctypes.c_double * 3)()
def propTle( tleid, ds50 : list ):
    ''' propagate initialized tle to a bunch of dates, return matrix'''
    def propDS50( tleid, date ):
        PA.Sgp4PropDll.Sgp4PropDs50UTC( tleid, date, mse, pos, vel, llh )
        return np.hstack( (date, float(mse.value), list(pos), list(vel)) )   # < -- note the use of list / float to get copies
    return np.vstack( [ propDS50(tleid,D) for D in ds50 ] )

# init the Sgp4 propagator on the tle
assert 0 == PA.Sgp4PropDll.Sgp4InitSat( tleid )
# now propagate
truth  = propTle( tleid, DS50_DATES )
# print(truth[0:2,:])


# # Iteration object : a copy of the seed that we will optimize
# - pick the type (0)
# - also pick the epoch

# In[ ]:


# get another data holder and populate it with the same data as the seed
TEST_TLE = PA.helpers.astrostd_named_fields( PA.TleDll, prefix='XA_TLE_') 
# use a conversion routine to convert to the array *without* loading
PA.TleDll.TleLinesToArray( PA.Cstr(L1,512), PA.Cstr(L2,512), TEST_TLE.data, XS_TLE )
# switch the type from 4 to zero... 
TEST_TLE['XA_TLE_EPHTYPE'] = 0.
# set the epoch to the start of our interval
TEST_TLE['XA_TLE_EPOCH'] = (WINDOW_START_DS50 + WINDOW_END_DS50) / 2


# In[ ]:


for k,v in TEST_TLE.toDict().items():
    print('{:20} {}'.format( k,v ).replace(' ','.') )


# ## Set the seed : set epoch and use the propagated SGP4-XP state closest to epoch as seed
# 
# - reset the epoch on the optimize object
# - find the P,V nearest the epoch
# - find the osculating elements there and set the TLE to those values

# In[ ]:


# reset the seed state to the new epoch, first find the closest date to the epoch
new_epoch = TEST_TLE['XA_TLE_EPOCH']
idx       = np.argmin( np.abs( truth[:,0] - new_epoch ) )
new_sv    = truth[ idx ]
P,V       = truth[idx,2:5], truth[idx,5:]

# we'll use the conversion in the astrostandards
XA_KEP    = PA.helpers.astrostd_named_fields( PA.AstroFuncDll,  prefix='XA_KEP_' )
PA.AstroFuncDll.PosVelToKep( (ctypes.c_double*3)(*P), (ctypes.c_double*3)(*V), XA_KEP.data )

# now set the values
TEST_TLE['XA_TLE_INCLI']  = XA_KEP['XA_KEP_INCLI']
TEST_TLE['XA_TLE_NODE']   = XA_KEP['XA_KEP_NODE']
TEST_TLE['XA_TLE_ECCEN']  = XA_KEP['XA_KEP_E']
TEST_TLE['XA_TLE_MNANOM'] = XA_KEP['XA_KEP_MA']
TEST_TLE['XA_TLE_OMEGA']  = XA_KEP['XA_KEP_OMEGA']
TEST_TLE['XA_TLE_OMEGA']  = XA_KEP['XA_KEP_OMEGA']
TEST_TLE['XA_TLE_MNMOTN'] = PA.AstroFuncDll.AToN( XA_KEP['XA_KEP_A'] )

# save those as our seedvals
seedvals = TEST_TLE.toDict()


# In[ ]:


# what fields will we optimize over?  This doubles as a field accessor list for the optimizer..
FIELDS = [
'XA_TLE_BTERM',
'XA_TLE_NDOT',
'XA_TLE_SP_BTERM',
'XA_TLE_INCLI',
'XA_TLE_NODE',
'XA_TLE_ECCEN',
'XA_TLE_OMEGA',
'XA_TLE_MNANOM',
'XA_TLE_MNMOTN',
]

def optFunction( X, FIELDS, seedvals, return_scalar=True ):
    # --------------------- set our seed values (all vals)
    # do this so that we are always setting off the main TLE
    for k,v in seedvals.items(): TEST_TLE[ k ] = v
    # --------------------- now use modified values
    # just in case there is a bad value in the above
    for k,v in zip(FIELDS,X) :   TEST_TLE[ k ] = v
    # --------------------- clear state
    PA.TleDll.TleRemoveAllSats()
    PA.Sgp4PropDll.Sgp4RemoveAllSats()
    # --------------------- init our test TLE from the modified data
    tleid = PA.TleDll.TleAddSatFrArray( TEST_TLE.data, XS_TLE )
    if tleid <= 0: return np.inf
    if PA.Sgp4PropDll.Sgp4InitSat( tleid ) != 0: return np.inf
    # --------------------- generate our test ephemeris
    test_tle = propTle( tleid, DS50_DATES )
    # use numpy to return the distance between our hypothesis and truth
    resids = test_tle[:,2:] - truth[:,2:]
    rms    = np.sqrt( np.sum( np.linalg.norm( resids, axis=1 ) ) / resids.shape[0] )
    print( rms , end='\r')
    if return_scalar:
        return rms
        # return np.sum( np.linalg.norm( resids[:,:3], axis=1 ) ) 
    else:
        np.linalg.norm( resids[:,:3], axis=1 ) 


# In[ ]:


st= time.time()

# get an initial simplex that is rich in entropy
# seed  = np.array( [ seedvals[k] for k in FIELDS ] )
# N     = len(seed)
# smplx = np.random.uniform(0,.5,size=(N+1,N))
# smplx += seed 

# -----------------------------  nelder mead -----------------------------
# if your seed is not near the final, nelder works great (at the expense of time)
ans = scipy.optimize.minimize( optFunction, 
                               [ seedvals[k] for k in FIELDS ],
                               args    = (FIELDS,seedvals,True),
                               method  = 'Nelder-Mead' )
                               # options = {'initial_simplex' : smplx } )

# -----------------------------  least_sq  -----------------------------
# if your seed IS near the final, follow the gradient
# ans = scipy.optimize.least_squares( optFunction, 
#                                [ seedvals[k] for k in FIELDS ],
#                                args = (FIELDS,seedvals, False) )

print('Optimization / conversion took {:8.3f} seconds'.format( time.time() - st ) )


# In[ ]:


print('Your original TLE was')
print('\n'.join( [L1,L2] ) ) 

# what was our seed value
for k,v in seedvals.items(): TEST_TLE[ k ] = v
print('Your seed TLE was')
print('\n'.join(arrayToTle( TEST_TLE )))

# now update with perturbed values
for k,v in zip(FIELDS,ans.x) :   TEST_TLE[ k ] = v
print('Your updated TLE is')
print('\n'.join(arrayToTle( TEST_TLE )))


# In[ ]:


for i,F in enumerate(FIELDS):
    print('{:15} {:12.8f} {:12.8f}'.format(F, seedvals[F], ans.x[i] ) )


# In[ ]:





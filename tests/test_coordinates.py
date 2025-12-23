import public_astrostandards as PA
import public_astrostandards_tools as PAT

# init the harness to the astrostandards (and setup 'aslog.txt')
PA.init_all()

# run the simple test
try: 
    PAT.coordinates.test()
except Exception as e:
    print(e)
    print(PA.get_last_errmsg())
import public_astrostandards as PA
import public_astrostandards_tools as PAT

import subprocess

# see what versions of the harness libs the OS is seeing
PA.get_versions()

# test the command to pull from github and output to console
# subprocess.run( ['python', '-m','public_astrostandards_tools.utils', '--pullgithub'] )

# test the GitHub auto-update (stored in lib dir / must be writeable)
# subprocess.run( ['python', '-m','public_astrostandards_tools.utils', '--updategithub'] )r

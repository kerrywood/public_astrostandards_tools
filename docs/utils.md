## utils

### Overview 

A collection of routines that is currently mostly focused on getting the time constants file.


### What is it useful for?

- Updating the time constants file for the astrostandards.

### FAQ

- Why does it default to looking for AsterismAI?
    * I host a service that pulls and auto-converts the timing data.  Convenience.


### Invocation
```
python.exe -m public_astrostandards_tools.utils

usage: utils [-h] [--updategithub] [--updatefile UPDATEFILE] [--pullgithub]

utilities to manage the install

options:
  -h, --help            show this help message and exit
  --updategithub, -ugh  update the current reduced time constants from GitHub
  --updatefile UPDATEFILE, -ufl UPDATEFILE
                        update the current reduced time constants from a file
  --pullgithub, -pgh    pull the reduced time constants file from GitHub and echo the text
```
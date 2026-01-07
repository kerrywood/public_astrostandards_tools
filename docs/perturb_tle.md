## perturb_tle

### Overview

`perturb_tle` does a localized approximation of a maneuver on a TLE (at the TLE epoch).  The user provides a TLE and a list of RIC maneuvers (in km/s).  

### How does it work?

The tool does a localized osculating-element approximation of the maneuver.  It:
- propagates the TLE to epoch and gets the resulting position and velocity
- applies the RIC vectors to the osculating state
- approximates the differences in orbital elements due to those perturbations
- generates an ensemble of TLE for those changes

### What is it useful for?
- generating some hypothesis TLE for a potential maneuver

### FAQ
- What if I want to change the TLE epoch as well so that I can maneuver then?
    * _Use the `egp` tool to change the epoch, then use this tool._

### Invocation 

```
python public_astrostandards_tools.perturb_tle

usage: TLE perturb-er [-h] --line1 LINE1 --line2 LINE2 [--RIC RIC] [--satno SATNO]

given a TLE, do an osculating approximation around epoch for your chosen perturbation

options:
  -h, --help            show this help message and exit
  --line1 LINE1, -l1 LINE1
                        line 1 of a TLE
  --line2 LINE2, -l2 LINE2
                        line2 of TLE
  --RIC RIC, -r RIC     RIC perturbations (list of lists; km/s ; JSON parse-able) : [[1,0,0],[.1,.1,.1]]
  --satno SATNO, -N SATNO
                        new TLE satno
```
## UDL EO fitter

### Overview 
 
 This is a convenience tool that will update a TLE with UDL EO obs.  No formatting of the UDL EO obs is required.  Just download them and run it.

### How does it work?

Straight up OD implementation: propagate hypothesis TLE and tune it until it fits the EO obs in the UDL file.  It does **not** weight obs (yet).  

### What is it useful for?
- updating a TLE with new EO obs from UDL

### FAQ
- What happens if you provide bad obs?
    * It'll do the best it can.  Enjoy your garbage.
- Will it always converge?
    * No.
- What does it set the epoch to?
    * The time of the final observation.
- Does it output covariance?
    * Not yet.  It uses Nelder-Mead for optimization so we could output the final simplex.  However, this isn't truly the covariance of the solution.  There is a known algorithm that needs to be implemented to get that (TODO).

### Invocation
```
python.exe -m public_astrostandards_tools.udl_eo_fitter

usage: UDL EO obs fitter [-h] --line1 LINE1 --line2 LINE2 --infile INFILE --outfile OUTFILE [--type TYPE] [--verbose] [--satno SATNO]

take a set of obs downloaded from UDL and an initial TLE, and fit it

options:
  -h, --help            show this help message and exit
  --line1 LINE1, -l1 LINE1
                        line 1 of a TLE
  --line2 LINE2, -l2 LINE2
                        line2 of TLE
  --infile INFILE, -F INFILE
                        load obs from this file
  --outfile OUTFILE, -O OUTFILE
                        store output from loaded jobs
  --type TYPE, -T TYPE  TLE type to fit (0,2,4)
  --verbose, -v         print debugging info
  --satno SATNO, -N SATNO
                        new TLE satno
```
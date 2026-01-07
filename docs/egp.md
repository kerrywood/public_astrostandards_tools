## EGP

### Overview 
 
 Given a TLE and a date range, move the epoch, or convert the type, or both.  The new epoch will always be in the middle of the given range.  Useful for type-4 to type-0 conversions or re-epoching.  As the user, you decide how big of a window you should fit, and what spacing.  

### How does it work?

This tool fuses a couple of underlying tools,  it:

- takes a TLE and propagates it to the user-defined window,
- fits a TLE (with a new epoch and potentially new type) to those ephemerides

### What is it useful for?
- changing a TLE type 
- getting a short-term approximation of type-0 TLE for a type-4 (if you have tooling that cannot accept a type-4, but need the additional accuracy for a future period)
- re-epoching a TLE

### FAQ
- What is `spacing?`
    * Spacing determines the underlying ephemeris spacing.  More points may be more accurate at the cost of compute time. 
- Does it always work?
    * No.  And no warranty is provided.  Good

### Invocation
```
python.exe -m public_astrostandards_tools.egp

usage: egp converter [-h] [--line1 LINE1] [--line2 LINE2] [--startdate STARTDATE] [--enddate ENDDATE] [--spacing SPACING] [--type TYPE] [--satno SATNO] [--test] [--infile INFILE] [--outfile OUTFILE]

take a type-4 TLE (or really any type) and convert it to a fit region

options:
  -h, --help            show this help message and exit
  --line1 LINE1, -l1 LINE1
                        line 1 of a TLE
  --line2 LINE2, -l2 LINE2
                        line2 of TLE
  --startdate STARTDATE, -sd STARTDATE
                        first date (in a format that pd.date_range can understand)
  --enddate ENDDATE, -ed ENDDATE
                        last date (in a format that pd.date_range can understand)
  --spacing SPACING, -sp SPACING
                        ephemeris spacing (in format pd.date_range can understand)
  --type TYPE, -ty TYPE
                        TLE type (0,2,4)
  --satno SATNO, -sno SATNO
                        new TLE satno
  --test, -ts           run a test (all other arguments ignored)
  --infile INFILE, -F INFILE
                        load jobs from a file (other arguments ignored)
  --outfile OUTFILE, -O OUTFILE
                        store output from loaded jobs
```
## EGP

### Overview 
 
Given some UDL EO observations, convert them to B3 type 9.

### How does it work?

This makes use of the ObsDll interface and the semi-standard UDL format.  

- EFG is converted from the sensor lat/lon/alt data and used as the sensor position.  
- UDL is assumed to have RA/DEC in J2K
- `satNo` is interpreted from `satNo` or `origObjectId`
- B3 sensor ID is hashed from the `idSensor` string

### What is it useful for?
Converting a bunch of UDL obs to B3.

### FAQ
- Can I input a list of sensor ID's to map to integers?
    * Not yet.  This is in the ToDo's.
- What if the `satNo` conversion fails because of bad data?
    * It will set the `satNo` in the B3 to 99999
- Can I send to stdout?
    * Yes, set `--outfile -`
- Can I send input from stdin?
    * No, not yet.  Must be a file
- What if the B3 conversion fails (due to bad data or dupes)?
    * You'll see `ERR` as the output for the B3.  Those are returned to you to handle.
- How does it treat the site ID and original tag in the B3?
    * For now, these are just the same as `satNo`.  A to-do is to allow users to specify this in case you want to store data there.
- How do I add weird data to the B3 when I'm done (e.g. vismag)?
    * Honestly, just modify the B3 when you're done.  Treat this as the first step of a conversion.


### Invocation
```
python.exe -m public_astrostandards_tools.udleo_to_b3 

usage: UDL EO obs to B3 [-h] --infile INFILE --outfile OUTFILE [--verbose]

take a set of obs from UDL (directly) and add a B3 column

options:
  -h, --help            show this help message and exit
  --infile INFILE, -F INFILE
                        load obs from this file
  --outfile OUTFILE, -O OUTFILE
                        store output from loaded jobs
  --verbose, -v         print debugging info
```
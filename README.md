# public_astrostandards_tools

Kerry Wood (kerry.wood@asterism.ai)

Updated January 2026

## Overview

`public_astrostandards_tools` is a set of tools that extends the [public_astrostandards](https://github.com/AsterismAI/public_astrostandards) harnesses for the United States Space Force (USSF) Standard Astrodynamic Algorithms Library (SAAL).  
(The SAAL are also known as the "astro-standards." )

Astro-standards has government reference design algorithms for coordinate transforms, data formats (e.g. TLE, B3), and common astrodynamic routines.  It is also the only available implementation for the SGP4-XP propagator.

## Philosophy

**Organization** : The idea is to keep the harnesses (`public_astrostandards`) simple and separate.  
Those harnesses are auto-generated and form a thin layer over the astrostandard shared libraries.  All additional tooling is in here. 

**Tooling** : the SAAL comes in two distributions: (1) the publicly available shared libraries on [space-track.org](https://space-track.org), and (2) a permission-only, ITAR-encumbered set of algorithms.  Set (1) is a subset of (2).  `public_astrostandards_tools` will provide open-source implementations and reference designs for useful tools.  Some of which may or may not be in set (2).

**Workflow** : the astrostandards are FORTRAN routines delivered as compiled shared libraries.  The `public_astrostandards` harnesses are generated via `ctypesgen` from the provided C-headers.  They're a more convenient way to access individual calls, but are still wonky.  `public_astrostandard_tools` are _Python-only code_ designed to turn the harnesses into a _data-science-like_ toolset.  Pandas is used extensively.  Get used to passing around dataframes!

**Portability** : there are other harness codes.  This one is text-only and wraps the DLL's manually.  Since this is just python text, you can move this code to systems that might otherwise be difficult to get to.

## Requirements / Pre-reqs

You'll need to setup the [public_astrostandards](https://github.com/AsterismAI/public_astrostandards) before using this code.  That includes setting up the shared libraries downloaded from [space-track.org](https://space-track.org).  If you can import that library and run `get_versions()`, you're all set.

**>>> YOU <<<**  are responsible for getting permission, acquiring, and setting up the shared libraries.

# Examples

There are different ways to use this package.  Command line tools for key features are included.  Jupyter notebooks are included as instructional examples.  Test code is a good reference as well.  All are linked below.

## Example notebooks

Data science examples are available in the [notebooks](./docs/notebooks).  

## Command line tools

Some utilities have been built as command line tools for convenience.  Always looking for suggestions on new additions.

[EGP](./docs/egp.md) : re-epoch a TLE or convert the type (or both)

[perturb tle](./docs/perturb_tle.md) : perturb a TLE at epoch with RIC changes

[UDL EO fitter](/docs/udl_eo_fitter.md) : using EO observations from the UDL, perform OD on a TLE

[UDL EO B3 converter](/docs/udlob_to_b3.md) : take raw UDL obs and convert them to B3

## Test code

There are lots of test routines that help me debug located in [./tests](./tests)

They are good references for those learning the code






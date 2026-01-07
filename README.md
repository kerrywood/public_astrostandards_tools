# public_astrostandards_tools

Kerry Wood (kerry.wood@asterism.ai)

Updated January 2026

## Overview

`public_astrostandards_tools` is a set of tools that extends the [public_astrostandards](https://github.com/AsterismAI/public_astrostandards) harnesses for the United States Space Force (USSF) Standard Astrodynamic Algorithms Library (SAAL).  The SAAL are also known as the "astro-standards." Astro-standards has government reference design algorithms for coordinate transforms, data formats (e.g. TLE, B3), and common astrodynamic routines.  Additionally, the library is the only source for the SGP4-XP propagator.

## Philosophy

**Organization** : The idea is to keep the harnesses (`public_astrostandards`) simple and separate, and add all additional tooling here. 

**Tooling** : the SAAL has two distributions: (1) the publicly available shared libraries on [space-track.org](https://space-track.org), and (2) a permission-only, ITAR-encumbered set of algorithms.  Set (1) is a subset of (2).  `public_astrostandards_tools` will provide open-source implementations and reference designs for useful tools that may or may not be in set (2).

**Workflow** : the astrostandards are FORTRAN routines delivered as compiled shared libraries.  The `public_astrostandards` harnesses are generated via `ctypesgen` from the provided C-headers.  `public_astrostandard_tools` are Python-only code designed to turn the code into a _data-science-like_ toolset.  Pandas is used extensively.  Get used to passing around dataframes!

**Portability** : there are other harness codes.  This one is text-only and wraps the DLL's manually.  You can move this code to systems that might otherwise be difficult to get to.

## Requirements / Pre-reqs

You'll need to setup the [public_astrostandards](https://github.com/AsterismAI/public_astrostandards) before using this library.  That includes setting up the shared libraries downloaded from [space-track.org](https://space-track.org).  If you can import that library and run `get_versions()`, you're all set.

**YOU ARE RESPONSIBLE FOR ACQUIRING AND SETTING UP THE ASTROSTANDARD SHARED LIBRARIES.**

## Example notebooks

Data science examples are available in the [notebooks](./docs/notebooks).  

## Command line tools

Some utilities have been built as command line tools for convenience.  Always looking for suggestions on new additions.

[EGP](./docs/test.md) : re-epoch a TLE or convert the type (or both)

[perturb tle](./docs/perturb_tle.md) : perturb a TLE at epoch with RIC changes

[UDL EO fitter](/docs/udl_eo_fitter.md) : using EO observations from the UDL, perform OD on a TLE








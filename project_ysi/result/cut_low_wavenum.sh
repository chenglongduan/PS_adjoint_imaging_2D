#!/bin/bash

# Original wavenumber spectrum
< test1.rsf sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" > spectra_test1.vpl

# Spectrum after cutting low wavenumber (equivalent to boost high wavenumber)
< multishots.rsf sfbandpass flo=0.003 | sfspectra all=y | sfscale axis=1 | sfgraph label2="Amplitude" unit2="" > spectra_lowcut_60shot.vpl

# Full image processing (water depth = 255 m) (<250m W=0; >260m W=1; interpolate W between 250m and 260m)
< multishots.rsf sfbandpass flo=0.003 | sfmutter half=n slope0=0 slopep=0 t0=250 tp=10 | sfgrey pclip=99 labelsz=10 labelfat=4 screenratio=0.6 wanttitle=n > multishots.vpl

# Save processed RTM
< multishots.rsf sfbandpass flo=0.003 | sfmutter half=n slope0=0 slopep=0 t0=250 tp=10 > multishots_bp_mute.rsf


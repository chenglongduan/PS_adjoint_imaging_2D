# PS_adjoint_imaging_2D

This code package can perform full wavefield seismic reflection imaging using adjoint P and S wave equations, which is parallelized using Message Passing Interface (MPI) techniques under C programming. Note that this code was originally developed to run on Intel HPC slurm system. Be aware of errors running on any other platforms.

This method [1] has been applied to volcanic structure imaging submitted to a journal paper [2].

Please cite the GJI methodology paper [1] and the submitted volcano application paper [2], if you use this code for your publications:

[1] Chenglong Duan, David Lumley, Hejun Zhu. "Estimation of micro-earthquake source locations based on full adjoint P and S wavefield imaging." Geophysical Journal International 226 (2021): 2116–2144. https://doi.org/10.1093/gji/ggab203

[2] Chenglong Duan, Wenkai Song, Brandon Schmandt, Jamie Farrell, David Lumley, Tobias Fischer, Lindsey Worthington, Fan-Chi Lin. "A sharply defined volatile-rich cap of magma reservoir beneath Yellowstone caldera." Submitted (2024).

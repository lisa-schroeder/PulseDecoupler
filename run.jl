#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      SETUP OF BASIS
#------------------------------------------------------------------------------------------------------------------------------------------------

# Packages 
using LinearAlgebra
using Plots
using FFTW
using Printf
using Rotations
using StaticArrays
using Dierckx
using LsqFit
using Pkg

# %Complex definition
const global i = sqrt(Complex(-1))

# Definition of individual spin operators
include("basis.jl")
include("converter.jl")
include("fid_2spins_heteronoclear.jl")
include("fid_2spins_homonuclear.jl")
include("fid_pulsesequences.jl")
include("generate_composite.jl")
include("initialization.jl")
include("plot.jl")
include("read_wave.jl")
include("sidebands.jl")
include("supercycles.jl")
include("waugh.jl")
include("write_wave.jl")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      SETUP OF SIMULATION PARAMETERS AND SIMULATIONS
#------------------------------------------------------------------------------------------------------------------------------------------------

nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params.txt")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOT RESULTS OF SIMULATION
#------------------------------------------------------------------------------------------------------------------------------------------------


# exemplary plots in plot.jl


#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      SAVE PULSE SHAPE IN XYT OR XYZT FORMAT
#------------------------------------------------------------------------------------------------------------------------------------------------


# info = String("# test pulse\n# rfmax = 4 kHz")
# save_file = "C:/Users/user/test.xyt"
# xyt_save(save_file, info, ux, uy, ut)
# xyzt_save(save_file, info, ux, uy, uz, ut)



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      FIT SIDEBANDS AND MAXIMUM TO J^2
#------------------------------------------------------------------------------------------------------------------------------------------------


# call functions from "sidebands.jl"
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
using QuadGK
using CSV
using DataFrames
using StatsBase
using Measures
using LaTeXStrings
using Test




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

nB1, noffs, fB1, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_cw.txt")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOT RESULTS OF SIMULATION
#------------------------------------------------------------------------------------------------------------------------------------------------


# exemplary plots in plot.jl

# Linewidth 6 Hz
# display(plot(offs./1000, projection[ceil(Int32,nB1/2),:,:]./2.505, legend=false, xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", linewidth=2, ylims=(0,100), dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./2.505, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))

# Linewidth 1.6 Hz
display(plot(offs./1000, projection[ceil(Int32,nB1/2),:,:]./10, legend=false, xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", linewidth=2, ylims=(0,100), dpi=500))
# savefig("filename.svg")
display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))

# display(contour(offs./1000, fB1, projection./10, xlabel=L"\nu_I"*" (kHz)", clims=(0, 100), ylabel=L"B_{1,rel}", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]', ylims=(-20,20), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))

# display(contour(offs./1000, fB1, Mz, label="Mz", xlabel=L"\nu_I"*" (kHz)", clims=(-1,1), ylabel=L"B_{1,rel}", color=cgrad([:navyblue, :white, :darkred], [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1], categorical = true), fill=true, linewidth=0, dpi=500))
# display(plot(offs./1000, abs.(lambda), legend=false, ylims=(0,0.5), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# display(plot(offs./1000, [Mx[1,:] My[1,:] Mz[1,:]], label=["Mx" "My" "Mz"], linewidth=2, xlabel="offset (kHz)", ylabel="Magnetization"))


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

# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), allspectr[1,ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# filename = "C:/Users/Lisa/Documents/KIT/WS202324/tests/caWURST2-M4P5"
# write_maxintensities_to_csv(maxintens, maxintensside, noffs, filename)
# Js, sidebands = read_in_biggest_sidebands(filename*"_sidebands.csv")
# plot_fit_sidebands(Js, sidebands, (40,260), (0,5), [0.002], noffs)
# Js, maxvals = read_in_biggest_sidebands(filename*"_maxintens.csv")
# plot_fit_maxintens(Js, maxvals, (40,260), (80,100), [-0.002], noffs)


# pulsename = "WALTZ16"
# pulsename = "caWURST2-M4P5"
# pulsename = "BROCODE_lvl4_linewidth6Hz"

# Js, sidebands = read_in_biggest_sidebands("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/all_same_rfav/Sidebands/all_offsets/"*pulsename*"_Sidebands.csv")
# plot_fit_sidebands(Js, sidebands, (40,260), (0,5), [0.002], noffs)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/tests/"*pulsename*"_sidebands_all.png")
# Js, maxvals = read_in_biggest_sidebands("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/all_same_rfav/Sidebands/all_offsets/"*pulsename*"_maxintens.csv")
# plot_fit_maxintens(Js, maxvals, (40,260), (88,100), [-0.002], noffs)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/tests/"*pulsename*"_maxvals_all.png")
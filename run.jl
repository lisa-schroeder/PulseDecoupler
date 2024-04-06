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

# using Pkg
# pkg"add FFMPEG"
# pkg"build FFMPEG"
# pkg"build Plots"

# %Complex definition
const global i = sqrt(Complex(-1))

# Definition of individual spin operators
include("averageH.jl")
include("basis.jl")
include("converter.jl")
include("fid_2spins_heteronoclear.jl")
include("fid_2spins_homonuclear.jl")
include("fid_pulsesequences.jl")
include("generate_composite.jl")
include("gradients.jl")
include("initialization.jl")
include("plot.jl")
include("read_wave.jl")
include("supercycles.jl")
include("waugh.jl")
include("write_wave.jl")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                SETUP OF SIMULATION PARAMETERS AND SIMULATIONS
#------------------------------------------------------------------------------------------------------------------------------------------------

nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params.txt")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOTTING OF THE RESULTS OF FID SIMULATIONS
#------------------------------------------------------------------------------------------------------------------------------------------------

# savingname="Figures/B1_Inhomogeneity/finished/WALTZ16"
# title="continuous wave 5 kHz"
# display(plot_pulse_bruker(filename, tpulse))
# display(plot_ampl_phase(acqux[2], acquy[2], acqut[2], ""))

# ampl, phase = xy_to_ampl_phase(acqux[2], acquy[2])
# ptime = get_ptgrid(acqut[2])
# plot(ptime.*1000, phase, label="xy-Phase", ylabel="Phase (°)", linetype=:steppre, fill=(0,0.2), xlabel="Time (ms)", linewidth=2, dpi=500)
# invis = ones(size(ptime)[1]).*(-100)
# plot!(ptime, invis, label= "rf-amplitude", color=:darkorange, ylims=(-5,365), linewidth=2, legend=:top, dpi=500)
# plot!(twinx(), ptime.*1000, ampl./1000, legend=false, label="Phase (°)", color=:darkorange, ylims=(0, 14), ylabel="Amplitude (kHz)", linewidth=2, dpi=500)
# display(plot(ptime.*1000, acquz[2]./1000,label="uz", linetype=:steppre, fill=(0, 0.5), legend=false, color=:green, ylabel="offset (kHz)", xlabel="time (ms)"))

# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/Bulus_Vortrag/WURST2-M4P5_decoupling_profile.png")
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_decoupled_peak.png")
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/Theory/Spectrum_Coupled_no_grid.png")

# display(plot_ux_uy(acqux[2][1:1000], acquy[2][1:1000], acqut[2][1:1000], (0, sum(acqut[2][1:1000])*1000), (-4,4), ""))

# savefig(savingname*"_ux_uy.png")
# display(plot_ux_uy(acqux[2], acquy[2], acqut[2], (0.0, sum(acqut[2])*1000/8), (-4.81,4.81), ""))
# display(plot_ux_uy_uz(acqux[2], acquy[2], acquz[2], acqut[2]))
# display(plot_pulse_phase(ux, uy, ut))
# display(plot_pulse_phase_ampl_offs_contour(ux, uy, uz, ut, spectr[ceil(Int32,nB1/2),:,:], npoints, dwell, offs))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32, nB1/2),:,:]'./10,  ylims=(-100, 100), clims=(0,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# display(annotate!(0, -50, "\\Xi = " * string(round(xi; digits = 1))))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :red, :yellow], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9,0.95, 0.98, 0.99], categorical = true), clims=(-10, 100), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-100, 100), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, clims=(0,100),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# display(contourf(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :red, :yellow], [-0.1, -0.05, 0.04, 0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9,0.95, 0.98, 0.99], categorical = true),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=palette([:darkblue, :white, :darkred], 21), clims=(-100,100),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig(savingname*"_decoupling_profile.png")
# display(surface(fftshift(fftfreq(npoints,1/dwell)),offs./1000, spectr[ceil(Int32,nB1/2),:,:]./2.505, xlims=(-100,100), zlims=(-5,100), zlabel="intensity (%)", xlabel=L"\nu_S"*" (Hz)", xguidefontsize=10, ylabel=L"\nu_I"*" (kHz)", yguidefontsize=10, dpi=500, camera=(40,30), c=cgrad([:white, :blue, :royalblue4, :darkblue]), cbar=false))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/WURST40-MLEV8_decProfile_sidebands.png")
# savefig(savingname*"_stacked_spectra.png")

# xaxis, yaxis, zaxis = collect(fftshift(fftfreq(npoints,1/dwell)))[2300:2700], offs, spectr[ceil(Int32,nB1/2),:,2300:2700]'
# plot(repeat(xaxis, 1, length(yaxis)), repeat(yaxis', length(xaxis), 1)./1000, zaxis, color =:darkblue, xlabel=L"\nu_S"*" (Hz)", ylabel=L"\nu_I"*" (kHz)", zlabel="intensity", legend=false, camera=(30,30))

# test_spectr_1 = fft_fid(test_single_fid[1,:,:,:], nB1, noffs, npoints)
# display(plot_1d(spectr[ceil(Int32,nB1/2),ceil(Int32, noffs/2),:]./9.20, title, (-250, 250), (-8, 10), npoints, dwell))
# println("Intens: ", findmax(spectr[1,1,2550:end])[1], " at freq ", collect(fftshift(fftfreq(npoints,1/dwell)))[findmax(abs.(spectr[1,1,2550:end]))[2]+2549])
# display(plot_1d(spectr[1,2,:], title, (-250, 250), (-8, 1000), npoints, dwell))
# println("Intens: ", findmax(spectr[1,2,2530:end])[1], " at freq ", collect(fftshift(fftfreq(npoints,1/dwell)))[findmax(abs.(spectr[1,2,2530:end]))[2]+2529])

# plot(collect(fftshift(fftfreq(npoints,1/dwell))).*0.002,test_spectr_1[ceil(Int32,nB1/2),ceil(Int32, noffs/2),:]./10, label="first FID", xlabel=L"\nu_S \cdot T_p", ylabel="intensity (%)", xlims=(-1.5,1.5), ylims=(-2,50), linewidth=2)
# plot!(collect(fftshift(fftfreq(npoints,1/dwell))).*0.002,test_spectr_2[ceil(Int32,nB1/2),ceil(Int32, noffs/2)+5,:]./10, label="second FID", xlabel=L"\nu_S \cdot T_p", ylabel="intensity (%)", xlims=(0,1.5), ylims=(-2,2), linewidth=2)
# plot!(collect(fftshift(fftfreq(npoints,1/dwell))).*0.002,spectr[ceil(Int32,nB1/2),ceil(Int32, noffs/2)+5,:]./10, label="averaged FID", xlabel=L"\nu_S \cdot T_p", ylabel="intensity (%)", xlims=(0,1.5), ylims=(-2,2), linewidth=2)
# savefig(savingname*"_suppressed_harmonics_offs_13kHz.png")
# sidebands = zeros(15, npoints)
# icount = 1
# for ioffs  = 9:23
#     sidebands[icount,:] = spectr[ceil(Int32,nB1/2),ioffs,:]
#     icount += 1
# end
# plot(collect(fftshift(fftfreq(npoints,1/dwell))).*0.002,sidebands'./10, xlabel=L"\nu_S \cdot T_p", xlims=(0,1.5), color=:darkblue, legend=false, ylabel="intensity (%)", ylims=(-1,2), dpi=500)
# savefig(savingname*"_sidebands_stacked_offs_bigger.png")
# display(contour(offs./1000,collect(fftshift(fftfreq(npoints,1/dwell))).*0.001, spectr[ceil(Int32,nB1/2),:,:]'./10, ylims=(-1.5,1.5), clims=(-1,1), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S \cdot T_p", dpi=500))
# savefig(savingname*"_decoupling_profile_sidebands.png")

# display(plot_1d(allspectr[1,ceil(Int32,nB1/2), 15,:]./10, "", (-2500,2500), (0.04, 100), npoints, dwell))
display(plot(fftshift(fftfreq(npoints,1/dwell)),spectr[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:]./10, ylims=(-5,60), xlims=(-250, 250), grid=false, dpi=500, thickness_scaling=1.2, linewidth=2, legend=false, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)"))
# display(plot(fftshift(fftfreq(npoints,1/dwell)),allspectr[1,1,1,:]./10, ylims=(0,0.5), xlims=(-1000, 1000), grid=false, dpi=500, thickness_scaling=1.2, linewidth=2, legend=false, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)"))
# test_temp = test ./ maximum(test) * maximum(allspectr[1,1,1,:])
# plot!(fftshift(fftfreq(npoints,1/dwell)), (allspectr[1,1,1,:]-test_temp)/10, xlims=(-2500,2500), ylims=(-0.2,0.6), linewidth=2, color=:green, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)", label="difference spectrum")
# plot!(fftshift(fftfreq(npoints,1/dwell)), test_temp./10, ylims=(-0.2,0.6), color=:orange, linewidth=2, label="perfectly decoupled")
# plot(fftshift(fftfreq(npoints,1/dwell)), allspectr[1,1,1,:]./10, linewidth=2, label="decoupled by pulse sequence")
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_perfectly_decoupled_difference.png")
# display(scatter(tim[480:520], real(fid[ceil(Int32,nB1/2),ceil(Int32,noffs/2),480:520]), legend=false))
# display(plot(tim, real(fid[ceil(Int32,nB1/2),ceil(Int32, noffs/2),:]), ylabel=L"I^m_{I1}", xlabel=L"\nu_{I2}"*" (kHz)", label="offset 0 Hz", linewidth=2, dpi=500))
# display(plot!(tim, real(fid[ceil(Int32,nB1/2),ceil(Int32, noffs/2)+12,:]), label="offset 2250 Hz", linewidth=2, dpi=500))
# display(plot!(tim, real(fid[ceil(Int32,nB1/2),ceil(Int32, noffs/2)+4,:]), label="offset 3000 Hz", linewidth=2, dpi=500))
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), clims=(0,1000), spectr[ceil(Int32,nB1/2),:,:]', ylims=(-20,20), xlabel=L"\nu_{I2}"*" (kHz)", ylabel=L"\nu_{I1}"*" (Hz)", dpi=500))
# fid_test = [fid[1,9,dig]/exp(-(dig-1)*0.0002/0.2)*cos(dig/1000)*1000/dig for dig in 1:npoints]
# display(plot(tim, real(fid_test), ylabel=L"I^m_{I1}", grid=false, legend=false, ylim=(-1.1,1.1), xlabel=L"\nu_{I2}"*" (kHz)", label="offset 0 kHz", linewidth=2, dpi=500))
# plot(tim.*1000,[real(fid[1,1,:]) real(waltz_fid[1,1,:])], thickness_scaling=1.2, linewidth=2, label=["WALTZ-16" "perfect decoupling"], xlabel="time (ms)", ylabel=(L"I^m"*" magnetization"), xlims=(0,2.4))
# plot(tim,real(fid[ceil(Int32,nB1/2),ceil(Int32, noffs/2),:]), thickness_scaling=1.2, legend=false, xlabel="time (s)", ylabel=(L"I^m"*" magnetization"))
# display(plot(tim,real(allfids[ceil(Int32,nB1/2),ceil(Int32, noffs/2),1,:]), thickness_scaling=1.2, legend=false, xlabel="time (s)", ylabel=("x magnetization")))
# display(plot_fid(tim, sfid[ceil(Int32,nB1/2),ceil(Int32, noffs/2),11,:]))
# display(plot_bandwidth("MLEV-4", spectr[ceil(Int32,nB1/2),:,:], Int(ceil(npoints/2))+1))      # plot spectrum at frequency 0
# display(plot_projection("", false, savingname, projection[ceil(Int32,nB1/2),:,:]./10, offs))  
# plot(offs./1000, proj_MLEV4[1,:], xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity", label="MLEV-4", linewidth=2)
# plot!(offs./1000, proj_MLEV16[1,:], xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity", label="MLEV-16", linewidth=2)
# plot!(offs./1000, proj_MLEV64[1,:], xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity", label="MLEV-64", linewidth=2)
# savefig(savingname*"_projection.png")
# display(annotate!(0, 400, "\\Xi = " * string(round(xi; digits = 1))))
# savefig(savingname*"_projection.png")                             # plot projection, maximum of each offset
# contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[1,:,:]', ylims=(-100, 100), xlabel="offset (kHz)", ylabel="freq (Hz)", title=(title))
# B1 = get_B1()*rfmax
# contour(B1./1000,fftshift(fftfreq(npoin ts,1/dwell)), spectr[:,ceil(Int32,noffs/2),:]', ylims=(-50, 50), xlabel="B1 max (kHz)", ylabel="freq (Hz)", title=(title))
# display(contour(offs./1000, B1./1000, projection, xlabel=L"\nu_I"*" (kHz)", clims=(0, 1000), ylabel="B1 max (kHz)", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0))
# savefig("Figures/B1_inhomogeneity/Version_2/WURST_MLEV8_B1_inhomogeneity.png")
# contour(offs./1000, B1, projection, xlabel="offs (kHz)", clims=(500, 1000), ylabel="B1 max (kHz)", c=cgrad(palette([:white, :grey16], 21), 10, categorical = true, scale = :(1-exp)), fill=true, linewidth=1, title=title*" - B1-inhomogeneity")
# contour(offs./1000, B1, projection, xlabel="offs (kHz)", clims=(500, 1000), ylabel="B1 max (kHz)", c=cgrad(palette([:white, :grey16], 31), 31, categorical = true, scale=:log), fill=true, linewidth=0, title=title*" - B1-inhomogeneity")
# savefig(savingname*"_grey_no_linewidth_log_starting_zero.png")
# contour(offs./1000, B1, projection, xlabel="offs (kHz)", clims=(0, 1000), ylabel="B1 max (kHz)", color=palette([:white, :grey16], 21), fill=true, linewidth=1, title=title*" - B1-inhomogeneity")
# contour(offs./1000, B1, projection, xlabel="offs (kHz)", clims=(800, 1200), ylabel="B1 max (kHz)", color=palette([:white, :maroon], 21), linewidth=1, title=title*" - B1-inhomogeneity")
# heatmap(offs./1000, B1, projection, xlabel="offs (kHz)", ylabel="B1 max (kHz)", color=palette([:white, :darkred], 21))
# display(plot(offs./1000, [Mx[1,:] My[1,:] Mz[1,:]], label=["Mx" "My" "Mz"], linewidth=2, xlabel="offset (kHz)", ylabel="Magnetization"))
# display(plot(offs./1000, Mz[1,:], label="Mz", xlabel="offset (kHz)", linewidth=2, ylims=(-1,0), ylabel="Magnetization"))
# display(contour(offs./1000, fB1./1000, Mz, label="Mz", xlabel="offset (kHz)", clims=(-1,1), ylabel="rf amplitude (kHz)", color=cgrad([:navyblue, :white, :darkred], [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1], categorical = true), fill=true, linewidth=0, dpi=500))
# display(contour(offs./1000, fB1./1000, Mz, label="Mz", xlabel="offset (kHz)", clims=(-1,1), ylabel="rf amplitude (kHz)", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500))
# savefig("Figures/HomoD/Rsnob_decoupling_profile_clim1000_dchunk025.png")
# savefig("Figures/HomoD/Rsnob_FID_dchunk1.png")
# savefig("Figures/Magnetization/180x-B1_closeup.png")
# Mzcomp180 = Mz[1,:]
# Mzcomp270 = Mz[1,:]
# Mzcomp240 = Mz[1,:]
# Mz180 = Mz[1,:]
# display(plot(offs./1000, [Mzcomp180 Mzcomp270], label=["180x" "90x-270y-90x"], linewidth=2, dpi=500, ylims=(-1,0), xlabel="offset (kHz)", ylabel="z-Magnetization"))
# savefig("Figures/Magnetization/Composite_Comparison.png")
# plot_sidebands_stacked_offs(spectr, noffs, offs, (0, 2500), (-0.5, 2), "", dwell, npoints)
# savefig(savingname*"_stacked_J250.png")
# plot_sidebands_stacked_J(allspectr, (-10, 900), (-10, 50), ceil(Int32, nB1/2), 6, allJs, title, dwell, npoints)
# savefig(savingname*"_stacked_offs0kHz.png")

# @show noffs


#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOTTING FOR MY MASTERS THESIS
#------------------------------------------------------------------------------------------------------------------------------------------------

pulsename = "WURST_exp_comp"

# plot_ampl_phase_2(acqux, acquy, acqut)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_ampl_phase.png")
# display(plot_ux_uy(acqux[2], acquy[2], acqut[2], (0, sum(acqut[2])*1000), (-10.15,10.15), ""))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_ux_uy.png")

# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad(:default, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_2percent_all100kHz.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad(:default, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-100, 100), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_100percent_all100kHz.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, clims=(0,100),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_default_2percent_all100kHz.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, clims=(-2,2),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_default_100percent_all100kHz.png")
# display(plot(offs./1000, projection[ceil(Int32,nB1/2),:,:]./2.505, legend=false, xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", linewidth=2, ylims=(0,100), dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_projection_all100kHz.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./2.505, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_newcolors_2percent_all100kHz.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./2.505, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-100, 100), ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_newcolors_100percent_all100kHz.png")

# xlimit = (-30,30)
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad(:default, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), xlims=xlimit, ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_2percent.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad(:default, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-100, 100), xlims=xlimit, ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_100percent.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, xlims=xlimit, clims=(0,100),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_default_2percent.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, xlims=xlimit, clims=(-2,2),ylims=(-100,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_default_100percent.png")
# display(plot(offs./1000, projection[ceil(Int32,nB1/2),:,:]./10, legend=false, xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", xlims=xlimit, ylims=(0,100), linewidth=2, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_projection.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-2, 2), xlims=xlimit, ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_newcolors_2percent.png")
# display(contour(offs./1000,fftshift(fftfreq(npoints,1/dwell)), spectr[ceil(Int32,nB1/2),:,:]'./10, color=cgrad([:darkblue, :blue, :white, :red, :darkred], [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1], categorical = true), clims=(-100, 100), xlims=xlimit, ylims=(-500,500), xlabel=L"\nu_I"*" (kHz)", ylabel=L"\nu_S"*" (Hz)", dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_dec_contour_defined_newcolors_100percent.png")

# spectrWURST2 = spectr
# spectrWURST6 = spectr
# spectrWURST7 = spectr
# spectrWURST8 = spectr

# spectr13723 = spectr
# proj1723 = projection[1,:]
# display(plot(offs./1000, [spectrBUSS5_v2[1,:,2049]./8.673856 spectrBUSS4[1,:,2049]./9.0210832 spectrBUSS3[1,:,2049]./9.0210832 spectrBUSS2[1,:,2049]./9.0210832], legend=:bottom, color=[:purple4 :steelblue2 :limegreen :darkred], label=[L"\nu_{rf,max}="*" 11478 Hz" L"\nu_{rf,max}="*" 12271 Hz" L"\nu_{rf,max}="*" 13018 Hz" L"\nu_{rf,max}="*" 13723 Hz"], xlabel=L"\nu_I"*" (kHz)", ylabel="intensity at "*L"\nu_S=0"*" Hz (%)", ylims=(0,100), linewidth=1, dpi=500))
# display(plot(offs./1000, [spectrWURST8[1,:,2049] spectrWURST7[1,:,2049] spectrWURST6[1,:,2049] spectrWURST2[1,:,2049]]./8.673856, legend=:topright, color=[:purple4 :steelblue2 :limegreen :darkred], label=[L"\nu_{rf,max}="*" 1724 Hz" L"\nu_{rf,max}="*" 2012 Hz" L"\nu_{rf,max}="*" 2414 Hz" L"\nu_{rf,max}="*" 3020 Hz"], xlabel=L"\nu_I"*" (kHz)", ylabel="intensity at "*L"\nu_S=0"*" Hz (%)", ylims=(0,100), linewidth=1, dpi=500))
# display(plot(offs./1000, [spectr1723[1,:,2501] spectr2011[1,:,2501] spectr2414[1,:,2501] spectr3019[1,:,2501]]./2.505, legend=:topright, color=[:purple4 :steelblue2 :limegreen :darkred], label=[L"\nu_{rf,max}="*" 1724 Hz" L"\nu_{rf,max}="*" 2012 Hz" L"\nu_{rf,max}="*" 2414 Hz" L"\nu_{rf,max}="*" 3020 Hz"], xlabel=L"\nu_I"*" (kHz)", ylabel="intensity at "*L"\nu_S=0"*" Hz (%)", ylims=(-10,100), linewidth=1, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_slice_nuS_0Hz.png")
# display(plot(offs./1000, [proj11477 proj12271 proj13017 proj13723]./2.505, color=[:purple4 :darkturquoise :greenyellow :darkred], xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", ylims=(0,100), linewidth=2, dpi=500))

# Simple Waugh
# display(plot(offs./1000, abs.(lambda), legend=false, ylims=(0,0.1), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Waugh_simple_lamda01_all100kHz.png")
# display(plot(offs./1000, abs.(lambda), legend=false, ylims=(0,1), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Waugh_simple_lamda1_all100kHz.png")

# xlimit = (-30,30)
# display(plot(offs./1000, abs.(lambda), legend=false, xlims=xlimit, ylims=(0,0.1), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Waugh_simple_lambda01.png")
# display(plot(offs./1000, abs.(lambda), legend=false, xlims=xlimit, ylims=(0,1), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Waugh_simple_lambda1.png")


# B1 inhomogeneity
# display(contour(offs./1000, B1./1000, projection./2.505, xlabel=L"\nu_I"*" (kHz)", clims=(0, 100), ylabel=L"B_{1}"*" (kHz)", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500))
# display(contour(offs./1000, fB1, projection./10, xlabel=L"\nu_I"*" (kHz)", clims=(0, 100), ylabel=L"B_{1,rel}", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_B1_inhomogeneity_all100kHz.png")
# xlimit = (-10,10)
# display(contour(offs./1000, fB1, projection./10, xlabel=L"\nu_I"*" (kHz)", xlims=xlimit, clims=(0, 100), ylabel=L"B_{1,rel}", color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_B1_inhomogeneity.png")


# Magnetization
# display(contour(offs./1000, fB1, Mz, label="Mz", xlabel=L"\nu_I"*" (kHz)", clims=(-1,1), ylabel=L"B_{1,rel}", color=cgrad([:navyblue, :white, :darkred], [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1], categorical = true), fill=true, linewidth=0, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Magnetization_all80kHz.png")
# display(contour(offs./1000, fB1, Mz, label="Mz", xlabel=L"\nu_I"*" (kHz)", xlims=xlimit, clims=(-1,1), ylabel=L"B_{1,rel}", color=cgrad([:navyblue, :white, :darkred], [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1], categorical = true), fill=true, linewidth=0, dpi=500))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Magnetization.png")

# sidebands: siehe unten

# # J dependence
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J10Hz.txt")
# allJproj = zeros(25, noffs)
# allJproj[1,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J20Hz.txt")
# allJproj[2,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J30Hz.txt")
# allJproj[3,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J40Hz.txt")
# allJproj[4,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J50Hz.txt")
# allJproj[5,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J60Hz.txt")
# allJproj[6,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J70Hz.txt")
# allJproj[7,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J80Hz.txt")
# allJproj[8,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J90Hz.txt")
# allJproj[9,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J100Hz.txt")
# allJproj[10,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J110Hz.txt")
# allJproj[11,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J120Hz.txt")
# allJproj[12,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J130Hz.txt")
# allJproj[13,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J140Hz.txt")
# allJproj[14,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J150Hz.txt")
# allJproj[15,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J160Hz.txt")
# allJproj[16,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J170Hz.txt")
# allJproj[17,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J180Hz.txt")
# allJproj[18,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J190Hz.txt")
# allJproj[19,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J200Hz.txt")
# allJproj[20,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J210Hz.txt")
# allJproj[21,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J220Hz.txt")
# allJproj[22,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J230Hz.txt")
# allJproj[23,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J240Hz.txt")
# allJproj[24,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]
# nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm = performsimulations("params_J250Hz.txt")
# allJproj[25,:] = [maximum(spectr[1, ioffs, :]) for ioffs = 1:noffs]

# using HDF5

# h5write("Figures/testing/allJproj_BUSS7kHz", "BUSS7kHz246ms4", allJproj)
# data = h5read("/tmp/test2.h5", "mygroup2/A", (2:3:15, 3:5))
# # where the last line reads back just A[2:3:15, 3:5] from the dataset.


# minJproj = [minimum(allJproj[iJ,11:91]) for iJ=1:size(allJproj)[1]]./2.504983791023052
# averageJproj = [sum(allJproj[iJ,11:91])./81 for iJ=1:size(allJproj)[1]]./2.504983791023052

# minJproj1kHz75 = minJproj
# averageJproj1kHz75 = averageJproj

# display(contour(offs./1000, [250,260,270,280,290,300,310,320,330,340,350], allJproj[1:11,:]./2.5, xlims=(-40,40,), color=palette([:darkblue, :white, :darkred], 21), fill=true, linewidth=0, dpi=500, clims=(90,100), xlabel=L"\nu_I"*" (kHz)", ylabel=L"J"*" (Hz)"))
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Jdependence_upto250Hz.png")



#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOTTING OF THE RESULTS OF AVERAGE HAMILTONIAN
#------------------------------------------------------------------------------------------------------------------------------------------------


# savingname = "C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/average_H/hard180"
# title="WURST average H"
# plot_xyz_one_spin(title, false, savingname, ptgrid, hamx, hamy, hamz, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_xyz_one_spin(title, false, savingname, dwelltime, hamdwellx, hamdwelly, hamdwellz, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_xyz_sum_one_spin(title, false, savingname, dwelltime, sumhamdwellx, sumhamdwelly, sumhamdwellz, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_xyz_sum_one_spin(title, false, savingname, dwelltime2, sumhamdwell2x, sumhamdwell2y, sumhamdwell2z, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_xyz_two_different_spins(title, false, savingname, ptgrid, hamx, hamy, hamz, ham2x, ham2y, ham2z, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_xyz_two_same_spins(title, false, savingname, ptgrid, hamx, hamy, hamz, ceil(Int32, nB1/2), ceil(Int32, noffs/2))
# plot_TF1D(title, false, savingname, hamx, hamy, hamz, ham2x, ham2y, ham2z, ceil(Int32, nB1/2), offs)
# plot_TF2D(title, false, savingname, hamx, hamy, hamz, ham2x, ham2y, ham2z, ceil(Int32, nB1/2), offs)
# plot_avH_B1(title, false, savingname, hamx, hamy, hamz, ham2x, ham2y, ham2z, rfmax, offs)
# spectr_avsum = fft_fid_digits(sumhamdwellz, npoints+2, nB1, noffs)
# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_avsum[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], xlabel="freq (Hz)", legend=false, ylabel="relative intensity (%)", title="FFT - sum(average H 1z) - 0 kHz offs ")
# spectr_av = fft_fid_digits(hamdwellz, npoints+2, nB1, noffs)
# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_av[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], xlabel="freq (Hz)", legend=false, ylabel="relative intensity (%)", title="FFT - average H 1z - 0 kHz offs ")
# contour(offs,fftshift(fftfreq(npoints+2,1/dwell)), spectr_av[ceil(Int32,nB1/2),:,:]', clims=(-100,100), ylims=(-1000, 1000), xlabel="offset (Hz)", ylabel="freq (Hz)", title="hard180 - FFT - sum(average H 1z)", color=palette([:darkblue, :white, :darkred], 21))
# spectr_3d = fft_non_equidistant(ptgrid, hamy, nB1, noffs, npoints)
# plot(fftshift(fftfreq(npoints,1/dwell)),spectr_3d[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], xlabel="freq (Hz)", legend=false, ylabel="relative intensity (%)", title="FFT - average H 1z - 0 kHz offs ")
# contour(offs,fftshift(fftfreq(npoints,1/dwell)), spectr_3d[ceil(Int32,nB1/2),:,:]', clims=(-100,100), xlabel="offset (Hz)", ylabel="freq (Hz)", title="hard180 - FFT - average H 1z", color=palette([:darkblue, :white, :darkred], 21))

# overlaying FFT of FID and FFT of avH
# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_av[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], label="FFT of av H", title="FFT of FID and sum(average H 1z) - 0 kHz offs ", linewidth=1.5)
# plot!(fftshift(fftfreq(npoints,1/dwell)),spectr[ceil(Int32, nB1/2),ceil(Int32, noffs/2),:], xlims=(-100, 5000), label="FFT of FID", xlabel="freq (Hz)", ylabel="intensity", linewidth=1.5)

# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_av[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], label="FFT of av H", title="FFT of FID and sum(average H 1z) - 0 kHz offs ", linewidth=1.5)
# plot!(fftshift(fftfreq(npoints,1/dwell)),spectr[ceil(Int32, nB1/2),ceil(Int32, noffs/2),:], xlims=(-100, 5000), ylims=(-500, 500), label="FFT of FID", xlabel="freq (Hz)", ylabel="intensity", linewidth=1.5)

# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_av[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], label="FFT of av H", title="FFT of FID and sum(average H 1z) - 0 kHz offs ", linewidth=1.5)
# plot!(fftshift(fftfreq(npoints,1/dwell)),spectr[ceil(Int32, nB1/2),ceil(Int32, noffs/2),:].*20, xlims=(-100, 5000), label="FFT of FID, enlarged 20 times", xlabel="freq (Hz)", ylabel="intensity", linewidth=1.5)

# plot(fftshift(fftfreq(npoints+2,1/dwell)),spectr_av[ceil(Int32,nB1/2),ceil(Int32,noffs/2),:], label="FFT of av H", title="FFT of FID and sum(average H 1z) - 0 kHz offs ", linewidth=1.5)
# plot!(fftshift(fftfreq(npoints,1/dwell)),spectr[ceil(Int32, nB1/2),ceil(Int32, noffs/2),:].*20, xlims=(-100, 5000), ylims=(-500,500), label="FFT of FID, enlarged 20 times", xlabel="freq (Hz)", ylabel="intensity", linewidth=1.5)


#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      PLOTTING OF THE RESULTS OF WAUGH CRITERIONS
#------------------------------------------------------------------------------------------------------------------------------------------------

# graphtitle = "GARP-MLEV4"

# Waugh exact
# I1 = [(1+dot(beffp[ioffs,:],beffm[ioffs,:]))/2 for ioffs = 1:noffs]
# I2 = [(1-dot(beffp[ioffs,:],beffm[ioffs,:]))/2 for ioffs = 1:noffs]
# beffy = plot(offs./1000, [dot.(beffp[:,1],beffm[:,1])], yguidefontrotation=270, ylabel="Beffy", legend=false, xlims=(-7,30))
# display(plot(offs./1000, Jeff, title=graphtitle*" - Jeff", legend = false, xlims=(-7,30)))
# plot_exact_waugh(J1, J2, I1, I2, graphtitle, offs)
# plot!(twinx(), offs./1000, abs.(J1), legend=false, xlims=(-7,30), xlabel="offset (kHz)", ylabel="J1")
# plot_beffs(beffp, beffm, graphtitle, offs)
# display(plot(offs./1000, [phip phim], label=["phi plus" "phi minus"], title=graphtitle, xlims=(-7,30)))


# display(plot(offs./1000, abs.(lambda), label="WALTZ-4", xlims=(-15, 15), ylims=(0,1), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))# lambda (Waugh simple)
# display(plot(offs./1000, abs.(lambda_WALTZ4), label="WALTZ-4", xlims=(-15, 15), ylims=(0,0.5), linewidth=2, ylabel=L"\lambda", xlabel=L"\nu_I"*" (kHz)"))
# display(plot!(offs./1000, abs.(lambda_WALTZ8), label="WALTZ-8", xlims=(-15, 15), ylims=(0,0.5), linewidth=2))
# display(plot!(offs./1000, abs.(lambda_WALTZ16), label="WALTZ-16", xlims=(-15, 15), ylims=(0,0.5), linewidth=2))
# savefig("Figures/Waugh/simple_WALTZ")

# lambda_WALTZ16 = lambda

# offs_lambda = get_offset(1001, 80000)

#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      SAVE PULSE SHAPE IN XYT OR XYZT FORMAT
#------------------------------------------------------------------------------------------------------------------------------------------------


# info = String("# M4P5F9\n# rfmax = $rfmax")
# save_file = "C:/Users/Lisa/Documents/KIT/WS202324/Pulse/xyzt_files/M4P5F9_rf10kHz.xyt"
# xyt_save(save_file, info, ux, uy, ut)
# xyzt_save(save_file, info, ux, uy, uz, ut)







#------------------------------------------------------------------------------------------------------------------------------------------------
#                                                      OLD VERSIONS - NOT USED ANYMORE
#------------------------------------------------------------------------------------------------------------------------------------------------



# using CSV
# using DataFrames
# using CurveFit
# using StatsBase
# using Printf

# function read_in_biggest_sidebands(filename)
#     data = DataFrame(CSV.File(filename))
#     Js = zeros(size(data)[1])
#     sidebands = zeros(size(data)[1],size(data)[2]-1)

#     for line in 1:size(data)[1]
#         Js[line] = data[line,1]
#         if !(isnan(data[line,2]))
#             sidebands[line,1] = data[line,2]
#         end
#         if !(isnan(data[line,3]))
#             sidebands[line,2] = data[line,3]
#         end
#     end

#     for line in 1:size(data)[1]
#         if isnan(data[line,2])
#             sidebands = sidebands[setdiff(1:end, line), :]
#             deleteat!(Js, line)
#         end
#     end

#     return Js, sidebands
# end

# function plot_fit_sidebands(Js, sidebands, xlimits, ylimits, offskhz, p0)

#     iJ = size(sidebands)[1]
#     offsstring = "offset " * string(offskhz) * " kHz"
#     scatter(Js, sidebands, label=["onresonant" offsstring], linewidth=2, xlabel="coupling constant "*L"J_{IS}", dpi=500, ylabel="intensity of biggest sideband in %")

#     @. modelside(x, p) = p[1]*x^2
#     fitside = curve_fit(modelside, Js[1:iJ], sidebands[1:iJ,1], p0)
#     fitcurve = [coef(fitside)[1]*x^2 for x in Js[1:iJ]]
#     SStot = sum((sidebands[1:iJ,1].-mean(sidebands[1:iJ,1])).^2)
#     SSres = sum((sidebands[1:iJ,1].-fitcurve).^2)
#     Rsq = 1-SSres/SStot
#     plot!(Js[1:iJ], fitcurve, xlims=xlimits, ylims=ylimits, color=palette(:tab10)[1], label="p="*(@sprintf "%.2E" coef(fitside)[1]), dpi=500, linewidth=2)
#     @show coef(fitside)
#     # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # fit_error = stderror(fitside)

#     @. modelside(x, p) = p[1]*x^2
#     fitside = curve_fit(modelside, Js[1:iJ], sidebands[1:iJ,2], p0)
#     fitcurve = [coef(fitside)[1]*x^2 for x in Js[1:iJ]]
#     SStot = sum((sidebands[1:iJ,2].-mean(sidebands[1:iJ,2])).^2)
#     SSres = sum((sidebands[1:iJ,2].-fitcurve).^2)
#     Rsq = 1-SSres/SStot
#     display(plot!(Js[1:iJ], fitcurve, legend=:topleft, xlims=xlimits, ylims=ylimits, color=palette(:tab10)[2], dpi=500, label="p="*(@sprintf "%.2E" coef(fitside)[1]), linewidth=2))
#     @show coef(fitside)
#     # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # fit_error = stderror(fitside)



#     # scatter(Js, sidebands./maxintens*100, label=["onresonant" "offset 15 kHz"], linewidth=2, dpi=500, xlabel="coupling constant "*L"J_{IS}", ylabel="intensity of biggest sideband in %")

#     # @. modelside(x, p) = p[1]*x^2+p[2]
#     # fitside = curve_fit(modelside, Js[1:iJ], sidebands[1:iJ,1], p0)
#     # fitcurve = [coef(fitside)[1]*x^2+coef(fitside)[2] for x in Js[1:iJ]]
#     # SStot = sum((sidebands[1:iJ,1].-mean(sidebands[1:iJ,1])).^2)
#     # SSres = sum((sidebands[1:iJ,1].-fitcurve).^2)
#     # Rsq = 1-SSres/SStot
#     # plot!(Js[1:iJ], fitcurve./maxintens*100, xlims=(50,xmax), ylims=(ymin,ymax), color=:blue, dpi=500, label=L"R^2="*string(round(Rsq,digits=2)), linewidth=2)
#     # @show coef(fitside)
#     # # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # # fit_error = stderror(fitside)

#     # @. modelside(x, p) = p[1]*x^2+p[2]
#     # fitside = curve_fit(modelside, Js[1:iJ], sidebands[1:iJ,2], p0)
#     # fitcurve = [coef(fitside)[1]*x^2+coef(fitside)[2] for x in Js[1:iJ]]
#     # SStot = sum((sidebands[1:iJ,2].-mean(sidebands[1:iJ,2])).^2)
#     # SSres = sum((sidebands[1:iJ,2].-fitcurve).^2)
#     # Rsq = 1-SSres/SStot
#     # display(plot!(Js[1:iJ], fitcurve./maxintens*100, xlims=(40,xmax), ylims=(ymin,ymax), color=:orange, dpi=500, label=L"R^2="*string(round(Rsq,digits=2)), linewidth=2))
#     # @show coef(fitside)
#     # # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # # fit_error = stderror(fitside)
# end


# function plot_fit_maxintens(Js, maxvals, xlimits, ylimits, offskhz, p0)

#     iJ = size(maxvals)[1]
#     @show iJ

#     offsstring = "offset " * string(offskhz) * " kHz"
    
#     scatter(Js, maxvals, label=["onresonant" offsstring], dpi=500, linewidth=2, xlabel="coupling constant "*L"J_{IS}", ylabel="maximum intensity in %")

#     @. modelmax(x, p) = 100+p[1]*x^2
#     fitmax = curve_fit(modelmax, Js[1:iJ], maxvals[1:iJ,1], p0)
#     fitcurve = [100+coef(fitmax)[1]*x^2 for x in Js[1:iJ]]
#     SStot = sum((maxvals[1:iJ,1].-mean(maxvals[1:iJ,1])).^2)
#     SSres = sum((maxvals[1:iJ,1].-fitcurve).^2)
#     Rsq = 1-SSres/SStot
#     plot!(Js[1:iJ], fitcurve, xlims=xlimits, ylims=ylimits, color=palette(:tab10)[1], dpi=500, label="p="*(@sprintf "%.2E" coef(fitmax)[1]), linewidth=2)
#     @show coef(fitmax)
#     # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # fit_error = stderror(fitside)

#     @. modelmax(x, p) = 100+p[1]*x^2
#     fitmax = curve_fit(modelmax, Js[1:iJ], maxvals[1:iJ,2], p0)
#     fitcurve = [100+coef(fitmax)[1]*x^2 for x in Js[1:iJ]]
#     SStot = sum((maxvals[1:iJ,2].-mean(maxvals[1:iJ,2])).^2)
#     SSres = sum((maxvals[1:iJ,2].-fitcurve).^2)
#     Rsq = 1-SSres/SStot
#     display(plot!(Js[1:iJ], fitcurve, legend=:bottomleft, xlims=xlimits, dpi=500, ylims=ylimits, color=palette(:tab10)[2], label="p="*(@sprintf "%.2E" coef(fitmax)[1]), linewidth=2))
    
#     @show coef(fitmax)
#     # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # fit_error = stderror(fitside)
# end


# pulsename = "MLEV64"
# pulsename = "WALTZ16"
# pulsename = "GARP1"
# pulsename = "M4P5F9"
# pulsename = "STUD-M4P9"
# pulsename = "caWURST2-M4P5"
# pulsename = "BUSS"
# pulsename = "BROCODE_lvl1_linewidth6Hz"
# pulsename = "BROCODE_lvl4_linewidth6Hz"

# Js, sidebands = read_in_biggest_sidebands("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/all_same_rfav/csv_files/Sidebands/"*pulsename*"_sidebands.csv")
# plot_fit_sidebands(allJs, sidebands, (40, 260), (0,3), 15, [0.002])
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Sidebands_fitted_v2.png")
# Js, sidebands = read_in_biggest_sidebands("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/all_same_rfav/Sidebands/csv_files/"*pulsename*"_Sidebands_subtracted.csv")
# plot_fit_sidebands(allJs, sidebands, (40, 260), (0,3), 15, [0.002])
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Sidebands_sub_fitted_v2.png")
# Js, maxintens = read_in_biggest_sidebands("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/all_same_rfav/Sidebands/csv_files/"*pulsename*"_MaxIntens.csv")
# plot_fit_maxintens(Js, maxintens.*2.5./2.505, (40,260), (85,100), 3, [-0.01])
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_max_fitted_v2.png")

# plot_sidebands_stacked_offs(spectr./10, noffs, offs, (0, 800), (-2, 8), "", dwell, npoints)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_sidebands_stacked.png")
# plot_sidebands_stacked_offs(spectr./10, noffs, offs, (0, 2500), (-4, 7), "", dwell, npoints)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_sidebands_stacked_v2.png")
# plot_sidebands_stacked_offs(spectr./10, noffs, offs, (0, 2500), (-4, 7), "", dwell, npoints)
# savefig("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_sidebands_stacked.png")
# display(plot(fftshift(fftfreq(npoints,1/dwell)), allspectr[1,1,1,:], legend=false, xlabel=L"\nu_I"*" (kHz)", ylabel="projection intensity (%)", linewidth=2, ylims=(-10,10), dpi=500))





# textside = "At J,Intensity onresonant (%),Intensity at 15 kHz offset (%)\n"
# textsidesub = "At J,Intensity onresonant (%),Intensity at 15 kHz offset (%)\n"
# textmax = "At J,Intensity onresonant (%),Intensity at 15 kHz offset (%)\n"

# for iJ = 1:size(maxintens[1,:])[1]
#     textside = textside * string(allJs[iJ], ",", maxintensside[1,iJ]/2.505, ",", maxintensside[2,iJ]/2.505, "\n")
#     textsidesub = textsidesub * string(allJs[iJ], ",", maxintenssidesub[1,iJ]/2.505, ",", maxintenssidesub[2,iJ]/2.505, "\n")
#     textmax = textmax * string(allJs[iJ], ",", maxintens[1,iJ]/2.5, ",", maxintens[2,iJ]/2.505, "\n")
# end


# write_to_file("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Sidebands.csv", textside)
# write_to_file("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_Sidebands_subtracted.csv", textsidesub)
# write_to_file("C:/Users/Lisa/Documents/KIT/WS202324/Simulationen/Figures/for_masters_thesis/"*pulsename*"_maxintens.csv", textmax)

# @show maxsideatfreq
# @show maxsideatfreqsub


# data = DataFrame(CSV.File("C:/Users/Lisa/Documents/KIT/WS202324/Experiments/Decoupling_MeOH/powerleveltoamplitude.csv"))

# scatter(collect(data[:,1]),collect(data[:,3])./1000, ylabel="rf amplitude (kHz)", xlabel="power level (W)", label="", dpi=500)

# @. modelpowerlev(x, p) = p[1]*x^0.5
# fitpowerlev = curve_fit(modelpowerlev, collect(data[:,1]), collect(data[:,3]), [1000])
# fitcurve = [coef(fitmax)[1]*x^0.5 for x =1:4]
# plot!(collect(1:80), collect(1:80).^0.5*1.786, color=palette(:tab10)[1], dpi=500, legend=true, label="rf amplitude"*L"=1.786\cdot "*"powerlevel"*L"^{0.5} \cdot"*"powerlevel", linewidth=2)
# @show coef(fitmax)
#     # StatsBase.L2dist(sidebands[1:10,1], fitcurve)
#     # fit_error = stderror(fitside)




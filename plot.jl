using Printf
using Measures
using LaTeXStrings
using Dierckx

"Fourier Transform of FID\n
output: spectr"
function fft_fid(fid, nB1, noffs, npoints)
  spectr = zeros(nB1, noffs, npoints)
  for iB1 = 1:nB1
    for ioffs = 1:noffs
        spectr[iB1, ioffs,:] = real(fftshift(fft(fid[iB1,ioffs,:])))
    end
  end
  return spectr
end

"Fourier Transform of FID with another amount of digits than npoints\n
output: spectr"
function fft_fid_digits(fid, ndigits, nB1, noffs)
  spectr = zeros(nB1, noffs, ndigits)
  for iB1 = 1:nB1
    for ioffs = 1:noffs
        spectr[iB1, ioffs,:] = real(fftshift(fft(fid[iB1,ioffs,:])))
    end
  end
  return spectr
end

"Fourier Transform of FID with non equidistant grid \n
spine interpolation is performed to achieve the equidistant grid\n
output: spectr"
function fft_non_equidistant(time, signal, nB1, noffs, npoints)
  spectr = zeros(nB1, noffs, npoints)
  ti_equi = range(0.0001, stop=time[end]*0.9999, length=npoints)
  for iB1 = 1:nB1
    for ioffs = 1:noffs
      spl = Spline1D(time, signal[iB1, ioffs,:], ti_equi; w=ones(length(time)), k=3, bc="nearest")
      sig_equi = spl(ti_equi)
      spectr[iB1,ioffs,:] = real(fftshift(fft(sig_equi)))
    end
  end
  return spectr
end

"Project spectrum to array containing maximum peak height for each offset\n
output: projection"
function get_projection(spectr, nB1, noffs)
  projection = zeros(size(spectr)[1], size(spectr)[2])
  for iB1 = 1:nB1
    for ioffs = 1:noffs
        projection[iB1, ioffs] = maximum(spectr[iB1, ioffs, :])
    end
  end
  return projection
end

"plot Fourier transformed spectrum of the FID\n
output: -"
function plot_1d(spectr, title, xlimit, ylimit, npoints, dwell)
  plot(fftshift(fftfreq(npoints,1/dwell)),spectr, xlabel=L"\nu"*" (Hz)", legend=false, ylabel="intensity (%)", title=(title), xlim=xlimit, ylim=ylimit)
end

"plot single FID \n
output: -"
function plot_fid(time, fid)
  plot(time,real(fid), thickness_scaling=1.2, legend=false, xlabel="time (s)", ylabel=("x magnetization"))
  #plot(time,real(fid), thickness_scaling=2, legend=false, GradientDescent=false, xlabel="time (s)", ylabel=("relative intensity (%)"))
end

" ### plot_projection(title, ifsave, savingname, projection, offs)
plot the projection and save it, if \"ifsave\" is true"
function plot_projection(title, ifsave, savingname, projection, offs)
  plot(offs./1000, projection, legend=false)
  xlabel!(L"\nu_S"*" (kHz)")
  ylabel!("projection intensity (%)")
  display(title!(title))
  if ifsave == true
    savefig(savingname*"_projection.png")
  end
end

"plot phase and amplitude of Bruker files\n
output: -"
function plot_pulse_bruker(filename, tpulse)
  ampl,phase = read_bruker_wave(filename)
  ampl = ampl*rfmax/100.0
  ptime = collect(tpulse/size(ampl)[1]:tpulse/size(ampl)[1]:tpulse)

  plot(ptime, phase, label="xy-Phase", legend=:topleft, ylabel="Phase (°)", xlabel="Time (s)")
  plot!(twinx(), ptime, ampl, label= "rf-amplitude", color=:red, ylabel="Amplitude (Hz)")
end

"plot phase and amplitude of pulse components ux, uy, ut\n
output: -"
function plot_ampl_phase(ux, uy, ut, title)
  ampl, phase = xy_to_ampl_phase(ux, uy)
  @show ampl, phase
  ptime = get_ptgrid(ut)
  pushfirst!(ptime, 0); pushfirst!(ampl,ampl[1]); pushfirst!(phase,phase[1])
  plot(ptime.*1000, phase, label="xy-Phase", legend=:topleft, linetype=:steppre, fill=(0,0.5), ylabel="Phase (°)", xlabel="Time (ms)")
  plot!(twinx(), ptime.*1000, ampl, label= "rf-amplitude", color=:red, legend=:topright, ylabel="Amplitude (Hz)", title=title)
end

"plot phase and amplitude of pulse components ux, uy, ut\n
output: -"
function plot_ampl_phase_2(acqux, acquy, acqut)
  ampl, phase = xy_to_ampl_phase(acqux[2], acquy[2])
  ptime = get_ptgrid(acqut[2])
  plot(ptime.*1000, phase, label="xy-Phase", ylabel="Phase (°)", linetype=:steppre, fill=(0,0.2), xlabel="Time (ms)", linewidth=2, dpi=500)
  invis = ones(size(ptime)[1]).*(-100)
  plot!(ptime, invis, label= "rf-amplitude", color=:darkorange, ylims=(-5,365), linewidth=2, legend=:top, dpi=500)
  plot!(twinx(), ptime.*1000, ampl./1000, legend=false, label="Phase (°)", color=:darkorange, ylims=(0, 14), ylabel="Amplitude (kHz)", linewidth=2, dpi=500)
end

"plot only phase of the pulse components ux, uy, ut\n
output: -"
function plot_pulse_phase(ux, uy, ut)
  ampl, phase = xy_to_ampl_phase(ux, uy)
  ptime = get_ptgrid(ut)
  plot(ptime, phase, label="xy-Phase", legend=:topleft, ylabel="Phase (°)", xlabel="Time (s)")
end

"plot 4 plots containing the phase, amplitude, offset and contour plot\n
output: -"
function plot_pulse_phase_ampl_offs_contour(ux, uy, uz, ut, spectr, npoints, dwell, offs)

  ampl, phase = xy_to_ampl_phase(ux, uy)
  ptime = get_ptgrid(ut)

  plot_phase = plot(ptime, phase, label="xy-Phase", legend=:topleft, ylabel="Phase (°)", xlabel="Time (s)")
  plot_ampl = plot(ptime, ampl, label= "rf-amplitude", color=:red, xlabel="Time (s)", ylabel="Amplitude (Hz)")
  plot_offs = plot(ptime, uz, label="offset", color=:green, xlabel="Time (s)", ylabel="Offset (Hz)")
  plot_contour = contour(offs,fftshift(fftfreq(npoints,1/dwell)), spectr', xlabel="offset (Hz)", ylabel="freq (Hz)")
  plot(plot_phase, plot_ampl, plot_offs, plot_contour, layout=4)
  
end

"plot ux and uy components of pulse\n
output: -"
function plot_ux_uy(plot_x, plot_y, plot_t, xlimit, ylimit, title)
  plot_x = copy(plot_x)
  plot_y = copy(plot_y)
  plot_t = copy(plot_t)
  ptime = [0.0]
  for dig in eachindex(plot_t)
    append!(ptime, ptime[dig]+plot_t[dig])
  end
  pushfirst!(plot_x, plot_x[1])
  pushfirst!(plot_y, plot_y[1])
  display(plot(ptime.*1000, [plot_x./1000 plot_y./1000], label=["ux" "uy"], linetype=:steppre, fill=(0, 0.5), xlabel="time (ms)", ylabel="amplitude (kHz)", ylim=ylimit, xlim=xlimit, title=title))
  
end

"plot ux, uy and ut components of pulse\n
output: -"
function plot_ux_uy_uz(plot_x, plot_y, plot_z, plot_t)
  plot_x = copy(plot_x)
  plot_y = copy(plot_y)
  plot_z = copy(plot_z)
  plot_t = copy(plot_t)
  ptime = [0.0]
  for dig in eachindex(plot_t)
    append!(ptime, ptime[dig]+plot_t[dig])
  end
  pushfirst!(plot_x, plot_x[1])
  pushfirst!(plot_y, plot_y[1])
  pushfirst!(plot_z, plot_z[1])
  plot(ptime.*1000, [plot_x./1000 plot_y./1000], legend=false, label=["ux" "uy"], linetype=:steppre, fill=(0, 0.5), ylabel="amplitude (kHz)")
  invis = zeros(size(ptime)[1])
  display(plot!(twinx(), ptime.*1000, [invis, invis, plot_z./1000],label=["ux" "uy" "uz"],  fill=(0, 0.5), ylims=(minimum(plot_z)/1000, maximum(plot_z)/1000), linetype=:steppre, ylabel="offset (kHz)", xlabel="time (ms)"))  
end


"functions for average Hamiltonian\n
output: -"
function plot_xyz_one_spin(title, ifsave, savingname, tim, hamx, hamy, hamz, iB1, ioffs)
  plotx = plot(tim.*1000, hamx[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{x}", ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yguidefontrotation=270, title=title)
  ploty = plot(tim.*1000, hamy[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{y}", ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm,legend=false, yguidefontrotation=270, )
  plotz = plot(tim.*1000, hamz[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{z}", ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, yguidefontrotation=270, xaxis="time (ms)")
  display(plot(plotx, ploty, plotz, layout=(3,1), left_margin=5mm))
  if ifsave==true
    savefig(savingname*"_avH_1spin.png")
  end
end


function plot_xyz_sum_one_spin(title, ifsave, savingname, tim, hamx, hamy, hamz, iB1, ioffs)
  plotx = plot(tim.*1000, hamx[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{x}", fill=(0, 0.5), bottom_margin=-3mm, legend=false, yguidefontrotation=270, title=title)
  ploty = plot(tim.*1000, hamy[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{y}", fill=(0, 0.5), bottom_margin=-3mm,legend=false, yguidefontrotation=270, )
  plotz = plot(tim.*1000, hamz[iB1,ioffs,:], yaxis=L"\mathrm{I}_\mathrm{z}", fill=(0, 0.5), legend=false, yguidefontrotation=270, xaxis="time (ms)")
  display(plot(plotx, ploty, plotz, layout=(3,1), left_margin=5mm))
  if ifsave==true
    savefig(savingname*"_sum_avH_1spin.png")
  end
end

function plot_xyz_two_same_spins(title, ifsave, savingname, tim, hamx, hamy, hamz, iB1, ioffs)
  plotx = plot(tim.*1000, hamx[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{I}_\mathrm{x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  ploty = plot(tim.*1000, hamy[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{I}_\mathrm{y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, title=title)
  plotz = plot(tim.*1000, hamz[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{I}_\mathrm{z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotxx = plot(tim.*1000, hamx[iB1,ioffs,:].*hamx[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{1x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  plotyy = plot(tim.*1000, hamy[iB1,ioffs,:].*hamy[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{1y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotzz = plot(tim.*1000, hamz[iB1,ioffs,:].*hamz[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{1z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotxy = plot(tim.*1000, hamx[iB1,ioffs,:].*hamy[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{1y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)", yaxis="average H")
  plotxz = plot(tim.*1000, hamx[iB1,ioffs,:].*hamz[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{1z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)")
  plotyz = plot(tim.*1000, hamz[iB1,ioffs,:].*hamy[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{1z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)")
  display(plot(plotx, ploty, plotz, plotxx, plotyy, plotzz, plotxy, plotxz, plotyz, layout=(3,3), size=(800,700)))
  if ifsave==true
    savefig(savingname*"_avH_2spins.png")
  end
end

function plot_xyz_two_different_spins(title, ifsave, savingname, tim, ham1x, ham1y, ham1z, ham2x, ham2y, ham2z, iB1, ioffs)
  plot1x = plot(tim.*1000, ham1x[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{1x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  plot1y = plot(tim.*1000, ham1y[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{1y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plot1z = plot(tim.*1000, ham1z[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{1z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  plot2x = plot(tim.*1000, ham2x[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plot2y = plot(tim.*1000, ham2y[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)", yaxis="average H")
  plot2z = plot(tim.*1000, ham2z[iB1,ioffs,:], annotations=(tim[end]*0.9*1000, 0.8, L"\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)")
  display(plot(plot1x, plot1y, plot1z, plot2x, plot2y, plot2z, layout=(3,2), size=(700,600), left_margin=7mm, plot_title=title))
  if ifsave==true
    savefig(savingname*"_1op.png")
  end
  plotxx = plot(tim.*1000, ham1x[iB1,ioffs,:].*ham2x[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  plotxy = plot(tim.*1000, ham1x[iB1,ioffs,:].*ham2y[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, title=title)
  plotxz = plot(tim.*1000, ham1x[iB1,ioffs,:].*ham2z[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotyx = plot(tim.*1000, ham1y[iB1,ioffs,:].*ham2x[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false, yaxis="average H")
  plotyy = plot(tim.*1000, ham1y[iB1,ioffs,:].*ham2y[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotyz = plot(tim.*1000, ham1y[iB1,ioffs,:].*ham2z[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), bottom_margin=-3mm, legend=false)
  plotzx = plot(tim.*1000, ham1z[iB1,ioffs,:].*ham2x[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)", yaxis="average H")
  plotzy = plot(tim.*1000, ham1z[iB1,ioffs,:].*ham2y[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)")
  plotzz = plot(tim.*1000, ham1z[iB1,ioffs,:].*ham2z[iB1,ioffs,:], annotations=(tim[end]*0.8*1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1], fill=(0, 0.5), legend=false, xaxis="time (ms)")
  display(plot(plotxx, plotxy, plotxz, plotyx, plotyy, plotyz, plotzx, plotzy, plotzz, layout=(3,3), size=(800,700)))
  if ifsave==true
    savefig(savingname*"_2ops.png")
  end
end

function plot_TF1D(title, ifsave, savingname, ham1x, ham1y, ham1z, ham2x, ham2y, ham2z, iB1, offs)
  hamoff1x = [sum(ham1x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1y = [sum(ham1y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1z = [sum(ham1z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff2x = [sum(ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff2y = [sum(ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff2z = [sum(ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  plot1x = plot(offs./1000, hamoff1x, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{1x}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false, yaxis="average H")
  plot1y = plot(offs./1000, hamoff1y, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{1y}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false)
  plot1z = plot(offs./1000, hamoff1z, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{1z}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false, yaxis="average H")
  plot2x = plot(offs./1000, hamoff2x, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false)
  plot2y = plot(offs./1000, hamoff2y, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, legend=false, xaxis="offset (kHz)", yaxis="average H")
  plot2z = plot(offs./1000, hamoff2z, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, legend=false, xaxis="offset (kHz)")
  display(plot(plot1x, plot1y, plot1z, plot2x, plot2y, plot2z, layout=(3,2), size=(700,600), left_margin=7mm, plot_title=title)) 
  if ifsave==true
    savefig(savingname*"_TF1D_1op.png")
  end
  hamoff1x2x = [sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1x2y = [sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1x2z = [sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1y2x = [sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1y2y = [sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1y2z = [sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1z2x = [sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1z2y = [sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  hamoff1z2z = [sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs]
  plotxx = plot(offs./1000, hamoff1x2x, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false, yaxis="average H")
  plotxy = plot(offs./1000, hamoff1x2y, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false, title=title)
  plotxz = plot(offs./1000, hamoff1x2z, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false)
  plotyx = plot(offs./1000, hamoff1y2x, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false, yaxis="average H")
  plotyy = plot(offs./1000, hamoff1y2y, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false)
  plotyz = plot(offs./1000, hamoff1y2z, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, bottom_margin=-3mm, legend=false)
  plotzx = plot(offs./1000, hamoff1z2x, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2x}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, legend=false, xaxis="offset (kHz)", yaxis="average H")
  plotzy = plot(offs./1000, hamoff1z2y, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2y}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, legend=false, xaxis="offset (kHz)")
  plotzz = plot(offs./1000, hamoff1z2z, annotations=(offs[end]*0.7/1000, 0.8, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2z}"), ylim=(-1, 1), yticks=[-1, 0, 1],  fill=(0, 0.5), linewidth=2, legend=false, xaxis="offset (kHz)")
  display(plot(plotxx, plotxy, plotxz, plotyx, plotyy, plotyz, plotzx, plotzy, plotzz, layout=(3,3), size=(800,700)))
  if ifsave==true
    savefig(savingname*"_TF1D_2ops.png")
  end
end

function plot_TF2D(title, ifsave, savingname, ham1x, ham1y, ham1z, ham2x, ham2y, ham2z, iB1, offs)
  hamoff1x2x = reduce(hcat,[[sum(ham1x[iB1,i2offs,:].*ham2x[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1y2y = reduce(hcat,[[sum(ham1y[iB1,i2offs,:].*ham2y[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1z2z = reduce(hcat,[[sum(ham1z[iB1,i2offs,:].*ham2z[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1x2y = reduce(hcat,[[sum(ham1x[iB1,i2offs,:].*ham2y[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1x2z = reduce(hcat,[[sum(ham1x[iB1,i2offs,:].*ham2z[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1y2x = reduce(hcat,[[sum(ham1y[iB1,i2offs,:].*ham2x[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1y2z = reduce(hcat,[[sum(ham1y[iB1,i2offs,:].*ham2z[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1z2x = reduce(hcat,[[sum(ham1z[iB1,i2offs,:].*ham2x[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  hamoff1z2y = reduce(hcat,[[sum(ham1z[iB1,i2offs,:].*ham2y[iB1,i1offs,:])/size(ham1x)[3] for i1offs = 1:noffs] for i2offs = 1:noffs])'
  plotxx = heatmap(offs./1000, offs./1000, hamoff1x2x, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="offset (kHz)")
  plotxy = heatmap(offs./1000, offs./1000, hamoff1x2y, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, title=title)
  plotxz = heatmap(offs./1000, offs./1000, hamoff1x2z, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plotyx = heatmap(offs./1000, offs./1000, hamoff1y2x, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="offset (kHz)")
  plotyy = heatmap(offs./1000, offs./1000, hamoff1y2y, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plotyz = heatmap(offs./1000, offs./1000, hamoff1y2z, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plotzx = heatmap(offs./1000, offs./1000, hamoff1z2x, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)", yaxis="offset (kHz)")
  plotzy = heatmap(offs./1000, offs./1000, hamoff1z2y, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)")
  plotzz = heatmap(offs./1000, offs./1000, hamoff1z2z, annotations=(offs[end]*0.65/1000, offs[end]*0.8/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)")
  bar = scatter([0,0], [0,1], zcolor=[0,3], color=palette([:firebrick, :white, :royalblue4], 21),  clims=(-1,1), xlims=(1,1.1), axis=false, label="", grid=false, left_margin = -10mm)
  display(plot(plotxx, plotxy, plotxz, bar, plotyx, plotyy, plotyz, bar, plotzx, plotzy, plotzz, bar, layout=(grid(3,4, widths=(0.3, 0.3, 0.3, 0.07))), size=(800,700), left_margin=2mm))
  if ifsave==true
    savefig(savingname*"_TF2D.png")
  end
end

"plot of average Hamiltonian dependent on B1 field inhomogeneity and offset\n
output: -"
function plot_avH_B1(title, ifsave, savingname, ham1x, ham1y, ham1z, ham2x, ham2y, ham2z, rfmax, offs)
  fB1 = get_B1()
  B1 = fB1 .* rfmax
  hamoff1x = reduce(hcat,[[sum(ham1x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1y = reduce(hcat,[[sum(ham1y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1z = reduce(hcat,[[sum(ham1z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff2x = reduce(hcat,[[sum(ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff2y = reduce(hcat,[[sum(ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff2z = reduce(hcat,[[sum(ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  plot1x = heatmap(offs./1000, B1./1000, hamoff1x', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{1x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="B1 field (kHz)")
  plot1y = heatmap(offs./1000, B1./1000, hamoff1y', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{1y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plot1z = heatmap(offs./1000, B1./1000, hamoff1z', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{1z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="B1 field (kHz)")
  plot2x = heatmap(offs./1000, B1./1000, hamoff2x', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plot2y = heatmap(offs./1000, B1./1000, hamoff2y', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)", yaxis="B1 field (kHz)")
  plot2z = heatmap(offs./1000, B1./1000, hamoff2z', annotations=(offs[end]*0.8/1000, (B1[end]-0.08*(B1[end]-B1[1]))/1000, L"\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)")
  bar = scatter([0,0], [0,1], zcolor=[0,3], color=palette([:firebrick, :white, :royalblue4], 21),  clims=(-1,1), xlims=(1,1.1), axis=false, label="", grid=false)
  display(plot(plot1x, plot1y, bar, plot1z, plot2x, bar, plot2y, plot2z, bar, layout=(grid(3,3, widths=(0.45, 0.45, 0.07))), size=(600,650), plot_title=title, left_margin=7mm))
  if ifsave==true
    savefig(savingname*"_avH_B1_1op.png")
  end
  hamoff1x2x = reduce(hcat,[[sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1x2y = reduce(hcat,[[sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1x2z = reduce(hcat,[[sum(ham1x[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1y2x = reduce(hcat,[[sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1y2y = reduce(hcat,[[sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1y2z = reduce(hcat,[[sum(ham1y[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1z2x = reduce(hcat,[[sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2x[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1z2y = reduce(hcat,[[sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2y[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  hamoff1z2z = reduce(hcat,[[sum(ham1z[iB1,ceil(Int32, noffs/2),:].*ham2z[iB1,ioffs,:])/size(ham1x)[3] for ioffs = 1:noffs] for iB1 = 1:nB1])'
  plotxx = heatmap(offs./1000, B1./1000, hamoff1x2x, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="B1 field (kHz)")
  plotxy = heatmap(offs./1000, B1./1000, hamoff1x2y, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plotxz = heatmap(offs./1000, B1./1000, hamoff1x2z, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1x}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm)
  plotyx = heatmap(offs./1000, B1./1000, hamoff1y2x, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, yaxis="B1 field (kHz)")
  plotyy = heatmap(offs./1000, B1./1000, hamoff1y2y, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, )
  plotyz = heatmap(offs./1000, B1./1000, hamoff1y2z, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1y}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, bottom_margin=-5mm, )
  plotzx = heatmap(offs./1000, B1./1000, hamoff1z2x, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2x}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)", yaxis="B1 field (kHz)")
  plotzy = heatmap(offs./1000, B1./1000, hamoff1z2y, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2y}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)")
  plotzz = heatmap(offs./1000, B1./1000, hamoff1z2z, annotations=(offs[end]*0.65/1000, (B1[end]-0.1*(B1[end]-B1[1]))/1000, L"\mathrm{2I}_\mathrm{1z}\mathrm{I}_\mathrm{2z}"), clim=(-1, 1), color=palette([:firebrick, :white, :royalblue4], 21), cbar=false, left_margin=-3mm, xaxis="offset (kHz)")
  display(plot(plotxx, plotxy, plotxz, bar, plotyx, plotyy, plotyz, bar, plotzx, plotzy, plotzz, bar, layout=(grid(3,4, widths=(0.3, 0.3, 0.3, 0.07))), size=(800,700), plot_title=title, left_margin=2mm))
  if ifsave==true
    savefig(savingname*"_avH_B1_2ops.png")
  end

end


"plot exact Waugh criterion\n
output: -"
function plot_exact_waugh(J1, J2, I1, I2, graphtitle, offs)
  J1plot = plot(offs./1000, abs.(J1), yguidefontrotation=270, legend=false, ylabel="J1", title=graphtitle, xlims=(-7,30))
  J2plot = plot(offs./1000, abs.(J2), yguidefontrotation=270, legend=false, ylabel="J2", xlims=(-7,30))
  I1plot = plot(offs./1000, abs.(I1), yguidefontrotation=270, legend=false, ylabel="I1", xlims=(-7,30))
  I2plot = plot(offs./1000, abs.(I2), yguidefontrotation=270, legend=false, ylabel="I2", xlims=(-7,30))
  display(plot(J1plot, I1plot, J2plot, I2plot, layout=(grid(4,1)), left_margin=7mm, size=(600,600)))
end

"plot effective B-field calculated during Waugh calculation\n
output: -"
function plot_beffs(beffp, beffm, graphtitle, offs)
  beffx = plot(offs./1000, [beffp[:,1] beffm[:,1]], yguidefontrotation=270, ylabel="Beffx", legend=false, title=graphtitle*" Beff", xlims=(-7,30))
  beffy = plot(offs./1000, [beffp[:,2] beffm[:,2]], yguidefontrotation=270, ylabel="Beffy", legend=false, xlims=(-7,30))
  beffz = plot(offs./1000, [beffp[:,3] beffm[:,3]], yguidefontrotation=270, ylabel="Beffz", label=["plus" "minus"], xlims=(-7,30))
  display(plot(beffx, beffy, beffz, layout=(grid(3,1)), left_margin=7mm))
end

"plot stacked sidebands calculated at different offsets but at the same coupling strenght J \n
only spectra with 90% of maximum intensity are taken into account\n
output: -"
function plot_sidebands_stacked_offs(spectr, noffs, offs, xlimits, ylimits, title, dwell, npoints)
  spectrdecoupled = Float64[]
  used_offs = []
  maxintens = maximum(spectr[ceil(Int32,nB1/2),:,:])
  firstflag = false
  for ioffs = 1:noffs
      if maximum(spectr[ceil(Int32,nB1/2), ioffs, :]) > (0.9*maxintens)
          if firstflag == false
              spectrdecoupled = spectr[ceil(Int32,nB1/2), ioffs, :]
              firstflag = true
          else
              spectrdecoupled = hcat(spectrdecoupled, spectr[ceil(Int32,nB1/2), ioffs, :])
          end
          append!(used_offs, offs[ioffs])
      end
  end

  @show used_offs
  # CList = reshape( range(colorant"blue", stop=colorant"blue",length=size(spectrdecoupled)[2]), 1, size(spectrdecoupled)[2])
  # display(plot(fftshift(fftfreq(npoints,1/dwell)), reverse(spectrdecoupled), label=reverse(used_offs'), legend=false, xlims=xlimits, ylims=ylimits, title=title, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)", linecolor=CList))
  display(plot(fftshift(fftfreq(npoints,1/dwell)), reverse(spectrdecoupled), label=reverse(used_offs'), legend=false, xlims=xlimits, ylims=ylimits, title=title, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)", color=:blue))

  
end

"plot stacked sidebands calculated at different coupling strengths J but at the same offset\n
output: -"
function plot_sidebands_stacked_J(allspectr, xlimits, ylimits, iB1, ioffs, allJs, title, dwell, npoints)

  CList = reshape( range(colorant"blue", stop=colorant"red",length=size(allJs)[1]), 1, size(allJs)[1])
  display(plot(fftshift(fftfreq(npoints,1/dwell)), reverse(allspectr[:,iB1, ioffs,:]'), label=reverse(allJs'), xlims=xlimits, ylims=ylimits, title=title, xlabel="freq (Hz)", ylabel="intensity", linecolor=CList))

  
end
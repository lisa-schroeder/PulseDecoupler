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
  display(plot(fftshift(fftfreq(npoints,1/dwell)), reverse(spectrdecoupled), label=reverse(used_offs'), legend=false, xlims=xlimits, ylims=ylimits, title=title, xlabel=L"\nu_S"*" (Hz)", ylabel="intensity (%)", color=:blue)) 
end


"plot stacked sidebands calculated at different coupling strengths J but at the same offset\n
output: -"
function plot_sidebands_stacked_J(allspectr, xlimits, ylimits, iB1, ioffs, allJs, title, dwell, npoints)
  CList = reshape( range(colorant"blue", stop=colorant"red",length=size(allJs)[1]), 1, size(allJs)[1])
  display(plot(fftshift(fftfreq(npoints,1/dwell)), reverse(allspectr[:,iB1, ioffs,:]'), label=reverse(allJs'), xlims=xlimits, ylims=ylimits, title=title, xlabel="freq (Hz)", ylabel="intensity", linecolor=CList))  
end


"### ampl_ phase_ to_xy(ampl, phase)
convert amplitude and phase to ux and uy"
function ampl_phase_to_xy(ampl, phase)
    x = zeros(length(ampl))
    y = zeros(length(ampl))
    for dig in eachindex(ampl)
        sin_d, cos_d = sincosd(phase[dig]) 
        x[dig] = ampl[dig]*cos_d
        y[dig] = ampl[dig]*sin_d
    end
    return x, y
end


"convert pulse components ux, uy and ut into amplitude and phase\n
output: ampl, phase"
function xy_to_ampl_phase(ux, uy)
    phase = Float64[]
    ampl = Float64[]
  
    for dig in eachindex(ux)
      append!(ampl, sqrt(ux[dig]^2+uy[dig]^2))
      if ux[dig] == 0 && uy[dig] == 0
        append!(phase, 0)
      elseif ux[dig] >= 0
        append!(phase, mod(atand(uy[dig]/ux[dig]), 360))
      else
        append!(phase, mod(180+atand(uy[dig]/ux[dig]), 360))
      end
    end
  
    return ampl, phase
end


"convert ux, uy and ut to angle and phase\n
output: angle, phase"
function xyt_to_angle_phase(ux, uy, ut)
    phase = Float64[]
    angle = Float64[]
  
    for dig in eachindex(ux)
      append!(angle, ut[dig]*sqrt(ux[dig]^2+uy[dig]^2)*360)
      if ux[dig] >= 0
        append!(phase, mod(atand(uy[dig]/ux[dig]), 360))
      else
        append!(phase, mod(180+atand(uy[dig]/ux[dig]), 360))
      end
    end
  
    return angle, phase
end

"convert arrays of angles and phases to arrays ux, uy and ut\n
output: ux, uy, ut"
function angle_phase_to_xyt(angle, phase, rfmax)

    ux = zeros(length(angle))
    uy = zeros(length(angle))
    ut = zeros(length(angle))
    for digit in eachindex(angle)
        ut[digit] = 1/(rfmax*360)*angle[digit]
    end

    ampl = ones(length(angle))
    ampl = ampl*rfmax
    ux, uy = ampl_phase_to_xy(ampl, phase)

    return ux, uy, ut
end

"average Hamiltonian of pulse is simulated\n
output: hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime"
function simulate_average_H(ux, uy, uz, ut, centered_inversion, haminit, thresh, fB1, nB1, noffs, dwell)

    ampl, phase, ptgrid, tdigit, digdwell = get_final_pulse_non_rep_avH(ux, uy, uz, ut, thresh, noffs, dwell)
    # println("starting calculation of average Hamiltonian")
    hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, dwelltime = get_average_H(ampl, phase, tdigit, centered_inversion, haminit, digdwell, fB1, nB1, noffs, dwell)
    pushfirst!(ptgrid, 0.0)

    return hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime

end



"get repeated pulse necessary for calculating the average Hamiltonian of a repeated pulse, until the aquisition is over, maximum angle must not be larger than threshold (e.g. 4°)\n
output: ampl, phase, ptgrid, tdigit, digdwell"
function get_final_pulse_rep_avH(ux, uy, uz, ut, thresh, npoints, dwell)

    # repeat pulse until it is longer than acquisition time
    ux, uy, uz, ut = repeat_pulse(ux, uy, uz, ut, npoints, dwell)

    # split up pulses at each dwellpoint
    ux, uy, uz, ut = correct_xyz_to_max_angle(ux, uy, uz, ut, thresh)
    ux, uy, uz, ut, digdwell = correct_xyz_to_dwell(ux, uy, uz, ut, npoints, dwell)
    ptgrid = get_ptgrid(ut)
    ampl, phase = xy_to_ampl_phase(ux, uy)

    tdigit = [ptgrid[1]]
    for dig in 2:size(ptgrid)[1]
        append!(tdigit, ptgrid[dig]-ptgrid[dig-1])
    end
    
    return ampl, phase, ptgrid, tdigit, digdwell
end

"get repeated pulse necessary for calculating the average Hamiltonian, maximum angle must not be larger than threshold (e.g. 4°)\n
output: ampl, phase, ptgrid, tdigit, digdwell"
function get_final_pulse_non_rep_avH(ux, uy, uz, ut, thresh, npoints, dwell)
    # split up pulses at each dwellpoint
    ux, uy, uz, ut = correct_xyz_to_max_angle(ux, uy, uz, ut, thresh)
    ux, uy, uz, ut, digdwell = correct_xyz_to_dwell(ux, uy, uz, ut, npoints, dwell)
    ptgrid = get_ptgrid(ut)
    ampl, phase = xy_to_ampl_phase(ux, uy)     # TODO not working for ux=0 and uy=0 at the same time

    # simulation parameters
    tdigit = [ptgrid[1]]
    for dig in 2:size(ptgrid)[1]
        append!(tdigit, ptgrid[dig]-ptgrid[dig-1])
    end
    
    return ampl, phase, ptgrid, tdigit, digdwell
end


"maximum rotation angle must not be larger than threshold (e.g. 4°)\n
output: ux, uy, uz, ut"
function correct_xyz_to_max_angle(ux, uy, uz, ut, thresh)

    fux = Float64[]          
    fuy = Float64[]          
    fuz = Float64[]          
    fut = Float64[]           # final time of each pulses used for simulation

    for dig in eachindex(ux)
        angle = ut[dig]*sqrt(ux[dig]^2+uy[dig]^2)*360
        if angle > thresh
            for rep = 1:Int(floor(angle/thresh))
                append!(fux, ux[dig])
                append!(fuy, uy[dig])
                append!(fuz, uz[dig])
                append!(fut, ut[dig]*thresh/angle)
            end
            append!(fux, ux[dig])
            append!(fuy, uy[dig])
            append!(fuz, uz[dig])
            append!(fut, ut[dig]*mod(angle, thresh)/angle)
        else
            append!(fux, ux[dig])
            append!(fuy, uy[dig])
            append!(fuz, uz[dig])
            append!(fut, ut[dig])
        end
    end
    return fux, fuy, fuz, fut
    
end

"calculate average Hamiltonian which already has corrected rotation angles\n
output: hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, dwelltime"
function get_average_H(ampl, phase, tdigit, centered_inversion, haminit, digdwell, fB1, nB1, noffs, dwell)
    ndigits = size(ampl)[1]

    # array definitions
    anglex = zeros(ndigits)
    angley = zeros(ndigits)
    offs = zeros(noffs)
    w = zeros(3)
    hamx = zeros(nB1,noffs,ndigits+1)
    hamy = zeros(nB1,noffs,ndigits+1)
    hamz = zeros(nB1,noffs,ndigits+1)
    hamdwellx = zeros(nB1,noffs,size(digdwell)[1]+1)              # save only at dwell points
    hamdwelly = zeros(nB1,noffs,size(digdwell)[1]+1)
    hamdwellz = zeros(nB1,noffs,size(digdwell)[1]+1)
    dwelltime = zeros(size(digdwell)[1])
    #hamoffx = zeros(nB1,noffs)
    #hamoffy = zeros(nB1,noffs)
    #hamoffz = zeros(nB1,noffs)
    hamx[:,:,1] = hamx[:,:,1].+haminit[1]
    hamy[:,:,1] = hamy[:,:,1].+haminit[2]
    hamz[:,:,1] = hamz[:,:,1].+haminit[3]
    hamdwellx[1] = haminit[1]
    hamdwelly[1] = haminit[2]
    hamdwellz[1] = haminit[3]

    # simulate offset dependence of pulse/pulse sequence
    for ip in eachindex(ampl)
        phasesin, phasecos = sincosd(phase[ip])
        anglex[ip] = 2*pi*ampl[ip]*phasecos*tdigit[ip]           # rotation angle of this digit around x
        angley[ip] = 2*pi*ampl[ip]*phasesin*tdigit[ip]
    end
    for iB1 = 1:nB1

        for ioffs=1:noffs
            global ham = haminit
            global fram = RotationVec(0,0,0)
            dwellcounter = 1
            digdwellcounter = 1
            for ip in eachindex(ampl)
                w[1] = fB1[iB1]*anglex[ip]
                w[2] = fB1[iB1]*angley[ip]
                w[3] = -2*pi*offs[ioffs]*tdigit[ip]         # 2*pi*offs*tdigit
                wf = fram * w
                local rot=RotationVec(wf[1],wf[2],wf[3])
                global ham = rot'*ham
                global fram = rot'*fram
                hamx[iB1,ioffs,ip+1] = ham[1]
                hamy[iB1,ioffs,ip+1] = ham[2]
                hamz[iB1,ioffs,ip+1] = ham[3]
                if digdwellcounter == digdwell[dwellcounter]
                    hamdwellx[iB1,ioffs,dwellcounter+1] = ham[1]
                    hamdwelly[iB1,ioffs,dwellcounter+1] = ham[2]
                    hamdwellz[iB1,ioffs,dwellcounter+1] = ham[3]
                    dwelltime[dwellcounter] = dwell*(dwellcounter-1)
                    dwellcounter += 1
                    digdwellcounter = 0
                end
                digdwellcounter += 1
                #hamoffx[iB1,ioffs]=hamoffx[iB1,ioffs]+ham[1]/length(ampl)
                #hamoffy[iB1,ioffs]=hamoffy[iB1,ioffs]+ham[2]/length(ampl)
                #hamoffz[iB1,ioffs]=hamoffz[iB1,ioffs]+ham[3]/length(ampl)
                if centered_inversion && ip==ceil(ndigits/2)
                    ham = -ham
                end
            end

        end

    end
    return hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, dwelltime
end

"simulate average Hamiltonian where the pulse is repeated until the aquisiton is finished\n
output: hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime, dtgrid"
function simulate_rep_average_H(ux, uy, uz, ut, centered_inversion, haminit, thresh, npoints, dwell, fB1, nB1, noffs)

    ampl, phase, ptgrid, tdigit, digdwell, dtgrid = get_final_pulse_rep_avH_changed_average(ux, uy, uz, ut, thresh, npoints, dwell)

    hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, dwelltime = get_average_H(ampl, phase, tdigit, centered_inversion, haminit, digdwell, fB1, nB1, noffs, dwell)
    pushfirst!(ptgrid, 0.0)
    append!(dwelltime, dwell*npoints)

    return hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime, dtgrid

end

"generate pulse to perform average Hamiltonian calculation with repeated pulse\n
output: ampl, phase, ptgrid, tdigit, digdwell, dtgrid"
function get_final_pulse_rep_avH_changed_average(ux, uy, uz, ut, thresh, npoints, dwell)

    # repeat pulse until it is longer than acquisition time
    ux, uy, uz, ut = repeat_pulse(ux, uy, uz, ut, npoints, dwell)

    # split up pulses at each dwellpoint
    ux, uy, uz, ut = correct_xyz_to_max_angle(ux, uy, uz, ut, thresh)
    ux, uy, uz, ut, digdwell, dtgrid = correct_xyz_to_dwell_changed_average(ux, uy, uz, ut, npoints, dwell)
    ptgrid = get_ptgrid(ut)
    ampl, phase = xy_to_ampl_phase(ux, uy)

    tdigit = [ptgrid[1]]
    for dig in 2:size(ptgrid)[1]
        append!(tdigit, ptgrid[dig]-ptgrid[dig-1])
    end

    dtgrid = dtgrid .- dwell
    dtgrid[1] = 0.0
    dtgrid[end] = dtgrid[end] + 0.5*dwell
    append!(dtgrid, dwell*npoints)

    return ampl, phase, ptgrid, tdigit, digdwell, dtgrid
end

"take the average of the hamiltonian of all points until each dwell time\n
output: sumhamdwellx, sumhamdwelly, sumhamdwellz"
function sum_ham_over_dw(hamdwellx, hamdwelly, hamdwellz, nB1, noffs)

    sumhamdwellx=zeros(nB1,noffs,size(hamdwellx)[3])
    sumhamdwelly=zeros(nB1,noffs,size(hamdwellx)[3])
    sumhamdwellz=zeros(nB1,noffs,size(hamdwellx)[3])

    for iB1 = 1:nB1
        for ioffs = 1: noffs
            for dw = 1:size(hamdwellx)[3]
                sumhamdwellx[iB1,ioffs,dw] = sum(hamdwellx[iB1,ioffs,1:dw])/dw
                sumhamdwelly[iB1,ioffs,dw] = sum(hamdwelly[iB1,ioffs,1:dw])/dw
                sumhamdwellz[iB1,ioffs,dw] = sum(hamdwellz[iB1,ioffs,1:dw])/dw
            end
        end
    end

    return sumhamdwellx, sumhamdwelly, sumhamdwellz
end

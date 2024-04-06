using Test


"### simulate_ fid_grad(cstrength, Liouville_space, ux, uy, uz, ut, gamma, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, alignment, gradspins, gradstrength, ngrads, dgrad, gradphase, tubesize)
starting with the pulse parameters of an adiabatic pulse, the spectra, fid and projection are simulated"
function simulate_fid_grad(cstrength, Liouville_space, acqux, acquy, acqgrad, acquz, acqut, gamma, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, ngrads, tubesize)

    acqgraduz = generate_gradient(acqgrad, ngrads, tubesize, gamma)

    facqux, facquy, facquz, facqut, digdwell = get_final_pulse(acqux, acquy, acquz, acqut, npoints, dwell)              # f for final
    facqgraduz = Vector{Float64}[]
    for igrad = 1:ngrads
        fgradux, fgraduy, fgraduz_temp, fgradut, graddigdwell = get_final_pulse(acqux, acquy, tgraduz[igrad], acqut, npoints, dwell)
        push!(facqgraduz, fgraduz_temp)
    end

    println("starting simulation")
    if Liouville_space != "not reduced"
        println("Only NOT reduced space is implemented for gradients")
    end

    fid, allfids = get_fid_Liouville_space_grad(cstrength, ngrads, gradspins, J, rhoinit, digdwell, nB1, noffs, npoints, T2_time, dwell, offs, facqux, facquy, facquz, facqgraduz, facqut)
    time = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done

    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)
    allspectr = zeros(Float64, nB1, noffs, ngrads, npoints)
    for igrad = 1:ngrads
        allspectr[:,:,igrad,:] = fft_fid(allfids[:,:,igrad,:], nB1, noffs, npoints)
    end

    return spectr, fid, allspectr, allfids, projection, time

end



"### get_ fid_ Liouville_ space_grad(cstrength, ngrads, gradspins, J, rhoinit, digdwell, nB1, noffs, npoints, T2_time, dwell, offs, ux, uy, uz, ut, graduz)
get fid calculated in spin operator formalism"
function get_fid_Liouville_space_grad(cstrength, ngrads, gradspins, J, rhoinit, digdwell, nB1, noffs, npoints, T2_time, dwell, offs, ux, uy, uz, graduz, ut)

    ffid = zeros(ComplexF64, nB1, noffs, npoints)
    allfids = zeros(ComplexF64, nB1, noffs, ngrads, npoints)

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # Kopplungsentwicklung zwischen Spin 1 und 2
    elseif cstrength == "strong"
        hj12 = 2*pi*J*(ix[:,:,1]*ix[:,:,2]+iy[:,:,1]*iy[:,:,2]+iz[:,:,1]*iz[:,:,2])
    else
        println("please enter \"strong\" or \"weak\" as coupling strength")
    end

    ampl, phase = xy_to_ampl_phase(ux, uy)

    for iB1 = 1:nB1
        for ioffs = 1:noffs
            fid = zeros(ComplexF64, ngrads, npoints)

            for igrad = 1:ngrads

                # simulating the FID
                rho = rhoinit
                ipulse = 1                  # points of pulse

                for ipoints = 1:npoints     # points of measurement
                    
                    # calculate fid at each acquisition point
                    fid[igrad, ipoints] = tr(rho*(im[:,:,1]))*exp(-(ipoints-1)*dwell/T2_time)

                    # TODO macht noch keinen Sinn, weil Puls ja so oft wiederholt wird bis dwell time vorbei. Bei Gradient aber nicht so sinnvoll?

                    for idig = 1:digdwell[ipoints]        # calculate rho for every digit between two acquisition points
                        phasesin, phasecos = sincosd(phase[ipulse])
                        hpulse = 2*pi*ampl[ipulse]*(phasecos*ix[:,:,2]+phasesin*iy[:,:,2])
                        hcs2 = -2*pi*(offs[ioffs]-uz[ipulse])*iz[:,:,2]    # chemical shift spin 2 including offset
                        hoffs = zeros(size(iz[:,:,1]))
                        for ispin in gradspins
                            hoffs = hoffs .- 2*pi*graduz[igrad][ipulse]*iz[:,:,ispin]
                        end
                        uevo = exp(-i*(hj12+hoffs+hcs2+hpulse)*ut[ipulse])
                        rho = uevo*rho*uevo'
                        ipulse += 1
                    end
                end
            end

            ffid[iB1, ioffs, :] = [sum(fid[:, ipoints])/ngrads for ipoints = 1:npoints]
            allfids[iB1, ioffs, :, :] = fid
        end
    end
    return ffid, allfids
end



"### generate_gradient(grad, ngrads, dgrad, tubesize, gamma)
generates the offset introduced by the gradient \n
ngrads: number of slices in which tube of tubesize is cut into"
function generate_gradient(grad, ngrads, tubesize, gamma)
    graduz = zeros(Float64, ngrads, size(grad))
    for igrad = 1:ngrads
        for idig = 1:size(grad)[1]
            graduz[igrad, idig] = gamma*grad[idig]*(-tubesize/2+(igrad-1)*tubesize/(ngrads-1))
        end
    end
    return graduz
end




"### generate_ gradient_old(alignment, gradstrength, ngrads, dgrad, gradphase, tubesize, ux_pulse, uy_pulse, uz_pulse, ut_pulse, gamma)
generates the offset introduced by the gradient \n
it starts and ends at the same time other pulses are too \n
alignement: pulse can be aligned centered, right or left \n
ngrads: number of slices in which tube of tubesize is cut
dgrad: duration of gradient
gradphase: phase of gradient"
function generate_gradient_old(alignment, gradstrength, ngrads, dgrad, gradphase, tubesize, ux_pulse, uy_pulse, uz_pulse, ut_pulse, gamma)
    graduz = Vector{Float64}[]
    gradut = Float64[]
    for igrad = 1:ngrads
        uz_temp = Float64[]
        ut_temp = Float64[]
        append!(uz_temp, gradphase*gamma*gradstrength*(-tubesize/2+(igrad-1)*tubesize/(ngrads-1)))            # offset of gradient only applied to spin 1 which is detected by FID
        append!(ut_temp, dgrad)

        # alignment
        if sum(ut_pulse) >= dgrad
            if alignment == "centered"
                pushfirst!(uz_temp, 0.0)
                pushfirst!(ut_temp, (sum(ut_pulse)-sum(ut_temp))/2)
                append!(uz_temp, 0.0)
                append!(ut_temp, ut_temp[1])
            elseif alignment == "left"
                append!(uz_temp, 0.0)
                append!(ut_temp, sum(ut_pulse)-sum(ut_temp))
            elseif alignment == "right"
                pushfirst!(uz_temp, 0.0)
                pushfirst!(ut_temp, sum(ut_pulse)-sum(ut_temp))
            else
                println("please enter for alignment: \"centered\", \"left\" or \"right\"")
            end
        end
        push!(graduz, uz_temp)
        gradut = ut_temp
    end

    if sum(ut_pulse) < dgrad
        if alignment == "centered"
            pushfirst!(ux_pulse, 0.0)
            pushfirst!(uy_pulse, 0.0)
            pushfirst!(uz_pulse, 0.0)
            pushfirst!(ut_pulse, (sum(gradut)-sum(ut_pulse))/2)
            append!(ux_pulse, 0.0)
            append!(uy_pulse, 0.0)
            append!(uz_pulse, 0.0)
            append!(ut_pulse, ut_pulse[1])
            @show ut_pulse
        elseif alignment == "left"
            append!(ux_pulse, 0.0)
            append!(uy_pulse, 0.0)
            append!(uz_pulse, 0.0)
            append!(ut_pulse, sum(gradut)-sum(ut_pulse))
        elseif alignment == "right"
            pushfirst!(ux_pulse, 0.0)
            pushfirst!(uy_pulse, 0.0)
            pushfirst!(uz_pulse, 0.0)
            pushfirst!(ut_pulse, sum(gradut)-sum(ut_pulse))
        else
            println("please enter for alignment: \"centered\", \"left\" or \"right\"")
        end
    end


    return ux_pulse, uy_pulse, uz_pulse, ut_pulse, graduz, gradut
end

# ux_pulse, uy_pulse, uz_pulse, ut_pulse, graduz, gradut = generate_gradient("centered", 0.5, 5, 0.003, 1, 0.03, [0.1,0.2,0.3], [0.1,0.2,0.3], [0.1,0.2,0.3], [0.005,0.005,0.005], 42000000)



"### correct_ timing_ two_pulses(ux1, uy1, uz1, ut1, ux2, uy2, uz2, ut2)
make timevector of the two pulses the same"
function correct_timing_two_pulses(ux, uy, uz, ut)
    fux1 = Float64[]
    fuy1 = Float64[]
    fuz1 = Float64[]
    fut1 = Float64[]
    fux2 = Float64[]
    fuy2 = Float64[]
    fuz2 = Float64[]
    fut2 = Float64[]
    ptgrid1 = get_ptgrid(ut[1])
    ptgrid2 = get_ptgrid(ut[2])

    ct = 0.0
    dig1 = 1
    dig2 = 1
    while dig1 <= size(ux[1])[1] && dig2 <= size(ux[2])[1]
        append!(fux1, ux[1][dig1])
        append!(fuy1, uy[1][dig1])
        append!(fuz1, uz[1][dig1])
        append!(fux2, ux[2][dig2])
        append!(fuy2, uy[2][dig2])
        append!(fuz2, uz[2][dig2])
        if ptgrid1[dig1] < ptgrid2[dig2]
            append!(fut1, ptgrid1[dig1]-ct)
            append!(fut2, ptgrid1[dig1]-ct)
            ct = ptgrid1[dig1]
            dig1 += 1
            # println(" next is 1")
        elseif ptgrid1[dig1] > ptgrid2[dig2]
            append!(fut1, ptgrid2[dig2]-ct)
            append!(fut2, ptgrid2[dig2]-ct)
            ct = ptgrid2[dig2]
            dig2 += 1
            # println(" next is 2")
        else 
            append!(fut1, ptgrid1[dig1]-ct)
            append!(fut2, ptgrid2[dig2]-ct)
            ct = ptgrid1[dig1]
            dig1 += 1
            dig2 += 1
            # println("both are next")
        end
        # @show ct, dig1, dig2, size(ux[1]), size(ux[2])
    end

    fux = [fux1 fux2]'
    fuy = [fuy1 fuy2]'
    fuz = [fuz1 fuz2]'
    fut = [fut1 fut2]'


    return fux, fuy, fuz, fut
end


"### correct_ timing_ two_pulses(ux, uy, uz, ut)
make timevector of the two pulses with bilevel the same"
function correct_timing_two_pulses_bilevel(ux, uy, uz, ut)

    fux, fuy, fuz, fut = Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}()
    for ispin = 1:size(ux)[1]
        push!(fux, [[]]); push!(fuy, [[]]); push!(fuz, [[]]); push!(fut, [[]])
    end

    for ilevel = 2:size(ux[1])[1]
        for ispin = 1:size(ux)[1]
            push!(fux[ispin], []); push!(fuy[ispin], []); push!(fuz[ispin], []); push!(fut[ispin], [])
        end
    end

    for ilevel = 1:size(ux[1])[1]
        uxtemp1, uytemp1, uztemp1, uttemp1, uxtemp2, uytemp2, uztemp2, uttemp2 = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
        ptgrid1 = get_ptgrid(ut[1][ilevel])
        ptgrid2 = get_ptgrid(ut[2][ilevel])

        ct = 0.0
        dig1 = 1
        dig2 = 1
        while dig1 <= size(ux[1][ilevel])[1] && dig2 <= size(ux[2][ilevel])[1]
            append!(uxtemp1, ux[1][ilevel][dig1])
            append!(uytemp1, uy[1][ilevel][dig1])
            append!(uztemp1, uz[1][ilevel][dig1])
            append!(uxtemp2, ux[2][ilevel][dig2])
            append!(uytemp2, uy[2][ilevel][dig2])
            append!(uztemp2, uz[2][ilevel][dig2])
            if ptgrid1[dig1] < ptgrid2[dig2]
                append!(uttemp1, ptgrid1[dig1]-ct)
                append!(uttemp2, ptgrid1[dig1]-ct)
                ct = ptgrid1[dig1]
                dig1 += 1
                # println(" next is 1")
            elseif ptgrid1[dig1] > ptgrid2[dig2]
                append!(uttemp1, ptgrid2[dig2]-ct)
                append!(uttemp2, ptgrid2[dig2]-ct)
                ct = ptgrid2[dig2]
                dig2 += 1
                # println(" next is 2")
            else 
                append!(uttemp1, ptgrid1[dig1]-ct)
                append!(uttemp2, ptgrid2[dig2]-ct)
                ct = ptgrid1[dig1]
                dig1 += 1
                dig2 += 1
                # println("both are next")
            end
        end

        fux[1][ilevel] = uxtemp1
        fuy[1][ilevel] = uytemp1
        fuz[1][ilevel] = uztemp1
        fut[1][ilevel] = uttemp1
        fux[2][ilevel] = uxtemp2
        fuy[2][ilevel] = uytemp2
        fuz[2][ilevel] = uztemp2
        fut[2][ilevel] = uttemp2
        
    end


    return fux, fuy, fuz, fut
end

function test_correct_timing_two_pulses_bilevel()
    ux = [[[1,1,2],[2,2,3]], [[3,4,3,4,3,4],[4,5,4,5,4,5]]]
    uy = [[[0,0,0],[0,0,0]], [[0,0,0,0,0,0],[0,0,0,0,0,0]]]
    uz = [[[0,0,0],[0,0,0]], [[0,0,0,0,0,0],[0,0,0,0,0,0]]]
    ut = [[[4,6,5],[4,6,5]], [[3,2,3,2,3,2],[3,2,3,2,3,2]]]
    
    fux, fuy, fuz, fut = correct_timing_two_pulses_bilevel(ux, uy, uz, ut)
    @show fux
    @show fut
    @test ([[[1,1,1,1,1,2,2],[2,2,2,2,2,3,3]], [[3,4,4,3,4,3,4],[4,5,5,4,5,4,5]]], [[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]], [[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]], [[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]], [[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]], [[[3,1,1,3,2,3,2],[3,1,1,3,2,3,2]], [[3,1,1,3,2,3,2],[3,1,1,3,2,3,2]]]) == correct_timing_two_pulses_bilevel(ux, uy, uz, ut)
end

# test_correct_timing_two_pulses_bilevel()

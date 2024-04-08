
"get fid calculated in spin operator formalism\n
output: fid"
function get_fid_pulsesequence(cstrength, rhoinit, sequx, sequy, sequz, sequt, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fid = zeros(ComplexF64, nB1, noffs, npoints)

    ampl, phase = zeros(size(sequx)), zeros(size(sequx))
    for ispin = 1:size(sequx)[1]
        ampl[ispin,:], phase[ispin,:] = xy_to_ampl_phase(sequx[ispin,:], sequy[ispin,:])
    end

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # Kopplungsentwicklung zwischen Spin 1 und 2
    elseif cstrength == "strong"
        hj12 = 2*pi*J*(ix[:,:,1]*ix[:,:,2]+iy[:,:,1]*iy[:,:,2]+iz[:,:,1]*iz[:,:,2])
    else
        println("please enter \"strong\" or \"weak\" as coupling strength")
    end


    for iB1 = 1:nB1
        for ioffs = 1:noffs

            # simulating the FID
            rho = rhoinit

            hcs = zeros(size(iz[:,:,1]))
            for ispin in offsetspins
                hcs = hcs - 2*pi*offs[ioffs]*iz[:,:,ispin]
            end
            # hcs12 = -2*pi*offs[ioffs]*(iz[:,:,1]+iz[:,:,2])    # chemical shift spin 1 and 2 including offset

            for idig = 1:size(ampl)[2]

                hpulse = zeros(size(iz[:,:,1]))
                for ispin in pulsespins
                    phasesin, phasecos = sincosd(phase[ispin,idig])
                    hpulse = hpulse + 2*pi*fB1[iB1]*ampl[ispin,idig]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                    hpulse = hpulse + 2*pi*sequz[ispin,idig]*iz[:,:,ispin]
                end

                uevo = exp(-i*(hj12+hcs+hpulse)*sequt[pulsespins[1],idig])
                rho = uevo*rho*uevo'
            end

            for ipoints = 1:npoints     # points of measurement
                
                # calculate fid at each acquisition point
                fid[iB1, ioffs, ipoints] = tr(rho*(im[:,:,fidspin]))*exp(-(ipoints-1)*dwell/T2_time)
                uevo = exp(-i*(hcs+hj12)*dwell)
                rho = uevo*rho*uevo'
            end
        end
    end

    return fid
end


"spectra, fid and projection are simulated\n
cstrength: coupling strength between the two spins (weak or strong)\n
Liouville space: can be \"reduced\" or \"not reduced\", according to Journal of Magnetic Resonance 201 (2009) 7–17 \n
output: spectr, fid, projection, time"
function simulate_fid_pulsesequence(cstrength, Liouville_space, rhoinit, sequx, sequy, sequz, sequt, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fux, fuy, fuz, fut = get_final_pulse_seq(sequx, sequy, sequz, sequt)
 
    if Liouville_space == "reduced"
        fid = get_fid_pulsesequence_reduced_Liouville_space(rhoinit, digdwell, fux, fuy, fuz, fut, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
        if cstrength == "strong"
            println("only weak coupling impemented for reduced Liouville space. Weak coupling was used instead")
        end
    elseif Liouville_space == "not reduced"
        fid = get_fid_pulsesequence(cstrength, rhoinit, fux, fuy, fuz, fut, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
    else
        println("algorithm not implemented. Please choose \"reduced\" or \"not reduced\"")
    end

    time = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done
    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)
    return spectr, fid, projection, time
end


"get parameters of final pulse for pulsesequence before acquiring the FID\n
if there are several pulsesequences on different nuclei, they are aligned so they end both right before the FID is aquired\n
output: sequx, sequy, sequz, sequt"
function get_final_pulse_seq(sequx, sequy, sequz, sequt)

    fux, fuy, fuz, fut = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()

    # align all pulses to right
    maxpulselength = maximum([sum(sequt[ispin2]) for ispin2 = 1:size(sequt)[1]])
    for ispin = 1:size(sequt)[1]
        if sequt[ispin] != [] && sum(sequt[ispin]) < maxpulselength
            pushfirst!(sequt[ispin], maxpulselength-sum(sequt[ispin]))
            pushfirst!(sequx[ispin], 0.0)
            pushfirst!(sequy[ispin], 0.0)
            pushfirst!(sequz[ispin], 0.0)
        end
    end

    fux, fuy, fuz, fut = correct_timing_two_pulses(sequx, sequy, sequz, sequt)

    return fux, fuy, fuz, fut
end


"starting with the pulse parameters of an adiabatic pulse, the spectra, fid and projection are simulated as BILEVEL\n
output: spectr, fid, projection, time"
function simulate_fid_pulsesequence_bilevel(cstrength, rhoinit, sequxm, sequym, sequzm, sequtm, nbilevel, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fid = zeros(nB1, noffs, npoints)
    sequxm, sequym, sequzm, sequtm = get_final_pulse_seq_bilevel(nbilevel, sequxm, sequym, sequzm, sequtm)

    for ifid = 1:size(sequxm[1])[1]
        println("starting simulation of FID ", ifid)
        fuxm = [sequxm[1][ifid] sequxm[2][ifid]]'            # only for 2 spin system
        fuym = [sequym[1][ifid] sequym[2][ifid]]'
        fuzm = [sequzm[1][ifid] sequzm[2][ifid]]'
        futm = [sequtm[1][ifid] sequtm[2][ifid]]'

        single_fid = get_fid_pulsesequence(cstrength, rhoinit, fuxm, fuym, fuzm, futm, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
        fid = single_fid/size(sequxm[1])[1] + fid
    end

    time = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done
    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)

    return spectr, fid, projection, time
end


"get parameters of final pulse, which is repeated until the end of the fid, and which is cut into two pieces at each dwell time\n
output: sequx, sequy, sequz, sequt"
function get_final_pulse_seq_bilevel(nbilevel::Int64, sequx, sequy, sequz, sequt)

    uxtemp, uytemp, uztemp, uttemp = Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}()
    # fux, fuy, fuz, fut = Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}()
    nspins = Int32(size(sequx)[1])
    maxpulselength = maximum([maximum([sum(sequt[ispin2][ilevel], init=0.0) for ilevel = 1:size(sequt[ispin2])[1]]) for ispin2 = 1:size(sequt)[1]])
    
    for ispin = 1:nspins
        push!(uxtemp, [[]]), push!(uytemp, [[]]), push!(uztemp, [[]]), push!(uttemp, [[]])
    end

    for ispin = 1:nspins
        if isempty(sequx[ispin][1])
            for ilevel = 1:nbilevel
                if ilevel == 1
                    sequx[ispin] = [[0]]; sequy[ispin] = [[0]]; sequz[ispin] = [[0]]; sequt[ispin] = [[maxpulselength]]
                else
                    push!(sequx[ispin], [0])
                    push!(sequy[ispin], [0])
                    push!(sequz[ispin], [0])
                    push!(sequt[ispin], [maxpulselength])
                end
            end
        else
            if size(sequx[ispin])[1] < nbilevel
                factor = ceil(Int32, nbilevel/size(sequx[ispin])[1])
                for ifac = 1:factor-1
                    push!(sequx[ispin], sequx[ispin][1]); push!(sequy[ispin], sequy[ispin][1]); push!(sequz[ispin], sequz[ispin][1]); push!(sequt[ispin], sequt[ispin][1])
                end
            end
            for ilevel = 1:nbilevel
                if sum(sequt[ispin][ilevel]) < maxpulselength
                    pushfirst!(sequt[ispin][ilevel], maxpulselength-sum(sequt[ispin][ilevel]))
                    pushfirst!(sequx[ispin][ilevel], 0.0)
                    pushfirst!(sequy[ispin][ilevel], 0.0)
                    pushfirst!(sequz[ispin][ilevel], 0.0)
                end
            end
        end
    end

    fux, fuy, fuz, fut = correct_timing_two_pulses_bilevel(sequx, sequy, sequz, sequt)

    return fux, fuy, fuz, fut
end


"simulate pulse sequence before fid is acquired and then perform heteronuclear decoupling during the acquisition\n
sequx, sequy, sequz, sequt: controls of pulsesequence before acquisition starts\n
acqux, acquy, acquz, acqut: controls of pulse during acquisition, used for heteronuclear decoupling\n
output: spectr, fid, projection, time"
function simulate_fid_pulsesequence_heteroD(cstrength, Liouville_space, sequx, sequy, sequz, sequt, acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fid = zeros(nB1, noffs, npoints)
    facqux, facquy, facquz, facqut, acqdigdwell = get_final_pulse_acq(acqux, acquy, acquz, acqut, npoints, dwell)
    fux, fuy, fuz, fut = get_final_pulse_seq(sequx, sequy, sequz, sequt)
 
 
    if Liouville_space == "reduced"
        fid = get_fid_pulsesequence__heteroD_reduced_Liouville_space(cstrength, rhoinit, acqdigdwell, fux, fuy, fuz, fut, facqux, facquy, facquz, facqut, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
        if cstrength == "strong"
            println("only weak coupling impemented for reduced Liouville space. Weak coupling was used instead")
        end
    elseif Liouville_space == "not reduced"
        fid = get_fid_pulsesequence_heteroD(cstrength, rhoinit, acqdigdwell, fux, fuy, fuz, fut, facqux, facquy, facquz, facqut, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
    else
        println("algorithm not implemented. Please choose \"reduced\" or \"not reduced\"")
    end

    time = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done
    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)

    return spectr, fid, projection, time
end


"simulate pulse sequence before fid is acquired and then perform bilevel heteronuclear decoupling during the acquisition\n
sequx, sequy, sequz, sequt: controls of pulsesequence before acquisition starts\n
acqux, acquy, acquz, acqut: controls of pulse during acquisition, used for heteronuclear decoupling\n
output: spectr, fid, projection, time"
function simulate_fid_pulsesequence_heteroD_bilevel(cstrength, sequxm, sequym, sequzm, sequtm, acquxm, acquym, acquzm, acqutm, nbilevel::Int64, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fid = zeros(nB1, noffs, npoints)
    acquxm, acquym, acquzm, acqutm, acqdigdwell = get_final_pulse_acq_bilevel(nbilevel, acquxm, acquym, acquzm, acqutm, npoints, dwell)
    sequxm, sequym, sequzm, sequtm = get_final_pulse_seq_bilevel(nbilevel, sequxm, sequym, sequzm, sequtm)

    
    for ifid = 1:size(sequxm[1])[1]     # TODO: same number of bilevel
        println("starting simulation of FID ", ifid)
        fuxm = [sequxm[1][ifid] sequxm[2][ifid]]'            # only for 2 spin system
        fuym = [sequym[1][ifid] sequym[2][ifid]]'
        fuzm = [sequzm[1][ifid] sequzm[2][ifid]]'
        futm = [sequtm[1][ifid] sequtm[2][ifid]]'
        facquxm = [acquxm[1][ifid] acquxm[2][ifid]]'            # only for 2 spin system
        facquym = [acquym[1][ifid] acquym[2][ifid]]'
        facquzm = [acquzm[1][ifid] acquzm[2][ifid]]'
        facqutm = [acqutm[1][ifid] acqutm[2][ifid]]'

        single_fid = get_fid_pulsesequence_heteroD(cstrength, rhoinit, acqdigdwell[ifid], fuxm, fuym, fuzm, futm, facquxm, facquym, facquzm, facqutm, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
        fid = single_fid/size(sequxm[1])[1] + fid
    end

    time = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done
    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)

    return spectr, fid, projection, time
end


"get fid calculated in spin operator formalism with a pulse sequence applied before the fid and heteronuclear coupling during the acquisition\n
output: fid"
function get_fid_pulsesequence_heteroD(cstrength, rhoinit, acqdigdwell, sequx, sequy, sequz, sequt, acqux, acquy, acquz, acqut, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)

    fid = zeros(ComplexF64, nB1, noffs, npoints)
    acqampl, acqphase = zeros(size(acqux)), zeros(size(acqux))
    ampl, phase = zeros(size(sequx)), zeros(size(sequx))
    for ispin = 1:size(acqux)[1]
        acqampl[ispin,:], acqphase[ispin,:] = xy_to_ampl_phase(acqux[ispin,:], acquy[ispin,:])
        ampl[ispin,:], phase[ispin,:] = xy_to_ampl_phase(sequx[ispin,:], sequy[ispin,:])
    end

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # Coupling evolution of spin1 and 2
    elseif cstrength == "strong"
        hj12 = 2*pi*J*(ix[:,:,1]*ix[:,:,2]+iy[:,:,1]*iy[:,:,2]+iz[:,:,1]*iz[:,:,2])
    else
        println("please enter \"strong\" or \"weak\" as coupling strength")
    end

    for iB1 = 1:nB1
        for ioffs = 1:noffs

            # simulating the FID
            rho = rhoinit
            ipulse = 1                  # points of pulse

            hcs = zeros(size(iz[:,:,1]))
            for ispin in offsetspins
                hcs = hcs - 2*pi*offs[ioffs]*iz[:,:,ispin]
            end
            # hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset          # TODO which chemical shift?
            # hcs12 = -2*pi*offs[ioffs]*(iz[:,:,1]+iz[:,:,2])    # chemical shift spin 1 and 2 including offset

            for idig = 1:size(ampl)[2]

                hpulse = zeros(size(iz[:,:,1]))
                for ispin in pulsespins
                    phasesin, phasecos = sincosd(phase[ispin,idig])
                    hpulse = hpulse + 2*pi*fB1[iB1]*ampl[ispin,idig]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                    hpulse = hpulse + 2*pi*sequz[ispin,idig]*iz[:,:,ispin]
                end
    
                uevo = exp(-i*(hj12+hcs+hpulse)*sequt[pulsespins[1],idig])
                rho = uevo*rho*uevo'
            end

            for ipoints = 1:npoints     # points of measurement
                
                # calculate fid at each acquisition point
                fid[iB1, ioffs, ipoints] = tr(rho*(im[:,:,fidspin]))*exp(-(ipoints-1)*dwell/T2_time)

                for idig = 1:acqdigdwell[ipoints]        # calculate rho for every digit between two acquisition points
                    hpulse = zeros(size(iz[:,:,1]))
                    for ispin in pulsespins
                        acqphasesin, acqphasecos = sincosd(acqphase[ispin,ipulse])
                        hpulse = hpulse + 2*pi*fB1[iB1]*acqampl[ispin,ipulse]*(acqphasecos*ix[:,:,ispin]+acqphasesin*iy[:,:,ispin])
                        hpulse = hpulse + 2*pi*acquz[ispin,ipulse]*iz[:,:,ispin]
                    end
                    uevo = exp(-i*(hj12+hcs+hpulse)*acqut[pulsespins[1],ipulse])
                    rho = uevo*rho*uevo'
                    ipulse += 1
                end
            end
        end
    end
    return fid
end


# -------------- FUNCTIONS FOR TESTING ---------------


# function test_get_final_pulse_seq_bilevel()
#     uxtest = [[Float64[]], [[1,1,2], [2,2,3]]]
#     uytest = [[Float64[]], [[0,0,0], [0,0,0]]]
#     uztest = [[Float64[]], [[0,0,0], [0,0,0]]]
#     uttest = [[Float64[]], [[5,5,5], [5,5,5]]]
#     @test ([[[0,0,0],[0,0,0]],[[1,1,2],[2,2,3]]], [[[0,0,0],[0,0,0]],[[0,0,0],[0,0,0]]], [[[0,0,0],[0,0,0]],[[0,0,0],[0,0,0]]], [[[5,5,5],[5,5,5]],[[5,5,5],[5,5,5]]]) == get_final_pulse_seq_bilevel(2, uxtest, uytest, uztest, uttest)

#     uxtest = [[[1,1,2], [2,2,3]], [[4,3,3]]]
#     uytest = [[[0,0,0], [0,0,0]], [[0,0,0]]]
#     uztest = [[[0,0,0], [0,0,0]], [[0,0,0]]]
#     uttest = [[[5,5,5], [5,5,5]], [[1,2,3]]]
#     @test ([[[1,1,1,2,2],[2,2,2,3,3]],[[0,0,4,3,3],[0,0,4,3,3]]], [[[0,0,0,0,0],[0,0,0,0,0]],[[0,0,0,0,0],[0,0,0,0,0]]], [[[0,0,0,0,0],[0,0,0,0,0]],[[0,0,0,0,0],[0,0,0,0,0]]], [[[5,4,1,2,3],[5,4,1,2,3]],[[5,4,1,2,3],[5,4,1,2,3]]]) == get_final_pulse_seq_bilevel(2, uxtest, uytest, uztest, uttest)
# end

# test_get_final_pulse_seq_bilevel()
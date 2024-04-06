


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

    @show sequxm
    @show sequtm
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
            # println("sequx of spin ", ispin, " is empty")
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
            # println("sequx of spin ", ispin, " is not empty")
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
            @show iB1, ioffs

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




# ------------------- NOT WORKING FUNCTIONS SO FAR -------------------------------

"starting with the pulse parameters of an adiabatic pulse, the spectra, fid and projection are simulated"
function simulate_fid_sequence_homoD(cstrength, Liouville_space, ifplot, sequx, sequy, sequz, sequt, homux, homuy, homuz, homut, rhoinit, olddchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

    fux, fuy, fuz, fut = get_final_pulse_seq(sequx, sequy, sequz, sequt)

    # dchunk length is defined to be multiple of dwelltime
    dchunk = round(olddchunk/(dwell))*dwell
    nchunks = ceil(Int32,(dwell*npoints-dchunk)/(2*dchunk)) +1                  # amount of chunks in total
    println(nchunks, " chunks with the chunklength of ", dchunk, " instead of ", olddchunk, " are used")

    # tmidchunk: time start of each chunk
    # tspulse: time start of each pulse
    tmidchunk = zeros(nchunks)                                        # is defined as the middle of the chunk, but for the first chunk only the second half starting at 0 is used
    tmidchunk[1] = 0.0
    tspulse = zeros(nchunks)
    tspulse[1] = -maximum([sum(homut[ispin]) for ispin = 1:size(homut)[1]])                                            # That is why only homonuclear Decoupling only works for one decoupling sequence at once
    for ichunk = 2:nchunks
        tmidchunk[ichunk] = tmidchunk[ichunk-1] + 2*dchunk
        tspulse[ichunk] = (2*tspulse[1]+tmidchunk[ichunk])/2
    end

    # fid: final FID which is used for FFT
    # sfid: single fid, only saved for easier debugging so far
    fid = zeros(ComplexF64, nB1, noffs, npoints)
    sfid = zeros(ComplexF64, nB1, noffs, nchunks, npoints)            # single fid, for each chunk

    tim = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done

    if Liouville_space != "not reduced"
        println("Only algorithm for \"not reduced\" Liouville space is possible, which is used instead.")
    end

    # rhobeg (rho-beginning): rho for each dwelltime until start of pulse in last FID without any pulse applied
    rhobeg = get_rho_without_pulse(tspulse[end], rhoinit, cstrength, tspulse[1], J, noffs, dwell, offs)

    # for each chunk one FID is simulated and the chunk is added to the final FID
    # chunk starts at tmidchunk-dchunk and ends at tmidchunk+dchunk
    for ifid = 1:nchunks
        digdwell = Float64[]
        fhomux, fhomuy, fhomuz, fhomut, preux, preuy, preuz, preut = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
        for ispin = 1:size(homux)[1]
            push!(fhomux, [])
            push!(fhomuy, [])
            push!(fhomuz, [])
            push!(fhomut, [])
            push!(preux, [])
            push!(preuy, [])
            push!(preuz, [])
            push!(preut, [])
        end
        for ispin = 1:size(homux)[1]
            fhomux[ispin], fhomuy[ispin], fhomuz[ispin], fhomut[ispin], preux[ispin], preuy[ispin], preuz[ispin], preut[ispin] = correct_pulse_to_interferogram(homux[ispin], homuy[ispin], homuz[ispin], homut[ispin], tspulse[ifid], tmidchunk[ifid], dchunk, dwell)      # only for one type of pulse working so far -> no timing needs to be changed
            fhomux[ispin], fhomuy[ispin], fhomuz[ispin], fhomut[ispin], digdwell = correct_xyz_to_dwell(fhomux[ispin], fhomuy[ispin], fhomuz[ispin], fhomut[ispin], npoints, dwell)
        end
        
        println("starting simulation of FID ", ifid)
        # TODO not working yet
        single_fid =  get_fid_sequence_homoD(cstrength, fB1, rhoinit, rhobeg, digdwell, fhomux, fhomuy, fhomuz, fhomut, preux, preuy, preuz, preut, tspulse[ifid], tmidchunk[ifid], dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

        # only keep necessary part of fid (chunk itself), rest is set to zero
        # sfid[:,ifid,:] = single_fid
        if ifid == 1
            fid[:,:,1:Int32(dchunk/dwell)] = single_fid[:,:,1:Int32(dchunk/dwell)]
        else
            if Int32(dchunk/dwell)*(2*ifid-1) <= size(fid)[2]
                fid[:,:,Int32(dchunk/dwell)*(2*ifid-3)+1:Int32(dchunk/dwell)*(2*ifid-1)] = single_fid[:,:,Int32(dchunk/dwell)*(2*ifid-3)+1:Int32(dchunk/dwell)*(2*ifid-1)]
            else
                fid[:,:,Int32(dchunk/dwell)*(2*ifid-3)+1:end] = single_fid[:,:,Int32(dchunk/dwell)*(2*ifid-3)+1:end]
            end
        end
        sfid[:,:,ifid,:] = fid

    end

    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)

    return spectr, fid, projection, tim, sfid
end



"### get_ fid_ homoD_ Liouville_space(cstrength, fB1, rhoinit, rhobeg, digdwell, sequx, sequy, sequz, sequt, preux, preuy, preuz, preut, tspulse, tmidchunk, dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time)
get fid of the chunk \n
apply pre controls first, start acquiring FID at t=0 when final controls (fut etc.) are used"
function get_fid_sequence_homoD(cstrength, fB1, rhoinit, rhobeg, digdwell, sequx, sequy, sequz, sequt, homux, homuy, homuz, homut, preux, preuy, preuz, preut, tspulse, tmidchunk, dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
    fid = zeros(ComplexF64, nB1, noffs, npoints)

    ampl, phase = zeros(size(sequx)), zeros(size(sequx))
    for ispin = 1:size(sequx)[1]
        ampl[ispin,:], phase[ispin,:] = xy_to_ampl_phase(sequx[ispin,:], sequy[ispin,:])
    end

    homampl, homphase, hompreampl, homprephase = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    for ispin = 1:size(homux)[1]
        push!(homampl, [])
        push!(homphase, [])
        push!(hompreampl, [])
        push!(homprephase, [])
    end

    for ispin = 1:size(homux)[1]
        homampl[ispin], homphase[ispin] = xy_to_ampl_phase(homux[ispin], homuy[ispin])
        hompreampl[ispin], homprephase[ispin] = xy_to_ampl_phase(preux[ispin], preuy[ispin])
    end

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # coupling evolution beetween spin 1 and 2
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

            for idig = 1:size(ampl)[2]              # TODO when to use which amplitude and phase??? --> needs to be corrected

                hpulse = zeros(size(iz[:,:,1]))
                for ispin in pulsespins
                    phasesin, phasecos = sincosd(phase[ispin,idig])
                    hpulse = hpulse + 2*pi*fB1[iB1]*ampl[ispin,idig]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                    hpulse = hpulse + 2*pi*sequz[ispin,idig]*iz[:,:,ispin]
                end

                uevo = exp(-i*(hj12+hcs+hpulse)*sequt[pulsespins[1],idig])
                rho = uevo*rho*uevo'
            end
            # now rho is at end of sequence = at beginning of FID and Homo Decoupling
            # how to handle homo Decoupling pulse starting before FID is taken, how long will pulse sequence be then?
            if tspulse <= 0
                startfid = 0
                startpulse = 0
                global rho = rhoinit
            else
                startfid = floor(Int64,(tmidchunk-dchunk)/dwell)
                startpulse = floor(Int64,tspulse/dwell)
                
                global rho = rhobeg[ioffs, startpulse+1,:,:]
            end

            # simulating the FID
            ipulse = 1                  # points of pulse

            hcs = zeros(size(iz[:,:,1]))
            for ispin in offsetspins
                hcs = hcs - 2*pi*offs[ioffs]*iz[:,:,ispin]
            end

            for iprepulse = 1:size(preux[1])[1]
                hpulse = zeros(size(iz[:,:,1]))
                for ispin in pulsespins
                    phasesin, phasecos = sincosd(prephase[ispin][iprepulse])
                    hpulse = hpulse + 2*pi*fB1[iB1]*preampl[ispin][iprepulse]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                    hpulse = hpulse + 2*pi*preuz[ispin][ipulse]*iz[:,:,ispin]
                end
                # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                # hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset
                uevo = exp(-i*(hj12+hcs+hpulse)*preut[pulsespins[1]][iprepulse])
                global rho = uevo*rho*uevo'
            end

            for ipoints = 1:minimum([ceil(Int64,(tmidchunk+dchunk)/dwell - startpulse), npoints-startpulse])
                
                # calculate fid at each acquisition point
                # fid[iB1, ioffs, startpulse+ipoints]=tr(rho*(im[:,:,1]+im[:,:,2]))*exp(-(startpulse+ipoints-1)*dwell/T2_time)
                fid[iB1, ioffs, startpulse+ipoints]=tr(rho*im[:,:,fidspin])*exp(-(startpulse+ipoints-1)*dwell/T2_time)

                for idig = 1:digdwell[ipoints]       # calculate rho for every digit between two acquisition points
                    hpulse = zeros(size(iz[:,:,1]))
                    for ispin in pulsespins
                        phasesin, phasecos = sincosd(phase[ispin][ipulse])
                        hpulse = hpulse + 2*pi*fB1[iB1]*ampl[ispin][ipulse]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                        hpulse = hpulse + 2*pi*homuz[ispin][ipulse]*iz[:,:,ispin]
                    end
                    # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                    # hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset
                    uevo = exp(-i*(hj12+hcs+hpulse)*homut[pulsespins[1]][ipulse])
                    global rho = uevo*rho*uevo'
                    ipulse += 1
                end
            end
        end
    end
    return fid
end




# -------------- FUNCTIONS FOR TESTING ---------------


function test_get_final_pulse_seq_bilevel()
    uxtest = [[Float64[]], [[1,1,2], [2,2,3]]]
    uytest = [[Float64[]], [[0,0,0], [0,0,0]]]
    uztest = [[Float64[]], [[0,0,0], [0,0,0]]]
    uttest = [[Float64[]], [[5,5,5], [5,5,5]]]
    @test ([[[0,0,0],[0,0,0]],[[1,1,2],[2,2,3]]], [[[0,0,0],[0,0,0]],[[0,0,0],[0,0,0]]], [[[0,0,0],[0,0,0]],[[0,0,0],[0,0,0]]], [[[5,5,5],[5,5,5]],[[5,5,5],[5,5,5]]]) == get_final_pulse_seq_bilevel(2, uxtest, uytest, uztest, uttest)

    uxtest = [[[1,1,2], [2,2,3]], [[4,3,3]]]
    uytest = [[[0,0,0], [0,0,0]], [[0,0,0]]]
    uztest = [[[0,0,0], [0,0,0]], [[0,0,0]]]
    uttest = [[[5,5,5], [5,5,5]], [[1,2,3]]]
    @test ([[[1,1,1,2,2],[2,2,2,3,3]],[[0,0,4,3,3],[0,0,4,3,3]]], [[[0,0,0,0,0],[0,0,0,0,0]],[[0,0,0,0,0],[0,0,0,0,0]]], [[[0,0,0,0,0],[0,0,0,0,0]],[[0,0,0,0,0],[0,0,0,0,0]]], [[[5,4,1,2,3],[5,4,1,2,3]],[[5,4,1,2,3],[5,4,1,2,3]]]) == get_final_pulse_seq_bilevel(2, uxtest, uytest, uztest, uttest)
end

# test_get_final_pulse_seq_bilevel()
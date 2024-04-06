include("fid_2spins_heteronoclear.jl")


"simulation of a spectrum where homonuclear decoupling of the interferogram type is applied\n
starting with the pulse parameters of an adiabatic pulse, the spectra, fid and projection are simulated\n
sfid: single fids of each chunk\n
output: spectr, fid, projection, time, sfid"
function simulate_fid_homoD(cstrength, Liouville_space, ifplot, homux, homuy, homuz, homut, rhoinit, olddchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

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
    # for ifid = 1:2
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
        single_fid =  get_fid_homoD_Liouville_space(cstrength, fB1, rhoinit, rhobeg, digdwell, fhomux, fhomuy, fhomuz, fhomut, preux, preuy, preuz, preut, tspulse[ifid], tmidchunk[ifid], dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

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

        if ifplot == true
            # plot each fid and pulse over time
            testptime = []
            testpulsex = []
            testpulsey = []
            if isempty(preut)
                testptime = get_ptgrid([tspulse[ifid]; fut])
                testpulsex = [0; fux]
                testpulsey = [0; fuy]
            else 
                testptime = get_ptgrid(preut)
                testptime = reverse(.- testptime)
                testpulsex = preux
                append!(testpulsex, fux)
                testpulsey = preuy
                append!(testpulsey, fuy)
                test2ptime = get_ptgrid(fut)
                append!(testptime, test2ptime)
            end

            savingname="Figures/HomoD/Hard180_FID_nr" * string(ifid) * ".png"

            testtime = [i * dwell for i in 0:(npoints-1)]       # times at which acquisition was done
            plot(testptime, [testpulsex testpulsey], xlim = (-0.01, dwell*npoints), label=["ux" "uy"])
            display(plot!(twinx(),testtime, real(single_fid[ceil(nB1, noffs/2),ceil(Int32, noffs/2),:]), color=:green, label="FID", ylim=(-2,2), xlim = (-0.01, dwell*npoints)))
            # savefig(savingname)
        end

    end

    spectr = fft_fid(fid, nB1, noffs, npoints)
    projection = get_projection(spectr, nB1, noffs)

    return spectr, fid, projection, tim, sfid
end


"pulse controls are generated for certain chunk: \n
FID is made up of two parts: \n
either prefut and fut: \n
fut starts at t=0, but controls are zero until pulse actually starts for this specific chunk, and after pulse is finished while chunk-FID is acquired \n
prefut starts at tspulse if it's <0, but everything after t=0 is still included in fut \n
or rhostart (rho calculated without pulse) and fut starting at last dwellpoint before pulse starts\n
output: ux, uy, uz, ut, preux, preuy, preuz, preut"
function correct_pulse_to_interferogram(ux, uy, uz, ut, tspulse, tschunk, dchunk, dwell)

    ptgrid = get_ptgrid(ut)
    preux = Float64[]                                       # pre ux contains all pulse digits which have to be simulated before FID is being measured
    preuy = Float64[]
    preuz = Float64[]
    preut = Float64[]
    fux = Float64[]                                         # fux contains all pulse digits which are simulated during FID
    fuy = Float64[]
    fuz = Float64[]
    fut = Float64[]
    pulse_duration = ptgrid[end]
    # ct = tschunk - pulse_duration                           # current time
    ct = tspulse                           # current time

    if ct < 0
        dig = 1
        while dig <= size(ux)[1] && ct+ut[dig] <= 0         # while pulse is before FID starts
            append!(preux, ux[dig])
            append!(preuy, uy[dig])
            append!(preuz, uz[dig])
            append!(preut, ut[dig])
            ct = ct + ut[dig]
            dig += 1
        end
        if dig <= size(ux)[1]                               
            append!(preux, ux[dig])
            append!(preuy, uy[dig])
            append!(preuz, uz[dig])
            append!(preut, abs(ct))
            ct = 0.0
            append!(fux, ux[dig])
            append!(fuy, uy[dig])
            append!(fuz, uz[dig])
            append!(fut, ut[dig] - preut[end])
            dig += 1
            ct = ct + fut[end]
        end
        while dig <= size(ux)[1]                            # pulse during FID is measured, still before chunk
            append!(fux, ux[dig])
            append!(fuy, uy[dig])
            append!(fuz, uz[dig])
            append!(fut, ut[dig])
            ct = ct + ut[dig]
            dig += 1
        end
        if tspulse+sum(ut) < tschunk                        # if pulse is finished before chunk starts
            append!(fux, 0.0, 0.0)
            append!(fuy, 0.0, 0.0)
            append!(fuz, 0.0, 0.0)
            append!(fut, tschunk-tspulse-sum(ut), dchunk)
        elseif tspulse+sum(ut) == tschunk                   # if pulse ends exactly when chunk starts
            append!(fux, 0.0)
            append!(fuy, 0.0)
            append!(fuz, 0.0)
            append!(fut, dchunk)
        end

    else                                                    # if pulse starts later than simulation of FID
        fux = Float64[0.0]                                  
        fuy = Float64[0.0]
        fuz = Float64[0.0]
        fut = Float64[tspulse-floor(tspulse/dwell)*dwell]       # only this changed in order to use rhostart
        if tspulse+sum(ut) < tschunk                        # if pulse is finished before chunk starts
            append!(fux, ux, 0.0, 0.0)
            append!(fuy, uy, 0.0, 0.0)
            append!(fuz, uz, 0.0, 0.0)
            append!(fut, ut, tschunk-tspulse-sum(ut), dchunk)
        elseif tspulse+sum(ut) == tschunk                   # if pulse ends exactly when chunk starts
            append!(fux, ux, 0.0)
            append!(fuy, uy, 0.0)
            append!(fuz, uz, 0.0)
            append!(fut, ut, dchunk)
        end
    end
    return fux, fuy, fuz, fut, preux, preuy, preuz, preut
end


"get fid of the single chunk \n
apply pre controls first, start acquiring FID at t=0 when final controls (fut etc.) are used\n
output: fid"
function get_fid_homoD_Liouville_space(cstrength, fB1, rhoinit, rhobeg, digdwell, homux, homuy, homuz, homut, preux, preuy, preuz, preut, tspulse, tmidchunk, dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
    fid = zeros(ComplexF64, nB1, noffs, npoints)

    ampl, phase, preampl, prephase = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    for ispin = 1:size(homux)[1]
        push!(ampl, [])
        push!(phase, [])
        push!(preampl, [])
        push!(prephase, [])
    end

    for ispin = 1:size(homux)[1]
        ampl[ispin], phase[ispin] = xy_to_ampl_phase(homux[ispin], homuy[ispin])
        preampl[ispin], prephase[ispin] = xy_to_ampl_phase(preux[ispin], preuy[ispin])
    end

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # coupling evolution beetween spin 1 and 2
    elseif cstrength == "strong"
        hj12 = 2*pi*J*(ix[:,:,1]*ix[:,:,2]+iy[:,:,1]*iy[:,:,2]+iz[:,:,1]*iz[:,:,2])
    else
        println("please enter \"strong\" or \"weak\" as coupling strength")
    end

    for iB1 = 1:nB1
        for ioffs=1:noffs

            if tspulse <= 0
                # startfid = 0
                startpulse = 0
                global rho = rhoinit
            else
                # startfid = floor(Int64,(tmidchunk-dchunk)/dwell)
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
                    hpulse = hpulse + 2*pi*preuz[ispin][iprepulse]*iz[:,:,ispin]
                end
                # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                # hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset
                uevo = exp(-i*(hj12+hcs+hpulse)*preut[pulsespins[1]][iprepulse])
                global rho = uevo*rho*uevo'
            end

            for ipoints = 1:minimum([ceil(Int64,(tmidchunk+dchunk)/dwell - startpulse), npoints-startpulse])
                
                # calculate fid at each acquisition point
                fid[iB1, ioffs, startpulse+ipoints] = tr(rho*im[:,:,fidspin])*exp(-(startpulse+ipoints-1)*dwell/T2_time)

                for idig = 1:digdwell[ipoints]       # calculate rho for every digit between two acquisition points
                    hpulse = zeros(size(iz[:,:,1]))
                    for ispin in pulsespins
                        phasesin, phasecos = sincosd(phase[ispin][ipulse])
                        hpulse = hpulse + 2*pi*fB1[iB1]*ampl[ispin][ipulse]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                        hpulse = hpulse + 2*pi*homuz[ispin][ipulse]*iz[:,:,ispin]
                    end
                    # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                    uevo = exp(-i*(hj12+hcs+hpulse)*homut[pulsespins[1]][ipulse])
                    global rho = uevo*rho*uevo'
                    ipulse += 1
                end
            end
        end
    end
    return fid
end


"simulate rho for each dwellpoint until time tend (time-end) is reached if no pulse is applied\n
output: rho"
function get_rho_without_pulse(tend, rhoinit, cstrength, fidstart, J, noffs, dwell, offs)
    rho = zeros(ComplexF64, noffs, ceil(Int64,tend/dwell)+1, size(rhoinit)[1], size(rhoinit)[2])

    if cstrength == "weak"
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # Kopplungsentwicklung zwischen Spin 1 und 2
    elseif cstrength == "strong"
        hj12 = 2*pi*J*(ix[:,:,1]*ix[:,:,2]+iy[:,:,1]*iy[:,:,2]+iz[:,:,1]*iz[:,:,2])
    else
        println("please enter \"strong\" or \"weak\" as coupling strength")
    end

    

    for ioffs=1:noffs
        # simulating the FID
        hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset
        uevo = exp(-i*(hj12+hcs2)*abs(fidstart))
        rho[ioffs, 1,:,:] = uevo * rhoinit * uevo'

        for ipoints = 2:size(rho)[2]     # points of measurement
            uevo = exp(-i*(hj12+hcs2)*dwell)
            rho[ioffs,ipoints,:,:] = uevo * rho[ioffs,ipoints-1,:,:] * uevo'
        end
    end

    return rho
end


"starting with the pulse parameters of an adiabatic pulse, the spectra, fid and projection are simulated where homonuclear and heteronuclear decoupling is applied\n
sfid: single fids of each chunk\n
output: spectr, fid, projection, time, sfid"
function simulate_fid_heteroD_homoD(cstrength, acqux, acquy, acquz, acqut, homux, homuy, homuz, homut, rhoinit, olddchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

     # dchunk length is defined to be multiple of dwelltime
    dchunk = round(olddchunk/(dwell))*dwell
    nchunks = ceil(Int32,(dwell*npoints-dchunk)/(2*dchunk)) +1                  # amount of chunks in total
    println(nchunks, " chunks with the chunklength of ", dchunk, " instead of ", olddchunk, " are used")
 
    facqux, facquy, facquz, facqut, acqdigdwell = get_final_pulse_acq_during_homoD(acqux, acquy, acquz, acqut, dchunk)
 
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

    # rhobeg (rho-beginning): rho for each dwelltime until start of pulse in last FID without any pulse applied
    rhobeg = get_rho_without_pulse(tspulse[end], rhoinit, cstrength, tspulse[1], J, noffs, dwell, offs)

    # for each chunk one FID is simulated and the chunk is added to the final FID
    # chunk starts at tmidchunk-dchunk and ends at tmidchunk+dchunk
    # for ifid = 1:2
    for ifid = 1:nchunks
        homdigdwell = Float64[]
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
            fhomux[ispin], fhomuy[ispin], fhomuz[ispin], fhomut[ispin], homdigdwell = correct_xyz_to_dwell(fhomux[ispin], fhomuy[ispin], fhomuz[ispin], fhomut[ispin], npoints, dwell)
        end
        println("starting simulation of FID ", ifid)
        single_fid =  get_fid_heteroD_homoD(cstrength, fB1, rhoinit, rhobeg, facqux, facquy, facquz, facqut, acqdigdwell, fhomux, fhomuy, fhomuz, fhomut, homdigdwell, preux, preuy, preuz, preut, tspulse[ifid], tmidchunk[ifid], dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)

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






"get parameters of final pulse, which is repeated until the end of the fid, and which is cut into two pieces at each dwell time\n
output: acqux, acquy, acquz, acqut, acqdigdwell"
function get_final_pulse_acq_during_homoD(acqux, acquy, acquz, acqut, dchunk)

    uxtemp, uytemp, uztemp, uttemp = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    facqux, facquy, facquz, facqut = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    acqdigdwell = []

    for ispin in eachindex(ux)
        if ux[ispin] == []
            push!(uxtemp, [0])
            push!(uytemp, [0])
            push!(uztemp, [0])
            push!(uttemp, [dchunk])
        else
            uxrep, uyrep, uzrep, utrep = repeat_pulse_during_homoD(acqux[ispin], acquy[ispin], acquz[ispin], acqut[ispin], dchunk)
            push!(uxtemp, uxrep)
            push!(uytemp, uyrep)
            push!(uztemp, uzrep)
            push!(uttemp, utrep)
        end
    end

    uxtemp, uytemp, uztemp, uttemp = correct_timing_two_pulses(uxtemp, uytemp, uztemp, uttemp)

    for ispin in eachindex(ux)
        uxdwell, uydwell, uzdwell, utdwell, acqdigdwell = correct_xyz_to_dwell(uxtemp[ispin,:], uytemp[ispin,:], uztemp[ispin,:], uttemp[ispin,:], npoints, dwell)
        push!(facqux, uxdwell)
        push!(facquy, uydwell)
        push!(facquz, uzdwell)
        push!(facqut, utdwell)
    end
    
    facqux = reduce(hcat, facqux)'
    facquy = reduce(hcat, facquy)'
    facquz = reduce(hcat, facquz)'
    facqut = reduce(hcat, facqut)'

    return facqux, facquy, facquz, facqut, acqdigdwell
end



"pulse is repeated until it is longer than the acquisition time\n
output: acqux, facquy, facquz, facqut, acqdigdwell"
function repeat_pulse_during_homoD(acqux, acquy, acquz, acqut, dchunk)
    pulse_duration = sum(acqut)

    if pulse_duration < dchunk      # if pulse is shorter than dwell time
        factor = ceil(Int, dchunk / pulse_duration)
        acqux = repeat(acqux, factor)
        acquy = repeat(acquy, factor)
        acquz = repeat(acquz, factor)
        acqut = repeat(acqut, factor)
        return acqux, acquy, acquz, acqut
    end
end




"get fid of one chunk where homonuclear and heteronuclear decoupling is applied \n
apply pre controls first, start acquiring FID at t=0 when final controls (fut etc.) are used\n
output: fid"
function get_fid_heteroD_homoD(cstrength, fB1, rhoinit, rhobeg, acqux, acquy, acquz, acqut, acqdigdwell, homux, homuy, homuz, homut, homdigdwell, preux, preuy, preuz, preut, tspulse, tmidchunk, dchunk, J, offs, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
    fid = zeros(ComplexF64, nB1, noffs, npoints)

    acqampl, acqphase = zeros(size(acqux)), zeros(size(acqux))
    for ispin = 1:size(acqux)[1]
        acqampl[ispin,:], acqphase[ispin,:] = xy_to_ampl_phase(acqux[ispin,:], acquy[ispin,:])
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

            if tspulse <= 0
                # startfid = 0
                startpulse = 0
                global rho = rhoinit
            else
                # startfid = floor(Int64,(tmidchunk-dchunk)/dwell)
                startpulse = round(Int64,tspulse/dwell)
                
                global rho = rhobeg[ioffs, startpulse+1,:,:]
            end

            # simulating the FID
            ihompulse = 1                  # points of pulse
            iacqpulse = 1

            hcs = zeros(size(iz[:,:,1]))
            for ispin in offsetspins
                hcs = hcs - 2*pi*offs[ioffs]*iz[:,:,ispin]
            end

            for iprepulse = 1:size(preux[1])[1]
                hpulse = zeros(size(iz[:,:,1]))
                for ispin in pulsespins
                    phasesin, phasecos = sincosd(homprephase[ispin][iprepulse])
                    hpulse = hpulse + 2*pi*fB1[iB1]*hompreampl[ispin][iprepulse]*(phasecos*ix[:,:,ispin]+phasesin*iy[:,:,ispin])
                    hpulse = hpulse + 2*pi*preuz[ispin][iprepulse]*iz[:,:,ispin]
                end
                # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                # hcs2 = -2*pi*offs[ioffs]*iz[:,:,2]    # chemical shift spin 2 including offset
                uevo = exp(-i*(hj12+hcs+hpulse)*preut[pulsespins[1]][iprepulse])
                global rho = uevo*rho*uevo'
            end

            for ipoints = 1:maximum([0,minimum([floor(Int64,(tmidchunk-dchunk)/dwell - startpulse), npoints-startpulse])])

                for idig = 1:homdigdwell[ipoints]       # calculate rho for every digit between two acquisition points
                    hpulse = zeros(size(iz[:,:,1]))
                    for ispin in pulsespins
                        homphasesin, homphasecos = sincosd(homphase[ispin][ihompulse])
                        hpulse = hpulse + 2*pi*fB1[iB1]*homampl[ispin][ihompulse]*(homphasecos*ix[:,:,ispin]+homphasesin*iy[:,:,ispin])
                        hpulse = hpulse + 2*pi*homuz[ispin][ihompulse]*iz[:,:,ispin]
                    end
                    # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                    uevo = exp(-i*(hj12+hcs+hpulse)*homut[pulsespins[1]][ihompulse])
                    global rho = uevo*rho*uevo'
                    ihompulse += 1
                end
            end


            for ipoints = maximum([1, minimum([floor(Int64,(tmidchunk-dchunk)/dwell - startpulse), npoints-startpulse])+1]):minimum([ceil(Int64,(tmidchunk+dchunk)/dwell - startpulse), npoints-startpulse])

                fid[iB1, ioffs, startpulse+ipoints] = tr(rho*im[:,:,fidspin])*exp(-(startpulse+ipoints-1)*dwell/T2_time)

                for idig = 1:acqdigdwell[ipoints]       # calculate rho for every digit between two acquisition points
                    hpulse = zeros(size(iz[:,:,1]))
                    for ispin in pulsespins
                        acqphasesin, acqphasecos = sincosd(acqphase[ispin,iacqpulse])
                        hpulse = hpulse + 2*pi*fB1[iB1]*acqampl[ispin,iacqpulse]*(acqphasecos*ix[:,:,ispin]+acqphasesin*iy[:,:,ispin])
                        hpulse = hpulse + 2*pi*acquz[ispin,iacqpulse]*iz[:,:,ispin]
                    end
                    # TODO: Pulse direkt auf anderer Frequenz einstrahlen, dazu Phase verändern. "Richtige" Phase muss in der Mitte sein (bei Bruker 0.5)
                    uevo = exp(-i*(hj12+hcs+hpulse)*acqut[pulsespins[1],iacqpulse])
                    global rho = uevo*rho*uevo'
                    iacqpulse += 1
                end
            end
        end
    end
    return fid
end





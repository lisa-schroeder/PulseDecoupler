using StaticArrays



"variables are initialized as given in paramfile, simulations defined in paramfile will be performed\\
results will be all global variables, dependent on the simulation\n
output: nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm"
function performsimulations(paramfilename)
    nspins, npoints, nB1, noffs, fB1, offs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, nbilevel::Int64, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, thresh, dchunk, eps, startmag, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH, bilevel_flag, iu, ix, iy, iz, ip, im, ia, ib = initialize(paramfilename)
    @time simulate(paramfilename, nspins, npoints, nB1, noffs, fB1, offs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, nbilevel, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, thresh, dchunk, eps, startmag, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH, bilevel_flag, iu, ix, iy, iz, ip, im, ia, ib)
    # B1 = fB1 .* maximum([maximum(rfmaxseq), maximum(rfmaxacq), maximum(rfmaxhom)])
    B1 = fB1 .* maximum(rfmaxacq[2])
    return nB1, noffs, B1, offs, dwell, npoints, allJs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm
end



"### initialize(paramfilename)
all variables are defined and their values written in paramfilename will be assigned and returned\n
output: nspins, npoints, nB1, noffs, fB1, offs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, nbilevel, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, thresh, dchunk, eps, startmag, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH, bilevel_flag, iu, ix, iy, iz, ip, im, ia, ib"
function initialize(paramfilename)
    println("Input file ", paramfilename, " will be used for analyzation setup")
    io = open(paramfilename, "r")
    data = readlines(paramfilename)

    simulate_fid_flag = false
    homo_fid_flag = false
    any_homo_flag = false
    hetero_fid_flag = false
    any_acq_flag = false
    sequence_fid_flag = false
    spins10_flag = false
    magnetization_simulation_flag = false
    avH_flag = false
    waugh_flag = false
    bilevel_flag = false
    sideband_flag = false
    gradient_flag = false
    xyzgt_flag = false
    pulsetype_flag = false
    nbilevel = 1::Int64
    nspins = 0
    # sequxtemp, uytemp, uztemp, seqgradtemp, uttemp, acquxtemp, acquytemp, acquztemp, acqgradtemp, acquttemp, homuxtemp, homuytemp, homuztemp, homgradtemp, homuttemp = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(),Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()

    for count in eachindex(data)

        currentline = data[count]
        if occursin("#", currentline)
            currentline = currentline[1:findfirst("#", data[count])[1]-1]
        end

        if occursin("nspins", currentline)
            ff = findfirst("=", currentline)[1]
            nspins = parse(Int64, currentline[ff+1:end])
            println("nspins:               ", nspins)
        end

    end

    npoints, nB1, noffs, pdigits, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, tpulseseq, tpulseacq, tpulsehom, bwdth, n, k, thresh, dchunk, eps, pulsefilename, pulsetype, supercycle, startmag, offsrange, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad,gradspins, gradstrength, tubesize, gammaH = initializeallparams(nspins)

    for count in eachindex(data)

        currentline = data[count]
        if occursin("#", currentline)
            currentline = currentline[1:findfirst("#", data[count])[1]-1]
        end

        if occursin("pulsetype", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            pulsetype[spinnumber] = lstrip(rstrip(currentline[ff+1:end]))
            pulsetype_flag = true
            println("pulsetype of spin ", spinnumber, ":  ", pulsetype[spinnumber])
        elseif occursin("simulate magnetization", currentline)
            magnetization_simulation_flag = true
        elseif occursin("average H", currentline)
            avH_flag = true
        elseif occursin("Waugh", currentline)
            waugh_flag = true
        elseif occursin("simulate sidebands", currentline)
            sideband_flag = true
        elseif occursin("xyzgt", currentline)
            xyzgt_flag = true
        end
        if occursin("Homo", currentline) && occursin("simulate FID", currentline)
            homo_fid_flag = true
        end
        if occursin("Homo", currentline) || occursin("hom", currentline)
            any_homo_flag = true
        end
        if occursin("Hetero", currentline) && occursin("simulate FID", currentline)
            hetero_fid_flag = true
        end
        if occursin("Hetero", currentline) || occursin("acq", currentline)
            any_acq_flag = true
        end
        if occursin("Sequence", currentline) && occursin("simulate FID", currentline)
            sequence_fid_flag = true
        end
        if occursin("10 spins", currentline)
            spins10_flag = true
        end
        if occursin("gradient", currentline)
            gradient_flag = true
        end
        if occursin("bilevel", currentline)
            bilevel_flag = true
            nbilevel = parse(Float64, split(currentline)[end])
            println("level of bilevel:     ", nbilevel)
        end
        if occursin("simulate FID", currentline)
            simulate_fid_flag = true
        end

    end

    @show xyzgt_flag, pulsetype_flag, hetero_fid_flag, any_acq_flag
    
    for count in eachindex(data)
        currentline = data[count]
        if occursin("#", currentline)
            currentline = currentline[1:findfirst("#", data[count])[1]-1]
        end
        if occursin("J =", currentline)
            ff = findfirst("=", currentline)[1]
            J = parse(Float64, currentline[ff+1:end])
            if sideband_flag == false
                println("J:                    ", J)
            end
        elseif occursin("dwell", currentline)
            ff = findfirst("=", currentline)[1]
            dwell = parse(Float64, currentline[ff+1:end])
            println("dwell:                ", dwell)
        elseif occursin("T2_time", currentline)
            ff = findfirst("=", currentline)[1]
            T2_time = parse(Float64, currentline[ff+1:end])
            println("T2_time:              ", T2_time)
        elseif occursin("npoints", currentline)
            ff = findfirst("=", currentline)[1]
            npoints = parse(Int64, currentline[ff+1:end])
            if simulate_fid_flag == true || avH_flag == true
                println("npoints:              ", npoints)
            end
        elseif occursin("nB1", currentline)
            ff = findfirst("=", currentline)[1]
            nB1 = parse(Int64, currentline[ff+1:end])
            if simulate_fid_flag == true || avH_flag == true
                println("nB1:                  ", nB1)
            end
        elseif occursin("theta", currentline)
            ff = findfirst("=", currentline)[1]
            theta = parse(Float64, currentline[ff+1:end])
            if simulate_fid_flag == true || avH_flag == true
                println("theta:                ", theta)
            end
        elseif occursin("noffs", currentline)
            ff = findfirst("=", currentline)[1]
            noffs = parse(Int64, currentline[ff+1:end])
            println("noffs:                ", noffs)
        elseif occursin("offsrange", currentline)
            ff = findfirst("=", currentline)[1]
            offsrange = split(currentline[ff+1:end], ",")
            offsrange = [parse(Float64, lstrip(rstrip(ioffs))) for ioffs in offsrange]
        elseif occursin("rfmaxseq", currentline)
            ff = findfirst("=", currentline)[1]
            rfmaxseq_temp = split(currentline[ff+1:end])
            rfmaxseq = [parse(Float64, irfmaxseq) for irfmaxseq in rfmaxseq_temp]
            if (!(xyzgt_flag) || pulsetype_flag) && sequence_fid_flag
                println("rfmaxseq:             ", rfmaxseq)
            end
        elseif occursin("rfmaxacq", currentline)
            ff = findfirst("=", currentline)[1]
            rfmaxacq_temp = split(currentline[ff+1:end])
            rfmaxacq = [parse(Float64, irfmaxacq) for irfmaxacq in rfmaxacq_temp]
            if (!(xyzgt_flag) || pulsetype_flag) && (hetero_fid_flag || any_acq_flag)
                println("rfmaxacq:             ", rfmaxacq)
            end
        elseif occursin("rfmaxhom", currentline)
            ff = findfirst("=", currentline)[1]
            rfmaxhom_temp = split(currentline[ff+1:end])
            rfmaxhom = [parse(Float64, irfmaxhom) for irfmaxhom in rfmaxhom_temp]
            if (!(xyzgt_flag) || pulsetype_flag) && (homo_fid_flag || any_homo_flag)
                println("rfmaxhom:             ", rfmaxhom)
            end
        elseif occursin("tpulseseq", currentline)
            ff = findfirst("=", currentline)[1]
            tpulseseq_temp = split(currentline[ff+1:end])
            tpulseseq = [parse(Float64, itpulseseq) for itpulseseq in tpulseseq_temp]./1000
            if (!(xyzgt_flag) || pulsetype_flag) && sequence_fid_flag
                println("tpulseseq:            ", tpulseseq)
            end
        elseif occursin("tpulseacq", currentline)
            ff = findfirst("=", currentline)[1]
            tpulseacq_temp = split(currentline[ff+1:end])
            tpulseacq = [parse(Float64, itpulseacq) for itpulseacq in tpulseacq_temp]./1000
            if (!(xyzgt_flag) || pulsetype_flag) && (hetero_fid_flag || any_acq_flag)
                println("tpulseacq:            ", tpulseacq)
            end
        elseif occursin("tpulsehom", currentline)
            ff = findfirst("=", currentline)[1]
            tpulsehom_temp = split(currentline[ff+1:end])
            tpulsehom = [parse(Float64, itpulsehom) for itpulsehom in tpulsehom_temp]./1000
            if (!(xyzgt_flag) || pulsetype_flag) && (homo_fid_flag || any_homo_flag)
                println("tpulsehom:            ", tpulsehom)
            end
        elseif occursin("pdigits", currentline)
            ff = findfirst("=", currentline)[1]
            pdigits_temp = split(currentline[ff+1:end])
            pdigits = [parse(Int32, ipdigits) for ipdigits in pdigits_temp]
            if pulsetype_flag
                println("pdigits:              ", pdigits)
            end
        elseif occursin("bwdth", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            bwdth[spinnumber] = parse(Float64, currentline[ff+1:end])
            if pulsetype[spinnumber] == "STUD"
                println("bwdth for STUD:       ", bwdth[spinnumber])
            end
        elseif occursin("n (WURST)", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            n[spinnumber] = parse(Float64, currentline[ff+1:end])
            if pulsetype[spinnumber] == "WURST"
                println("n for WURST:          ", n[spinnumber])
            end
        elseif occursin("k (WURST) ", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            k[spinnumber] = parse(Float64, currentline[ff+1:end])
            if pulsetype[spinnumber] == "WURST"
                println("k for WURST:          ", k[spinnumber])
            end
        elseif occursin("thresh", currentline) && avH_flag
            ff = findfirst("=", currentline)[1]
            thresh = parse(Float64, currentline[ff+1:end])
            println("thresh for AvH:       ", thresh)
        elseif occursin("dchunk", currentline) && homo_fid_flag
            ff = findfirst("=", currentline)[1]
            dchunk = parse(Float64, currentline[ff+1:end]) / J
            println("dchunk for HomoD:     ", dchunk)
        elseif occursin("eps", currentline) && waugh_flag
            ff = findfirst("=", currentline)[1]
            eps = parse(Float64, currentline[ff+1:end])
            println("eps for simple Waugh: ", eps)
        elseif occursin("filename", currentline) && bilevel_flag == false
            ff = findfirst("=", currentline)[1]
            pulsefilename = lstrip(rstrip(currentline[ff+1:end]))
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            println("filename for spin ", spinnumber, ":  ", pulsefilename)
            if occursin("xyzgt", pulsefilename)
                sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber], acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = read_xyzgt(pulsefilename)
            else
                if hetero_fid_flag || any_acq_flag
                    acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber] = bruker_to_xyzgt(pulsefilename, rfmaxacq[spinnumber], tpulseacq[spinnumber])
                    # sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                elseif homo_fid_flag || any_homo_flag
                    homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = bruker_to_xyzgt(pulsefilename, rfmaxhom[spinnumber], tpulsehom[spinnumber])
                    # sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber], acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                else
                    sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber] = bruker_to_xyzgt(pulsefilename, rfmaxseq[spinnumber], tpulseseq[spinnumber])
                    # acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                end
            end
        elseif occursin("filename", currentline) && bilevel_flag == true
            ff = findfirst("=", currentline)[1]
            pulsefilename = split(currentline[ff+1:end])
            pulsefilename = [lstrip(rstrip(ifile)) for ifile in pulsefilename]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            println("filenames for bilevelpulse for spin ", spinnumber, ": ", pulsefilename)
            if occursin("xyzgt", pulsefilename[1])     # pulsefilenames muss either all be bruker or xyzgt format!!
                println("bilevel files, xyzgt format")
                sequxm[spinnumber], sequym[spinnumber], sequzm[spinnumber], seqgradm[spinnumber], sequtm[spinnumber], acquxm[spinnumber], acquym[spinnumber], acquzm[spinnumber], acqgradm[spinnumber], acqutm[spinnumber], homuxm[spinnumber], homuym[spinnumber], homuzm[spinnumber], homgradm[spinnumber], homutm[spinnumber] = read_xyzgt_bilevel(pulsefilename)
            else
                println("bilevel files, bruker format")
                if hetero_fid_flag || any_acq_flag
                    acquxm[spinnumber], acquym[spinnumber], acquzm[spinnumber], acqgradm[spinnumber], acqutm[spinnumber] = bruker_to_xyzgt_bilevel(pulsefilename, rfmaxacq[spinnumber], tpulseacq[spinnumber])
                elseif homo_fid_flag || any_homo_flag
                    homuxm[spinnumber], homuym[spinnumber], homuzm[spinnumber], homgradm[spinnumber], homutm[spinnumber] = bruker_to_xyzgt_bilevel(pulsefilename, rfmaxhom[spinnumber], tpulsehom[spinnumber])
                else
                    sequxm[spinnumber], sequym[spinnumber], sequzm[spinnumber], seqgradm[spinnumber], sequtm[spinnumber] = bruker_to_xyzgt_bilevel(pulsefilename, rfmaxseq[spinnumber], tpulseseq[spinnumber])
                end
            end
        elseif occursin("pulsetype", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            if hetero_fid_flag || any_acq_flag
                acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqut[spinnumber] = pulse_for_supercycle(pulsetype[spinnumber], rfmaxacq[spinnumber], pdigits[spinnumber], tpulseacq[spinnumber], bwdth[spinnumber], n[spinnumber], k[spinnumber])
                # sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                acqgrad = zeros(size(acqux[spinnumber]))
            elseif homo_fid_flag || any_homo_flag
                homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homut[spinnumber] = pulse_for_supercycle(pulsetype[spinnumber], rfmaxhom[spinnumber], pdigits[spinnumber], tpulsehom[spinnumber], bwdth[spinnumber], n[spinnumber], k[spinnumber])
                # sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], seqgrad[spinnumber], sequt[spinnumber], acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                homgrad = zeros(size(homux[spinnumber]))
            else
                sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], sequt[spinnumber] = pulse_for_supercycle(pulsetype[spinnumber], rfmaxseq[spinnumber], pdigits[spinnumber], tpulseseq[spinnumber], bwdth[spinnumber], n[spinnumber], k[spinnumber])
                # acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqgrad[spinnumber], acqut[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homgrad[spinnumber], homut[spinnumber] = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                seqgrad = zeros(size(sequx[spinnumber]))
            end
            if bilevel_flag == true && (pulsetype[spinnumber] == "WURST" || pulsetype[spinnumber] == "STUD" )
                acquxm, acquym, acquzm, acqutm = adiabatic_to_bilevel(nbilevel, nspins, acqux, acquy, acquz, acqut)
            end    
        elseif occursin("supercycle", currentline)
            ff = findfirst("=", currentline)[1]
            spinnumber = parse(Int32, split(currentline[1:ff-1])[end])
            supercycle_temp = split(currentline[ff+1:end])
            supercycle[spinnumber] = [string(lstrip(rstrip(icyc))) for icyc in supercycle_temp]
            if hetero_fid_flag || any_acq_flag
                acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqut[spinnumber] = supercycles(supercycle[spinnumber], acqux[spinnumber], acquy[spinnumber], acquz[spinnumber], acqut[spinnumber])
            elseif homo_fid_flag || any_homo_flag
                homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homut[spinnumber] = supercycles(supercycle[spinnumber], homux[spinnumber], homuy[spinnumber], homuz[spinnumber], homut[spinnumber])
            else
                sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], sequt[spinnumber] = supercycles(supercycle[spinnumber], sequx[spinnumber], sequy[spinnumber], sequz[spinnumber], sequt[spinnumber])
            end
            println("supercycle of spin ", spinnumber, ": ", supercycle[spinnumber])
        elseif occursin("rhoinit", currentline)
            ff = findfirst("=", currentline)[1]
            startmag = currentline[ff+1:end]
        # elseif occursin("cs =", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     cs = split(currentline[ff+1:end])
        #     cs = [parse(Float64, ics) for ics in cs]
        #     println("cs:                   ", cs)
        # elseif occursin("sc =", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     sc_split = split(currentline[ff+1:end], ";")
        #     sc = [split(sc_split[isc]) for isc in eachindex(sc_split)]
        #     sc = [parse(Float64, sc[isc1][isc2]) for isc1 in eachindex(sc) for isc2 in eachindex(sc[1])]
        #     sc = reshape(sc, size(sc_split)[1], size(sc_split)[1])'
        #     println("sc:                   ", sc)
        # elseif occursin("dc =", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     dc_split = split(currentline[ff+1:end], ";")
        #     dc = [split(dc_split[idc]) for idc in eachindex(dc_split)]
        #     dc = [parse(Float64, dc[idc1][idc2]) for idc1 in eachindex(dc) for idc2 in eachindex(dc[1])]
        #     dc = reshape(dc, size(dc_split)[1], size(dc_split)[1])'
        #     println("dc:                   ", dc)
        # elseif occursin("spingroups", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     spingroups_temp = split(currentline[ff+1:end])
        #     spingroups = [parse(Int32, ispingroup) for ispingroup in spingroups_temp]
        #     println("spingroups:           ", spingroups)
        # elseif occursin("pulsespingroup", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     pulsespingroup_temp = split(currentline[ff+1:end])
        #     pulsespingroup = [parse(Int32, ipulsespingroup) for ipulsespingroup in pulsespingroup_temp]
        #     println("pulsespingroup:       ", pulsespingroup)
        # elseif occursin("fidspingroup", currentline) && spins10_flag
        #     ff = findfirst("=", currentline)[1]
        #     fidspingroup = parse(Int32, currentline[ff+1:end])
        #     println("fidspingroup:         ", fidspingroup)
        elseif occursin("pulsespins", currentline)
            ff = findfirst("=", currentline)[1]
            pulsespins = [parse(Int32, ipulsespins) for ipulsespins in split(currentline[ff+1:end])]
            println("pulsespins:           ", pulsespins)
        elseif occursin("offsetspins", currentline)
            ff = findfirst("=", currentline)[1]
            offsetspins = [parse(Int32, ioffsetspins) for ioffsetspins in split(currentline[ff+1:end])]
            println("offsetspins:          ", offsetspins)
        elseif occursin("fidspin", currentline)
            ff = findfirst("=", currentline)[1]
            fidspin = parse(Int32, currentline[ff+1:end])
            println("fidspin:              ", fidspin)
        elseif occursin("haminit", currentline) && avH_flag
            if occursin("x", currentline)
                haminit = [1.0, 0.0, 0.0]
                println("haminit:              x")
            elseif occursin("y", currentline)
                haminit = [0.0, 1.0, 0.0]
                println("haminit:              y")
            elseif occursin("z", currentline)
                haminit = [0.0, 0.0, 1.0]
                println("haminit:              z")
            end
        elseif occursin("centered inversion", currentline) && avH_flag
            ff = findfirst("=", currentline)[1]
            centered_inversion = currentline[ff+1:end]
            if centered_inversion == "true"
                centered_inversion = true
                println("centered inversion used")
            else
                centered_inversion = false
                println("centered inversion not used")
            end
        elseif occursin("allJs", currentline) && sideband_flag
            ff = findfirst("=", currentline)[1]
            allJs = split(currentline[ff+1:end])
            allJs = [parse(Float64,iJ) for iJ in allJs]
            println("allJs:                ", allJs)
        elseif occursin("dgrad", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            dgrad = parse(Float64, currentline[ff+1:end])/1000
            println("dgrad:                ", dgrad)
        elseif occursin("gradstrength", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            gradstrength = parse(Float64, currentline[ff+1:end])
            println("gradstrength:         ", gradstrength)
        elseif occursin("ngrads", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            ngrads = parse(Int32, currentline[ff+1:end])
            println("ngrads:               ", ngrads)
        elseif occursin("gradphase", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            gradphase = parse(Int32, currentline[ff+1:end])
            println("gradphase:            ", gradphase)
        elseif occursin("tubesize", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            tubesize = parse(Float64, currentline[ff+1:end])
            println("tubesize:             ", tubesize)
        elseif occursin("gammaH", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            gammaH = parse(Float64, currentline[ff+1:end]) * 1000000
            println("gammaH:               ", gammaH)
        elseif occursin("gradspins", currentline) && gradient_flag
            ff = findfirst("=", currentline)[1]
            gradspins_temp = split(currentline[ff+1:end])
            gradspins = [parse(Int32, igradspins) for igradspins in gradspins_temp]
            println("gradspins:            ", gradspins)
        end
    end
    iu, ix, iy, iz, ip, im, ia, ib = basis(nspins)
    close(io)

    offs = get_offset_asymm(offsrange[1], offsrange[2], noffs)
    println("offsrange:            ", offsrange[1], " to ", offsrange[2], " Hz")
    fB1 = get_B1(nB1, theta)

    return nspins, npoints, nB1, noffs, fB1, offs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, nbilevel, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, thresh, dchunk, eps, startmag, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH, bilevel_flag, iu, ix, iy, iz, ip, im, ia, ib
end

"simulations are performed as defined in the paramfilename\n
output: global variables depending on the performed simulations "
function simulate(paramfilename, nspins, npoints, nB1, noffs, fB1, offs, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, nbilevel::Int64, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, thresh, dchunk, eps, startmag, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH, bilevel_flag, iu, ix, iy, iz, ip, im, ia, ib)
    println("Simulations file ", paramfilename, " will be used for definition of simulations")

    io = open(paramfilename, "r")
    data = readlines(paramfilename)

    for count in eachindex(data)
        currentline = data[count]
        if occursin("#", currentline)
            currentline = currentline[1:findfirst("#", data[count])[1]-1]
        end
        
        if occursin("simulate FID", currentline)
            if occursin("reduced", currentline)
                println("starting magnetization along x is used")
                rhoinit = [0.0, 0.0, 0.0, 1.0]
            else
                if occursin("x", startmag)
                    rhoinit = sum(ix[:,:,:], dims=3)[:,:,1]
                    println("starting magnetization along x is used")
                elseif occursin("y", startmag)
                    rhoinit = sum(iy[:,:,:], dims=3)[:,:,1]
                    println("starting magnetization along y is used")
                elseif occursin("z", startmag)
                    rhoinit = sum(iz[:,:,:], dims=3)[:,:,1]
                    println("starting magnetization along z is used")
                end
            end
            if occursin("Hetero", currentline) && !(occursin("Homo", currentline)) && !(occursin("Sequence", currentline))
                if occursin("reduced", currentline)
                    if occursin("weak", currentline) && bilevel_flag == false
                        println("FID is simulated with weak coupling in the reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD("weak", "reduced", acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("strong", currentline) && bilevel_flag == false
                        println("FID is simulated with strong coupling in the reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD("strong", "reduced", acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("weak", currentline) && bilevel_flag == true
                        println("bilevel FID is simulated with weak coupling in the reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD_bilevel("weak", "reduced", acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, nspins, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("strong", currentline) && bilevel_flag == true
                        println("bilevel FID is simulated with strong coupling in the reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD_bilevel("strong", "reduced", acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, nspins, offs, fidspin, pulsespins, offsetspins)
                    end
                else
                    if occursin("weak", currentline) && bilevel_flag == false
                        println("FID is simulated with weak coupling in the NOT reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD("weak", "not reduced", acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("strong", currentline) && bilevel_flag == false
                        println("FID is simulated with strong coupling in the NOT reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD("strong", "not reduced", acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("weak", currentline) && bilevel_flag == true
                        println("bilevel FID is simulated with weak coupling in the NOT reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD_bilevel("weak", "not reduced", acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, nspins, offs, fidspin, pulsespins, offsetspins)
                    elseif occursin("strong", currentline) && bilevel_flag == true
                        println("bilevel FID is simulated with strong coupling in the NOT reduced Liouville space")
                        global spectr, fid, projection, tim = @time simulate_fid_heteroD_bilevel("strong", "not reduced", acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, nspins, offs, fidspin, pulsespins, offsetspins)
                    end
                end
            elseif occursin("Homo", currentline) && !(occursin("Hetero", currentline)) && !(occursin("Sequence", currentline))
                if occursin("weak", currentline) && bilevel_flag == false
                    println("Homo FID is simulated with weak coupling in the NOT reduced Liouville space")
                    global spectr, fid, projection, tim, sfid = @time simulate_fid_homoD("weak", "not reduced", false, homux, homuy, homuz, homut, rhoinit, dchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == false
                    println("Homo FID is simulated with strong coupling in the NOT reduced Liouville space")
                    global spectr, fid, projection, tim, sfid = @time simulate_fid_homoD("strong", "not reduced", false, homux, homuy, homuz, homut, rhoinit, dchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
                end
            elseif occursin("Hetero", currentline) && occursin("Homo", currentline) && !(occursin("Sequence", currentline))
                if occursin("weak", currentline) && bilevel_flag == false
                    println("FID is simulated with weak coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, tim, sfid = @time simulate_fid_heteroD_homoD("weak", acqux, acquy, acquz, acqut, homux, homuy, homuz, homut, rhoinit, olddchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == false
                    println("FID is simulated with strong coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, tim, sfid = @time simulate_fid_heteroD_homoD("strong", acqux, acquy, acquz, acqut, homux, homuy, homuz, homut, rhoinit, dchunk, J, offs, fB1, nB1, noffs, npoints, dwell, T2_time, fidspin, pulsespins, offsetspins)
                end
            elseif occursin("Sequence", currentline) && !(occursin("Hetero", currentline)) && !(occursin("Homo", currentline))
                if occursin("weak", currentline) && bilevel_flag == false
                    println("FID is simulated with weak coupling in the NOT reduced Liouville space with a pulsesequence before the FID")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence("weak", "not reduced", rhoinit, sequx, sequy, sequz, sequt, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == false
                    println("FID is simulated with strong coupling in the NOT reduced Liouville space with a pulsesequence before the FID")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence("strong", "not reduced", rhoinit, sequx, sequy, sequz, sequt, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("weak", currentline) && bilevel_flag == true
                    println("bilevel FID is simulated with weak coupling in the NOT reduced Liouville space with a pulsesequence before the FID")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_bilevel("weak", rhoinit, sequxm, sequym, sequzm, sequtm, nbilevel, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == true
                    println("bilevel FID is simulated with strong coupling in the NOT reduced Liouville space with a pulsesequence before the FID")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_bilevel("strong", rhoinit, sequxm, sequym, sequzm, sequtm, nbilevel, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                end
            elseif occursin("Sequence", currentline) && occursin("Hetero", currentline) && !(occursin("Homo", currentline))
                if occursin("weak", currentline) && bilevel_flag == false
                    println("FID is simulated with weak coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_heteroD("weak", "not reduced", sequx, sequy, sequz, sequt, acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == false
                    println("FID is simulated with strong coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_heteroD("strong", "not reduced", sequx, sequy, sequz, sequt, acqux, acquy, acquz, acqut, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("weak", currentline) && bilevel_flag == true
                    println("bilevel FID is simulated with weak coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_heteroD_bilevel("weak", sequxm, sequym, sequzm, sequtm, acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                elseif occursin("strong", currentline) && bilevel_flag == true
                    println("bilevel FID is simulated with strong coupling in the NOT reduced Liouville space with a pulsesequence before the FID and heteronuclear decoupling")
                    global spectr, fid, projection, time = @time simulate_fid_pulsesequence_heteroD_bilevel("strong", sequxm, sequym, sequzm, sequtm, acquxm, acquym, acquzm, acqutm, nbilevel, rhoinit, J, fB1, nB1, noffs, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
                end
            end
        end

        if occursin("simulate FID 10 spins", currentline)
            if occursin("weak", currentline)
                println("FID is simulated with weak coupling for up to 10 spins")       # TODO not working yet
                global spectr, fid, projection, tim = @time simulate_fid_10spins("weak", nspins, spingroups, pulsespingroup, fidspingroup, cs, sc, dc, acqux, acquy, acquz, acqut, rhoinit, fB1, nB1, noffs, npoints, dwell, T2_time)
            elseif occursin("strong", currentline)
                println("FID is simulated with strong coupling for up to 10 spins")       # TODO not working yet
                global spectr, fid, projection, tim = @time simulate_fid_10spins("strong", nspins, spingroups, pulsespingroup, fidspingroup, cs, sc, dc, acqux, acquy, acquz, acqut, rhoinit, fB1, nB1, noffs, npoints, dwell, T2_time)
            end
        elseif occursin("simulate FID gradients", currentline)
            # if occursin("reduced", currentline)
            #     if occursin("weak", currentline)
            #         println("FID with gradients is simulated with weak coupling in the reduced Liouville space")
            #         global spectr, fid, allfids, projection, time = simulate_fid_grad("weak", "reduced", ux, uy, uz, ut, gammaH, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, alignment, gradspins, gradstrength, ngrads, dgrad, gradphase, tubesize)
            #     elseif occursin("strong", currentline)
            #         println("FID with gradients is simulated with strong coupling in the reduced Liouville space")
            #         global spectr, fid, allfids, projection, time = simulate_fid_grad("strong", "reduced", ux, uy, uz, ut, gammaH, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, alignment, gradspins, gradstrength, ngrads, dgrad, gradphase, tubesize)
            #     end
            # else
                if occursin("weak", currentline)
                    println("FID with gradients is simulated with weak coupling in the NOT reduced Liouville space")       # TODO not working yet
                    global spectr, fid, allspectr, allfids, projection, tim = @time simulate_fid_grad("weak", "not reduced", acqux, acquy, acqgrad, acquz, acqut, gamma, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, ngrads, tubesize)
                elseif occursin("strong", currentline)
                    println("FID with gradients is simulated with strong coupling in the NOT reduced Liouville space")       # TODO not working yet
                    global spectr, fid, allspectr, allfids, projection, tim = @time simulate_fid_grad("strong", "not reduced", acqux, acquy, acqgrad, acquz, acqut, gamma, J, rhoinit, npoints, dwell, noffs, offs, T2_time, nB1, ngrads, tubesize)
                end
            # end
        elseif occursin("simulate sidebands", currentline)
            if !bilevel_flag
                println("FIDs are simulated with weak coupling in the reduced Liouville space for different Js")       # TODO not working yet
                global maxintens, maxatfreq, maxintensatoffs, maxintensside, maxsideatfreq, maxsideatoffs, maxintenssidesub, maxsideatfreqsub, maxsidesubatoffs, allspectr = @time simulate_sidebands_J(acqux, acquy, acquz, acqut, allJs, nB1, noffs, fB1, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins)
            else
                println("FIDs are simulated with weak coupling in the reduced Liouville space for different Js")       # TODO not working yet
                global maxintens, maxatfreq, maxintensside, maxsideatfreq, maxintenssidesub, maxsideatfreqsub, allspectr = @time simulate_sidebands_J_bilevel(nbilevel, acquxm, acquym, acquzm, acqutm, allJs, nB1, noffs, fB1, npoints, dwell, T2_time, offs, fidspin, pulsespins, offsetspins, nspins)
            end
        end
        if occursin("xi", currentline)      # only for acq at the moment
            println("Calculation of xi")
            if bilevel_flag == false
                global rfav_power = get_rf_average_power(acqux, acquy, acqut, pulsespins)
            else
                global rfav_power = get_rf_average_power_bilevel(acquxm, acquym, acqutm, pulsespins, nbilevel)
            end
            threshold = parse(Float64, split(currentline)[2])
            global xi, bandwidth = determine_xi(threshold, rfav_power, projection, offs)
        end
        if occursin("simulate average H", currentline)
            println("Average Hamiltonian is simulated")
            global hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime, = @time simulate_average_H(sequx, sequy, sequz, sequt, centered_inversion, haminit, thresh, fB1, nB1, noffs, dwell)
            # global ham2x, ham2y, ham2z, hamdwell2x, hamdwell2y, hamdwell2z, ptgrid, dwelltime, = @time simulate_average_H(sequx2, uy2, uz2, ut2, rfmax2, centered_inversion2, haminit, thresh)
        elseif occursin("simulate repeated pulse average H", currentline)
            println("Average Hamiltonian of repeated pulse is simulated")
            global hamx, hamy, hamz, hamdwellx, hamdwelly, hamdwellz, ptgrid, dwelltime, dtgrid, = @time simulate_rep_average_H(sequx, sequy, sequz, sequt, centered_inversion, haminit, thresh, npoints, dwell, fB1, nB1, noffs)
            # global ham2x, ham2y, ham2z, hamdwell2x, hamdwell2y, hamdwell2z, ptgrid, dwelltime, dtgrid, = simulate_rep_average_H(sequx2, uy2, uz2, ut2, rfmax2, centered_inversion2, haminit, thresh)
        end
        if occursin("sum H", currentline)
            println("average over whole dwell time is calculated")
            global sumhamdwellx, sumhamdwelly, sumhamdwellz = sum_ham_over_dw(hamdwellx, hamdwelly, hamdwellz, nB1, noffs)
            # global sumhamdwell2x, sumhamdwell2y, sumhamdwell2z = sum_ham_over_dw(hamdwell2x, hamdwell2y, hamdwell2z)
        end
        if occursin("Waugh simple", currentline)
            pulsenumber = parse(Int32,split(currentline)[end])
            if occursin("acq", currentline)
                println("Simulation of simple Waugh criterion of acq pulse of spin ", pulsenumber)
                global lambda, phip, phim = @time simulate_waugh_simple(acqux[pulsenumber], acquy[pulsenumber], acquz[pulsenumber], acqut[pulsenumber], eps, noffs, offs)
            elseif occursin("hom", currentline)
                println("Simulation of simple Waugh criterion of hom pulse of spin ", pulsenumber)
                global lambda, phip, phim = @time simulate_waugh_simple(homux[pulsenumber], homuy[pulsenumber], homuz[pulsenumber], homut[pulsenumber], eps, noffs, offs)
            else
                println("Simulation of simple Waugh criterion of seq pulse of spin ", pulsenumber)
                global lambda, phip, phim = @time simulate_waugh_simple(sequx[pulsenumber], sequy[pulsenumber], sequz[pulsenumber], sequt[pulsenumber], eps, noffs, offs)
            end
        end
        if occursin("Waugh exact", currentline)
            println("Simulation of exact Waugh criterion")
            global J1, J2, I1, I2, phip, phim, beffp, beffm, toobig = @time simulate_waugh_exact(true, acqux, acquy, acquz, acqut, noffs, offs, J)
        end
        if occursin("simulate magnetization", currentline)
            pulsenumber = parse(Int32,split(currentline)[end])
            println("x, y and z Magnetization of pulse of spin ", pulsenumber, " is simulated, starting with magnetization along z")
            if occursin("acq", currentline)
                global Mx, My, Mz =  magnetization(acqux[pulsenumber], acquy[pulsenumber], acqut[pulsenumber], nB1, noffs, fB1, offs)
            elseif occursin("hom", currentline)
                global Mx, My, Mz =  magnetization(homux[pulsenumber], homuy[pulsenumber], homut[pulsenumber], nB1, noffs, fB1, offs)
            else
                global Mx, My, Mz =  magnetization(sequx[pulsenumber], sequy[pulsenumber], sequt[pulsenumber], nB1, noffs, fB1, offs)
            end
        end

    end

    # NOT IMPLEMENTED YET
    # simulate_fid_grad("weak", "not reduced", ux, uy, uz, ut, graduz, gradut, gamma, rhoinit)


    close(io)
    return 
end


"types of all necessary variables are defined\n
output: npoints, nB1, noffs, pdigits, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, tpulseseq, tpulseacq, tpulsehom, bwdth, n, k, thresh, dchunk, eps, pulsefilename, pulsetype, supercycle, startmag, offsrange, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH"
function initializeallparams(nspins)
    npoints::Int32, nB1::Int32, noffs::Int32, fidspin::Int32, fidspingroup::Int32, ngrads::Int32, gradphase::Int32 = 0, 0, 0, 0, 0, 0, 0, 0
    J::Float64, dwell::Float64, T2_time::Float64, theta::Float64, thresh::Float64, dchunk::Float64, eps::Float64, dgrad::Float64, gradstrength::Float64, tubesize::Float64, gammaH::Float64 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    pulsefilename::String, startmag::String = "", ""
    offsrange = @MVector zeros(Float64, 2)
    haminit = @MVector zeros(Float64, 3)
    rfmaxseq, rfmaxacq, rfmaxhom = zeros(Float64, nspins), zeros(Float64, nspins), zeros(Float64, nspins)
    tpulseseq, tpulseacq, tpulsehom = zeros(Float64, nspins), zeros(Float64, nspins), zeros(Float64, nspins)
    pdigits = zeros(Int32, nspins)
    bwdth = zeros(Float64, nspins)
    n = zeros(Float64, nspins)
    k = zeros(Float64, nspins)
    centered_inversion = false
    pulsetype = String[]
    cs, sc, dc, spingroups, pulsespingroup, pulsespins, offsetspins, allJs = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(),Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm =  Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}(), Vector{Vector{Vector{Float64}}}()
    supercycle = Vector{Vector{String}}()
    for ispin = 1:nspins
        push!(sequx, []); push!(sequy, []); push!(sequz, []); push!(seqgrad, []); push!(sequt, [])
        push!(acqux, []); push!(acquy, []); push!(acquz, []); push!(acqgrad, []); push!(acqut, [])
        push!(homux, []); push!(homuy, []); push!(homuz, []); push!(homgrad, []); push!(homut, [])
        push!(sequxm, []); push!(sequym, []); push!(sequzm, []); push!(seqgradm, []); push!(sequtm, [])
        push!(acquxm, []); push!(acquym, []); push!(acquzm, []); push!(acqgradm, []); push!(acqutm, [])
        push!(homuxm, []); push!(homuym, []); push!(homuzm, []); push!(homgradm, []); push!(homutm, [])
        push!(supercycle, [])
        push!(pulsetype, "")
    end
    
    gradspins = Int32[]
    
    return npoints, nB1, noffs, pdigits, sequx, sequy, sequz, seqgrad, sequt, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut, sequxm, sequym, sequzm, seqgradm, sequtm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm, J, dwell, T2_time, theta, rfmaxseq, rfmaxacq, rfmaxhom, tpulseseq, tpulseacq, tpulsehom, bwdth, n, k, thresh, dchunk, eps, pulsefilename, pulsetype, supercycle, startmag, offsrange, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup, fidspin, pulsespins, offsetspins, ngrads, gradphase, dgrad, gradspins, gradstrength, tubesize, gammaH
end




# ----------------- OLD FUNCTIONS ------------------


# function order_by_spinnumber(spinnumber, nspins, uxold, uyold, uzold, gradold, utold, acquxold, acquyold, acquzold, acqgradold, acqutold, homuxold, homuyold, homuzold, homgradold, homutold)
    
#     ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(),Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()


#     for ispin = 1:nspins
#         if !(ispin in spinnumber)
#             append!(spinnumber, ispin)
#             push!(uxold, [])
#             push!(uyold, [])
#             push!(uzold, [])
#             push!(gradold, [])
#             push!(utold, [])
#             push!(acquxold, [])
#             push!(acquyold, [])
#             push!(acquzold, [])
#             push!(acqgradold, [])
#             push!(acqutold, [])
#             push!(homuxold, [])
#             push!(homuyold, [])
#             push!(homuzold, [])
#             push!(homgradold, [])
#             push!(homutold, [])
#         end
#         currentindex = findfirst(isequal(ispin), spinnumber)
#         push!(ux, uxold[currentindex])
#         push!(uy, uyold[currentindex])
#         push!(uz, uzold[currentindex])
#         push!(grad, gradold[currentindex])
#         push!(ut, utold[currentindex])
#         push!(acqux, acquxold[currentindex])
#         push!(acquy, acquyold[currentindex])
#         push!(acquz, acquzold[currentindex])
#         push!(acqgrad, acqgradold[currentindex])
#         push!(acqut, acqutold[currentindex])
#         push!(homux, homuxold[currentindex])
#         push!(homuy, homuyold[currentindex])
#         push!(homuz, homuzold[currentindex])
#         push!(homgrad, homgradold[currentindex])
#         push!(homut, homutold[currentindex])
#     end

#     return ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut
# end





# "### initializeallresults()
# types of all necessary variables are defined"
# function initializeallresults()
#     nspins::Int32, npoints::Int32, nB1::Int32, noffs::Int32, pdigits::Int32, fidspingroup::Int32 = 0, 0, 0, 0, 0, 0
#     J::Float64, dwell::Float64, T2_time::Float64, theta::Float64, rfmax::Float64, tpulse::Float64, bwdth::Float64, n::Float64, k::Float64, thresh::Float64, dchunk::Float64, eps::Float64 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#     pulsefilename::String, pulsetype::String, startmag::String = "", "", ""
#     offsrange = zeros(Float64, 2)   # TODO StaticArray
#     haminit = zeros(Float64, 3)
#     centered_inversion = false
#     supercycles = String[]
#     cs, sc, dc, spingroups, pulsespingroup, allJs = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

#     # misisng: supercycle, cs, sc, dc, spingroups, pulsespingroup, fidspingroup, allJs

#     return nspins, npoints, nB1, noffs, pdigits, J, dwell, T2_time, theta, rfmax, tpulse, bwdth, n, k, thresh, dchunk, eps, pulsefilename, pulsetype, startmag, offsrange, haminit, centered_inversion, cs, sc, dc, spingroups, pulsespingroup, allJs, fidspingroup
# end
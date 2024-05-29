
"read in ux, uy, uz and ut of xyzt file format\n
output: ux, uy, uz, ut"
function read_xyzt(filename)

    ux = Float64[]
    uy = Float64[]
    uz = Float64[]
    ut = Float64[]

    fi = open(filename,"r")
    data = readlines(filename)
    for line in data
        if !occursin("#", line) 
            a,b,c,d = map(x->parse(Float64,x),collect(eachsplit(line)))
            append!(ux,a)
            append!(uy,b)
            append!(uz,c)
            append!(ut,d)
        end
    end

    close(fi)

    return ux, uy, uz, ut
end


"read in ux, uy and ut of xyt file format\n
output: ux, uy, ut"
function read_xyt(filename)

    ux = Float64[]
    uy = Float64[]
    ut = Float64[]

    open(filename,"r")
    data = readlines(filename)
    for line in data
        if !occursin("#", line) 
            a,b,c = map(x->parse(Float64,x),collect(eachsplit(line)))
            append!(ux,a)
            append!(uy,b)
            append!(ut,c)
        end
    end

    return ux, uy, ut
end


"read in ux, uy,grad and ut of xyzgt file format\n
output: ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut"
function read_xyzgt(filename)

    ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

    seqflag = false
    acqflag = false
    homodflag = false

    fi = open(filename,"r")
    data = readlines(filename)
    for line in data
        if !occursin("#", line) 

            if occursin("SEQ=start", line)
                seqflag = true
            elseif occursin("SEQ=end", line)
                seqflag = false
            elseif occursin("ACQ=start", line)
                acqflag = true
            elseif occursin("ACQ=end", line)
                acqflag = false
            elseif occursin("HOMOD=start", line)
                homodflag = true
            elseif occursin("HOMOD=end", line)
                homodflag = false
            else

                a,b,c,d,e = map(x->parse(Float64,x),collect(eachsplit(line)))

                if seqflag
                    append!(ux,a); append!(uy,b); append!(uz,c); append!(grad,d); append!(ut,e)
                end
                if acqflag
                    append!(acqux,a); append!(acquy,b); append!(acquz,c); append!(acqgrad,d); append!(acqut,e)
                end
                if homodflag
                    append!(homux,a); append!(homuy,b); append!(homuz,c); append!(homgrad,d); append!(homut,e) 
                end
            end
        end
    end

    close(fi)

    return ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut
end


"read in bruker file and return x, y, z and time of each pulse as vector of arrays for multiple FIDs\n
output: uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm"
function read_xyzgt_bilevel(filename)
    
    uxm, uym, uzm, gradm, utm = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    acquxm, acquym, acquzm, acqgradm, acqutm = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    homuxm, homuym, homuzm, homgradm, homutm = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]


    for ifid = 1:size(filename)[1]
        ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut = read_xyzgt(filename[ifid])
        push!(uxm, ux); push!(uym, uy); push!(uzm, uz); push!(gradm, grad); push!(utm, ut)
        push!(acquxm, acqux); push!(acquym, acquy); push!(acquzm, acquz); push!(acqgradm, acqgrad); push!(acqutm, acqut)
        push!(homuxm, homux); push!(homuym, homuy); push!(homuzm, homuz); push!(homgradm, homgrad); push!(homutm, homut)
    end
    return uxm, uym, uzm, gradm, utm, acquxm, acquym, acquzm, acqgradm, acqutm, homuxm, homuym, homuzm, homgradm, homutm
end


"read in ampl and phase from Bruker file format\n
output: ampl, phase"
function read_bruker_wave(filename)

    ampl = Float64[]
    phase = Float64[]

    fi = open(filename,"r")
    data = readlines(filename)
    for line in data
        if !occursin("#", line) 
            a,b = map(x->parse(Float64,x),collect(eachsplit(line,",")))
            append!(ampl,a)
            append!(phase,b)
        end
    end

    close(fi)

    return ampl, phase
end


"read in bruker file and return x, y, z and time of each pulse\n
output: ux, uy, uz, grad, ut"
function bruker_to_xyzgt(filename, rfmax, tpulse)

    ampl, phase = read_bruker_wave(filename)
    ampl = ampl*rfmax/100.0
    ux, uy = ampl_phase_to_xy(ampl, phase)
    ut = ones(size(ampl)[1])
    ut = ut .* (tpulse ./ size(ampl)[1])
    uz = zeros(size(ux)[1])
    grad = zeros(size(ux)[1])

    return ux, uy, uz, grad, ut
end

"read in bruker file and return x, y, z and time of each pulse as vector of arrays for multiple FIDs\n
output: fux, fuy, fuz, fgrad, fut"
function bruker_to_xyzgt_bilevel(filename, rfmax, tpulse)
    
    fux, fuy, fuz, fgrad, fut = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]

    for ifid = 1:size(filename)[1]
        ampl, phase = read_bruker_wave(filename[ifid])
        ampl = ampl*rfmax/100.0
        ux, uy = ampl_phase_to_xy(ampl, phase)
        ut = ones(size(ampl)[1])
        ut = ut .* (tpulse ./ size(ampl)[1])
        uz = zeros(size(ux)[1])
        grad = zeros(size(ux)[1])
        push!(fux, ux); push!(fuy, uy); push!(fuz, uz); push!(fgrad, grad); push!(fut, ut)
    end
    return fux, fuy, fuz, fgrad, fut
end



function test_read_xyzgt()
    write_xyzgt_wave("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test.xyzgt", "# test pulse \n# having gradient, homo and hetero Decoupling", [1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [2], [1], [0], [0], [0.03])
    write_xyzgt_wave("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_heteroD.xyzgt", "# test pulse \n# having gradient, homo Decoupling", [1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [], [], [], [], [], [2], [1], [0], [0], [0.03])
    write_xyzgt_wave("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_homoD.xyzgt", "# test pulse \n# having gradient, hetero Decoupling", [1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [], [], [], [], [])
    write_xyzgt_wave("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_sequence.xyzgt", "# test pulse \n# having gradient, homo and hetero Decoupling, no sequence before acquisition", [], [], [], [], [], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [2], [1], [0], [0], [0.03])

    @test ([1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [2], [1], [0], [0], [0.03]) == read_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test.xyzgt")
    @test ([1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [], [], [], [], [], [2], [1], [0], [0], [0.03]) == read_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_heteroD.xyzgt")
    @test ([1,2,3], [4,5,6], [0,0,0], [7,7,7], [0.5, 0.5, 0.5], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [], [], [], [], []) == read_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_homoD.xyzgt")
    @test ([], [], [], [], [], [8,8,8], [1,1,1], [0,0,0], [9,9,9], [0.1, 0.1, 0.2], [2], [1], [0], [0], [0.03]) == read_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/test_no_sequence.xyzgt")
end

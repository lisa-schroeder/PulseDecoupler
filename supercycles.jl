
# Packages 
using LinearAlgebra
using Plots
using FFTW
using Statistics
include("converter.jl")


"Pulses can be generated, e.g. GARP, WALTZ4, WALTZ8, WALTZ16, MLEV, DIPSI1, DIPSI2, DIPSI3, SUSAN, STUD, WURST, F1, F2, ..., F10\n
output: ux, uy, uz, ut"
function pulse_for_supercycle(cycle, rfmax, pdigits, tpulse, bwdth, n, k)
    angle = [0.0]
    phase = [0.0]
    ampl = [0.0]
    ux = Float64[]
    uy = Float64[]
    uz = Float64[]
    ut = Float64[]
    # all pulses generated using angle and phase
    if cycle == "GARP" || cycle == "WALTZ4" || cycle == "WALTZ8" || cycle == "WALTZ16" || cycle == "MLEV" || cycle == "DIPSI1" || cycle == "DIPSI2" || cycle == "DIPSI3" || cycle == "SUSAN"
        if cycle == "GARP"
            angle = [30.5,  55.2, 257.8, 268.3, 69.3,  62.2, 85.0,  91.8, 134.5, 256.1, 66.4,  45.9, 25.5,  72.7, 119.5, 138.2, 258.4,  64.9, 70.9,  77.2, 98.2, 133.6, 255.9,  65.6, 53.4]
            phase = [ 0.0, 180.0,   0.0, 180.0,  0.0, 180.0,  0.0, 180.0,   0.0, 180.0,  0.0, 180.0,  0.0, 180.0,   0.0, 180.0,   0.0, 180.0,  0.0, 180.0,  0.0, 180.0,   0.0, 180.0,  0.0]
        elseif cycle == "WALTZ4"
            angle = [90, 180, 270]
            phase = [0, 180, 0]
        elseif cycle == "WALTZ8"
            angle = [180, 360, 180, 270, 90]
            phase = [180, 0, 180, 0, 180]
        elseif cycle == "WALTZ16"
            angle = [270, 360, 180, 270, 90, 180, 360, 180, 270]
            phase = [180, 0, 180, 0, 180, 0, 180, 0, 180]
        elseif cycle == "MLEV"
            angle = [90, 270, 90]
            phase = [0, 90, 0]
        elseif cycle == "DIPSI1"
            angle = [365, 295, 65, 305, 350]
            phase = [0, 180, 0, 180, 0]
        elseif cycle == "DIPSI2"
            angle = [320, 410, 290, 285, 30, 245, 375, 265, 370]
            phase = [0, 180, 0, 180, 0, 180, 0, 180, 0]
        elseif cycle == "DIPSI3"
            angle = [245, 395, 250, 275, 30, 230, 360, 245, 370, 340, 350, 260, 270, 30, 225, 365, 255, 395]
            phase = [0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0]
        elseif cycle == "SUSAN"
            angle = [27.7, 59.7, 37.6, 17.7, 41.1, 80.0, 43.7, 34.3, 68.1, 81.7, 60.5, 49.3, 0.1, 37.1, 110.3, 163.4, 66.2, 110.5, 0.81, 145.5, 148.0]
            phase = [0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0,180, 0, 180, 0, 180, 0]
        end

        ux, uy, ut = angle_phase_to_xyt(angle, phase, rfmax)
        uz = zeros(size(ux)[1])
    end

    # all pulses generated using ampl and phase
    if cycle == "STUD" || cycle == "WURST"
        if cycle == "STUD"
            ampl = Float64[]
            phase = Float64[]
            ut = ones(pdigits)*tpulse/pdigits
            for idig = 1:pdigits
                append!(ampl, rfmax*sech(5.3(2*idig/pdigits-1)))
                append!(phase, (360*tpulse/10.6)*(bwdth/2)*log(cosh(5.3(2*idig/pdigits-1))))
            end

        elseif cycle == "WURST"
            ampl = rfmax*(1 .- abs.(sin.(([1:pdigits;].*pi./pdigits).-pi/2)).^n)
            phase = (mod.((360 .* k/2)*(([0.0:tpulse/pdigits:tpulse;].-tpulse/2).^2),360))
            ut = ones(pdigits)*tpulse/pdigits
        end

        ux, uy = ampl_phase_to_xy(ampl, phase)
        uz = zeros(size(ux)[1])
    end

    # TPG-n
    if cycle == "F1" || cycle == "F2" || cycle == "F3" || cycle == "F4" || cycle == "F5" || cycle == "F6" || cycle == "F7" || cycle == "F8" || cycle == "F9" || cycle == "F10"
        if cycle == "F1"
            angle = [180]
            uz = [0]
        elseif cycle == "F2"
            angle = [170, 170]
            uz = [-0.92, 0.92]
        elseif cycle == "F3"
            angle = [180, 180, 180]
            uz = [-1.55, 0, 1.55]
        elseif cycle == "F4"
            angle = [166, 198, 198, 166]
            uz = [-2.2, -0.66, 0.66, 2.2]
        elseif cycle == "F5"
            angle = [175, 190, 195, 190, 175]
            uz = [-2.9, -1.4, 0, 1.4, 2.9]
        elseif cycle == "F6"
            angle = [166, 195, 182, 182, 195, 166]
            uz = [-3.6, -2.06, -0.66, 0.66, 2.06, 3.6]
        elseif cycle == "F7"
            angle = [168, 190, 192, 174, 192, 190, 168]
            uz = [-4.27, -2.71, -1.37, 0, 1.37, 2.71, 4.27]
        elseif cycle == "F8"
            angle = [172, 184, 186, 181, 181, 186, 184, 172]
            uz = [-5.09, -3.5, -2.1, -0.71, 0.71, 2.1, 3.5, 5.09]
        elseif cycle == "F9"
            angle = [174, 191, 185, 177, 188, 177, 185, 191, 174]
            uz = [-5.69, -4.19, -2.79, -1.41, 0, 1.41, 2.79, 4.19, 5.69]
        elseif cycle == "F10"
            angle = [181, 185, 181, 180, 181, 181, 180, 181, 185, 181]
            uz = [-6.45, -5.05, -3.55, -2.12, -0.70, 0.70, 2.12, 3.55, 5.05, 6.45]
        end

        phase = zeros(size(angle)[1])
        ux, uy, ut = angle_phase_to_xyt(angle, phase, rfmax)
        uz = uz .* rfmax
    end

    if !(cycle == "GARP" || cycle == "WALTZ4" || cycle == "WALTZ8" || cycle == "WALTZ16" || cycle == "MLEV" || cycle == "DIPSI1" || cycle == "DIPSI2" || cycle == "DIPSI3" || cycle == "SUSAN" || cycle == "STUD" || cycle == "WURST" || cycle == "F1" || cycle == "F2" || cycle == "F3" || cycle == "F4" || cycle == "F5" || cycle == "F6" || cycle == "F7" || cycle == "F8" || cycle == "F9" || cycle == "F10")
        println("this pulse ", cycle, " is not implemented yet")
    end

    return ux, uy, uz, ut
end


"generate ux, uy, uz and ut arrays, starting from a given pulse in ux, uy, uz and ut and the name of the supercycle\n
output: ux, uy, uz, ut"
function supercycles(cycle, ux, uy, uz, ut)
    ampl, phase = xy_to_ampl_phase(ux, uy)
    fphase = Float64[]
    fampl = Float64[]
    fut = Float64[]
    fuz = Float64[]

    # get phases and amplitudes of the supercycles
    cycphase = [get_cycphase(cyc) for cyc in cycle]

    # dependend on the amount of nested supercycles: implemented for up to three
    # calculate final phase, amplitude, uz and ut
    if size(cycle)[1] == 1
        for dig in eachindex(cycphase[1])
            append!(fphase, mod.(phase .+ cycphase[1][dig], 360))
            append!(fampl, ampl)
            append!(fuz, uz)
            append!(fut, ut)
        end
    elseif size(cycle)[1] == 2
        for dig1 in eachindex(cycphase[1])
            for dig2 in eachindex(cycphase[2])
                append!(fphase, mod.(phase .+ cycphase[1][dig1] .+ cycphase[2][dig2], 360))
                append!(fampl, ampl)
                append!(fuz, uz)
                append!(fut, ut)
            end
        end
    elseif size(cycle)[1] == 3
        for dig1 in eachindex(cycphase[1])
            for dig2 in eachindex(cycphase[2])
                for dig3 in eachindex(cycphase[3])
                    append!(fphase, mod.(phase .+ cycphase[1][dig1] .+ cycphase[2][dig2] .+ cycphase[3][dig3], 360))
                    append!(fampl, ampl)
                    append!(fuz, uz)
                    append!(fut, ut)
                end
            end
        end
    else
        println("only up to 3 nested supercycles implemented")
    end

    fux, fuy = ampl_phase_to_xy(fampl, fphase)
    return fux, fuy, fuz, fut
end


"supercycles are applied, e.g. MLEV4, MLEV8, MLEV16, MLEV64, WALTZ4, WALTZ8, WALTZ16, P5, P7, P9\n
output: cycphase"
function get_cycphase(cycle)
    if cycle == "MLEV4" || cycle == "WALTZ4"
        cycphase = [0, 0, 180, 180]
    elseif cycle == "WALTZ8" || cycle == "WALTZ16"
        cycphase = [0, 180, 180, 0]
    elseif cycle == "MLEV8"
        cycphase = [0, 0, 180, 180, 180, 0, 0, 180]
    elseif cycle == "MLEV16"
        cycphase = [0, 0, 180, 180, 180, 0, 0, 180, 180, 180, 0, 0, 0, 180, 180, 0]           
    elseif cycle == "MLEV64"
        cycphase = [0, 0, 180, 180, 180, 0, 0, 180, 180, 180, 0, 0, 0, 180, 180, 0, 0, 0, 0, 180, 180, 180, 0, 0, 180, 180, 180, 0, 0, 0, 180, 180, 180, 180, 0, 0, 0, 180, 180, 0, 0, 0, 180, 180, 180, 0, 0, 180, 180, 180, 180, 0, 0, 0, 180, 180, 0, 0, 0, 180, 180, 180, 0, 0]           # MLEV-64
    elseif cycle == "P5"
        cycphase = [0, 150, 60, 150, 0]
    elseif cycle == "P7"
        cycphase = [0, 105, 300, 255, 300, 105, 0]
    elseif cycle == "P9"
        cycphase = [0, 15, 180, 165, 270, 165, 180, 15, 0]
    else
        println("supercycle not implemented yet")
    end
    return cycphase
end
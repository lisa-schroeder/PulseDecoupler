

"calculate rotation angle of pulse\n
uevop: evolution propagator for offset + J/2\n
uevom: evolution propagator for offset - J/2\n
output: uevop, uevom"
function get_waugh_rotation_simple(ux, uy, uz, ut, ioffs)

    ampl, phase = xy_to_ampl_phase(ux, uy)

    uevo = zeros(ComplexF64, size(ampl)[1], 4,4)

    for idig in eachindex(ampl)     # points of measurement
        phasesin, phasecos = sincosd(phase[idig])
        hpulse = 2*pi*ampl[idig]*(phasecos*ix[:,:,2]+phasesin*iy[:,:,2])
        hcs2 = -2*pi*(ioffs-uz[idig])*iz[:,:,2]    # chemical shift spin 2 including offset
        uevo[idig,:,:] = exp(-i*(hcs2+hpulse)*ut[idig])
    end

    phip, phim = 0.0, 0.0
    Ap, Bp, Cp, Dp, Am, Bm, Cm, Dm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    fuevop, fuevom = [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]
    
    for idig in reverse(eachindex(ampl))
        Dp = real(uevo[idig,1,1])
        Cp = imag(uevo[idig,1,1])
        Bp = real(uevo[idig,1,2])
        Ap = imag(uevo[idig,1,2])
        Dm = real(uevo[idig,3,3])
        Cm = imag(uevo[idig,3,3])
        Bm = real(uevo[idig,3,4])
        Am = imag(uevo[idig,3,4])
        if idig == size(ampl)[1]
            fuevop = [Ap, Bp, Cp, Dp]
            fuevom = [Am, Bm, Cm, Dm]
        else
            fuevop = [[Dp, Cp, -Bp, -Ap] [-Cp, Dp, Ap, -Bp] [Bp, -Ap, Dp, -Cp] [Ap, Bp, Cp, Dp]] * fuevop
            fuevom = [[Dm, Cm, -Bm, -Am] [-Cm, Dm, Am, -Bm] [Bm, -Am, Dm, -Cm] [Am, Bm, Cm, Dm]] * fuevom
        end
    end

    if fuevop[4] <= 1
        phip = 2*acos(fuevop[4])
    end
    if fuevom[4] <= 1
        phim = 2*acos(fuevom[4])
    end

    return fuevop, fuevom
end


"simple Waugh criterion is simulated\n
output: lambda"
function simulate_waugh_simple(ux, uy, uz, ut, eps, noffs, offs)

    println("starting simulation")

    phipall = zeros(noffs)
    phimall = zeros(noffs)
    
    tcyc = sum(ut)
    lambda = zeros(noffs)

    for ioffs in eachindex(offs)

        fuevop, fuevom = get_waugh_rotation_simple(ux, uy, uz, ut, offs[ioffs])
        fuevopeps, fuevomeps = get_waugh_rotation_simple(ux, uy, uz, ut, offs[ioffs]+eps)

        phip = 2*asin(sqrt(fuevop[1]^2+fuevop[2]^2+fuevop[3]^2))
        phipeps = 2*asin(sqrt(fuevopeps[1]^2+fuevopeps[2]^2+fuevopeps[3]^2))
        phim = 2*asin(sqrt(fuevom[1]^2+fuevom[2]^2+fuevom[3]^2))
        phimeps = 2*asin(sqrt(fuevomeps[1]^2+fuevomeps[2]^2+fuevomeps[3]^2))

        phipall[ioffs] = phip
        phimall[ioffs] = phim

        lambda[ioffs] = (phip-phipeps)/(tcyc*eps)
        
    end

    return lambda
end


"write text to file \"filename\"\n
output: written in file"
function write_to_file(filename, text)

    io = open(filename,"w")
    println(io, text) 
    close(io)

    return
end


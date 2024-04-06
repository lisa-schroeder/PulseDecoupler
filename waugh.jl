using QuadGK
using StatsBase
# using DoubleFloats

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
        @show ioffs

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



"calculate rotation angle and effective B-field \n
phip: rotation angle for offset + J/2\n
phim: rotation angle for offset - J/2\n
beffp: effective B-field for offset + J/2\n
beffm: effective B-field for offset - J/2\n
uevop: evolution propagator for offset + J/2\n
uevom: evolution propagator for offset - J/2\n
output: phip, phim, beffp, beffm, uevop, uevom"
function get_waugh_rotation_exact(ifJ, J, ux, uy, uz, ut, ioffs)

    # J=145
    ampl, phase = xy_to_ampl_phase(ux, uy)

    uevo = zeros(ComplexF64, size(ampl)[1], 4,4)
    beffp = zeros(3)
    beffm = zeros(3)
    toobig = 0

    if ifJ == true
        hj12 = 2*pi*J*(iz[:,:,1]*iz[:,:,2])           # Kopplungsentwicklung zwischen Spin 1 und 2
    else
        hj12 = 0*(iz[:,:,1]*iz[:,:,2])
    end

    for idig in eachindex(ampl)     # points of measurement
        phasesin, phasecos = sincosd(phase[idig])
        hpulse = 2*pi*ampl[idig]*(phasecos*ix[:,:,2]+phasesin*iy[:,:,2])
        hcs2 = -2*pi*(ioffs-uz[idig])*iz[:,:,2]    # chemical shift spin 2 including offset
        uevo[idig,:,:] = exp(-i*(hj12+hcs2+hpulse)*ut[idig])
    end

    phip, phim = 0.0, 0.0
    Ap, Bp, Cp, Dp, Am, Bm, Cm, Dm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    fuevop, fuevom = [0.0, 0.0, 0., 0.0], [0.0, 0.0, 0., 0.0]
    
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
    else
        print("a plus too big at offset ", ioffs, " Hz: ", fuevop[4])
        # print(ioffs)
        # print(" Hz: ")
        # println(fuevop[4])
        toobig =1
    end
    if fuevom[4] <= 1
        phim = 2*acos(fuevom[4])
    else
        print("a minus too big at offset ", ioffs, " Hz: ", fuevom[4])
        # print(ioffs)
        # print(" Hz: ")
        # println(fuevom[4])
        toobig = 1
    end

    # @show phip, phim
    sinp = sin(phip/2)
    sinm = sin(phim/2)
    beffp[1] = fuevop[1]/sinp
    beffm[1] = fuevom[1]/sinm
    beffp[2] = fuevop[2]/sinp
    beffm[2] = fuevom[2]/sinm
    beffp[3] = fuevop[3]/sinp
    beffm[3] = fuevom[3]/sinm

    return phip, phim, beffp, beffm, fuevop, fuevom
end



"exact Waugh criterion is simulated\n
output: J1, J2, I1, I2"
function simulate_waugh_exact(ifJ, ux, uy, uz, ut, noffs, offs, J)

    println("starting simulation")

    phipall = zeros(noffs)
    phimall = zeros(noffs)
    beffpall = zeros(noffs,3)
    beffmall = zeros(noffs,3)


    tcyc = sum(ut)
    J1 = zeros(noffs)
    J2 = zeros(noffs)
    I1 = zeros(noffs)
    I2 = zeros(noffs)
    toobig = zeros(noffs)

    for ioffs in eachindex(offs)

        phip, phim, beffp, beffm, fuevop, fuevom = get_waugh_rotation_exact(ifJ, J, ux, uy, uz, ut, offs[ioffs])
        toobig[ioffs] = toobig_temp

        phipall[ioffs] = phip
        phimall[ioffs] = phim

        beffpall[ioffs,:] = beffp
        beffmall[ioffs,:] = beffm

        # Jeff[ioffs] = (phip-phim)/tcyc
        # lambda[ioffs] = Jeff[ioffs]/J
        J1[ioffs] = (phip-phim)/(2*pi*tcyc)
        J2[ioffs] = (phip+phim)/(2*pi*tcyc)
        I1[ioffs] = (1+dot(beffp,beffm))/2
        I2[ioffs] = (1-dot(beffp,beffm))/2
        # @show offs[ioffs], Jeff[ioffs], J1[ioffs], J2[ioffs], I1[ioffs], I2[ioffs], lambda[ioffs]
        # @show Jeff[ioffs], lambda[ioffs]
    end


    return J1, J2, I1, I2
end


"write text to file \"filename\"\n
output: written in file"
function write_to_file(filename, text)

    io = open(filename,"w")
    println(io, text) 
    close(io)

    return
end


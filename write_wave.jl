
"write ampl, phase and info to file according to Bruker file format\n
output: written in file"
function write_bruker_wave(filename,ampl,phase,info)

    io = open(filename,"w")
    println(io,"##TITLE=") 
    println(io,"##JCAMP-DX=") 
    println(io,"##DATA TYPE= Shape Data")
    println(io,"##ORIGIN= ")
    println(io,"##OWNER= ")
    println(io,"##DATE= ")
    println(io,"##TIME= ")
    println(io,"##\$SAPE_PARAMETER=")
    println(io,"##MINX=")
    println(io,"##MAXX=")
    println(io,"##MINY=")
    println(io,"##MAXY=")
    println(io,"##\$SHAPE_EXMODE=")
    println(io,"##\$SHAPE_TOTROT=")
    println(io,"##\$SHAPE_BWFAC=")
    println(io,"##\$SHAPE_INTEGFAC=")
    println(io,"##\$SHAPE_MODE=")
    for count = 1:length(info)
        println(io,"## ",info[count])
    end
    println(io,"##NPOINTS= ",length(ampl))         
    println(io,"##XYPOINTS=")
    ma = maximum(ampl)
    for count = 1:length(ampl)
        println(io,"  ",lpad(ampl[count]*100/ma,20),",  ",lpad(phase[count],20))   
    end
    write(io,"##END= ")
    close(io)
    
    
    return
end
    

"save pulse components ux, uy and ut as xyt file\n
output: written in file"
function write_xyt_wave(filename, info, ux, uy, ut)
  
  io = open(filename,"w")
  @printf(io, "%s\n", info) 
  for dig in eachindex(ux)
      @printf(io, "%20.13f   %20.13f   %20.13e\n", ux[dig], uy[dig], ut[dig])   
  end
  close(io)

end

"save pulse components ux, uy, offset and ut as xyzt file\n
output: written in file"
function write_xyzt_wave(filename, info, ux, uy, uz, ut)
  
  io = open(filename,"w")
  @printf(io, "%s\n", info) 
  for dig in eachindex(ux)
      @printf(io, "%20.13f   %20.13f   %20.13f   %20.13e\n", ux[dig], uy[dig], uz[dig], ut[dig])   
  end
  close(io)

end

"save pulse components ux, uy, offset and ut as xygt file
grad: strength of gradient in T/M of each time step (not graduz!!!)\n
output: written in file"
function write_xyzgt_wave(filename, info, ux, uy, uz, grad, ut, acqux, acquy, acquz, acqgrad, acqut, homux, homuy, homuz, homgrad, homut)
  
    io = open(filename,"w")
    @printf(io, "%s\n", info) 
    @printf(io, "#    x                      y                      z                      g                  t                  \n")

    if ux != []
        @printf(io, "SEQ=start\n")
        for dig in eachindex(ux)
            @printf(io, "%20.13f   %20.13f   %20.13f   %20.13f   %20.13e\n", ux[dig], uy[dig], uz[dig], grad[dig], ut[dig])   
        end
        @printf(io, "SEQ=end\n")
    end

    if acqux != []
        @printf(io, "ACQ=start\n")
        for dig in eachindex(acqux)
            @printf(io, "%20.13f   %20.13f   %20.13f   %20.13f   %20.13e\n", acqux[dig], acquy[dig], acquz[dig], acqgrad[dig], acqut[dig])   
        end
        @printf(io, "ACQ=end\n")
    end

    if homux != []
        @printf(io, "HOMOD=start\n")
        for dig in eachindex(homux)
            @printf(io, "%20.13f   %20.13f   %20.13f   %20.13f   %20.13e\n", homux[dig], homuy[dig], homuz[dig], homgrad[dig], homut[dig])   
        end
        @printf(io, "HOMOD=end\n")
    end
  close(io)

end

"converts bruker file format to xygt file format\n
so far only either seq, acq oder homod is implemented in the converter\n
tpulse in ms\n
rfmax in Hz\n
info is a string, containing additional information for the reader\n
output: written in file"
function convert_bruker_to_xyzgt(brukerfilename, xygtfilename, rfmax, tpulse, seq, acq, homod, info)
    ampl, phase = read_bruker_wave(brukerfilename)
    ux, uy = ampl_phase_to_xy(ampl, phase)
    ux = ux.*rfmax/100
    uy = uy.*rfmax/100
    uz = zeros(size(ux)[1])
    grad = zeros(size(ux)[1])
    ut = ones(size(ux)[1]).*(tpulse/1000/size(ux)[1])
    if seq
        write_xyzgt_wave(xygtfilename, "# " * info * "\n# rfmax = " * string(rfmax) * " Hz, tpulse = " * string(tpulse) * " ms", ux, uy, uz, grad, ut, [], [], [], [], [], [], [], [], [], [])
    elseif acq
        write_xyzgt_wave(xygtfilename, "# " * info * "\n# rfmax = " * string(rfmax) * " Hz, tpulse = " * string(tpulse) * " ms", [], [], [], [], [], ux, uy, uz, grad, ut, [], [], [], [], [])
    elseif homod
        write_xyzgt_wave(xygtfilename, "# " * info * "\n# rfmax = " * string(rfmax) * " Hz, tpulse = " * string(tpulse) * " ms", [], [], [], [], [], [], [], [], [], [], ux, uy, uz, grad, ut)
    end
end


# convert_bruker_to_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/cw.bruker", "C:/Users/Lisa/Documents/KIT/WS202324/Pulse/cw-acq.xyzgt", 10000, 0.05, false, true, false, "cw")
# convert_bruker_to_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/BUSS.bruker", "C:/Users/Lisa/Documents/KIT/WS202324/Pulse/BUSS.xyzgt", 14000, 123.2, false, true, false, "BUSS")
# convert_bruker_to_xyzgt("C:/Users/Lisa/Documents/KIT/WS202324/Pulse/Rsnob_Lisa.bruker", "C:/Users/Lisa/Documents/KIT/WS202324/Pulse/Rsnob_shifted.xyzgt", 234, 10, false, false, true, "Rsnob shifted by 2250 Hz")
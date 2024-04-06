"method to get the greatest commond divisor a for two floats\n
output: a"
function float_gcd(a, b, rtol = 1e-05, atol = 1e-08)
    t = min(abs(a), abs(b))
    while abs(b) > rtol * t + atol
        a, b = b, a % b
    end
    return a
end
     

"method to generate composite pulses, information is defined in filename_info: \n
angle and phase defined in each line, different composite pulses are seperated by #\n
output: written in file"
function generate_all_composite(filename_info, amplitude)

    open(filename_info, "r")
    lines=readlines(filename_info)
    
    for count = 1:size(lines)[1] 
        if occursin("#", lines[count]) && count != size(lines)[1] 
            durations=empty([1.0], Float64)
            phase=empty([1.0], Float64)    
            angle=empty([1.0], Float64)    
            count_2 = count + 1
            while !occursin("#", lines[count_2])
                a,b=map(x-> parse(Float64,x), collect(eachsplit(lines[count_2], ",")))
                append!(angle, a)
                a = 1/amplitude/360*a*1000000
                append!(durations, a)
                append!(phase, b)
                count_2 = count_2 + 1
            end
            
            total_durations = sum(durations)

            # info describes which pulses make up the composite pulse
            info = ""
            for n in eachindex(phase)
                if n == 1
                    info = info * string(Int(round(angle[n]))) * "(" * string(Int(round(phase[n]))) * ")"
                else
                    info = info * "-" * string(Int(round(angle[n]))) * "(" * string(Int(round(phase[n]))) * ")"
                end
                
            end

            #x=[cos(ang) for ang in angle/180*pi]
            #y=[sin(ang) for ang in angle/180*pi]
            gcd_1 = 0.0
            if size(durations)[1] > 1
                gcd_1 = float_gcd(durations[1], durations[2])
            else
                gcd_1 = durations[1]
            end

            for i in 1:size(durations)[1]
                gcd_1 = float_gcd(gcd_1, durations[i])
            end
            
            # phase_2 is required that the number of data points corresponds to the durations of the pulses
            point = [dur/gcd_1 for dur in durations]
            println(point)
            phase_2=empty([1.0], Float64)
            for len in eachindex(durations)
                for i in 1:point[len]*100
                    append!(phase_2, phase[len])
                end
            end

            ampl = ones(1,size(phase_2)[1])
            #@show(phase_2, point, durations, angle, phase, total_durations, [info])
            #@show(ampl, phase_2)
            
            # generate the wave file in bruker format
            write_bruker_wave(info, ampl, phase_2, [info, "durations = " * string(total_durations) * " us", "RF" * string(amplitude/1000)])
            
        end
        
    end
    
end
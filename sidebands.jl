
using CSV
using DataFrames
using StatsBase
using Printf


"reading in intensity of sidebands or center peak of the csv file
Js: all coupling constants J in Hz\n
sidebands: intensity as written in the csv file\n
output: Js, sidebands"
function read_in_biggest_sidebands(filename)
    data = DataFrame(CSV.File(filename))
    Js = zeros(size(data)[1])
    sidebands = zeros(size(data)[1],size(data)[2]-2)

    for line in 1:size(data)[1]
        @show line
        Js[line] = data[line,1]
        for row = 1:size(data)[2]-2
            @show row
            sidebands[line,row] = data[line,row+1]
        end
    end

    return Js, sidebands
end


"plotting of intensity of sidebands\n
fitting to the function pJ^2\n
sidebands: intensity of largest sideband for each J and each offset\n
xlimits, ylimits: limits of plot\n
p0: estimated fitting parameter\n
output: -"
function plot_fit_sidebands(Js, sidebands, xlimits, ylimits, p0, noffs)

    iJ = size(sidebands)[1]
    scatter(Js, sidebands, linewidth=2, color=palette([:darkblue, :red, :darkblue], noffs), legend=false, xlabel="coupling constant "*L"J_{IS}", dpi=500, ylabel="intensity of largest sideband (%)")

    @. modelside(x, p) = p[1]*x^2
    for ioffs = 1:noffs
        fitside = LsqFit.curve_fit(modelside, Js[1:iJ], sidebands[1:iJ,ioffs], p0)
        fitcurve = [coef(fitside)[1]*x^2 for x in Js[1:iJ]]
        # SStot = sum((sidebands[1:iJ,ioffs].-mean(sidebands[1:iJ,ioffs])).^2)
        # SSres = sum((sidebands[1:iJ,ioffs].-fitcurve).^2)
        # Rsq = 1-SSres/SStot
        plot!(Js[1:iJ], fitcurve, xlims=xlimits, ylims=ylimits, dpi=500, color=palette([:darkblue, :red, :darkblue], noffs)[ioffs], linewidth=2)
        # @show coef(fitside)
    end

    display(plot!())

end


"plotting of intensity of center peak\n
fitting to the function 100-pJ^2\n
maxvals: intensity of center peak for each J and each offset\n
xlimits, ylimits: limits of plot\n
p0: estimated fitting parameter\n
output: -"
function plot_fit_maxintens(Js, maxvals, xlimits, ylimits, p0, noffs)

    iJ = size(maxvals)[1]
    @show iJ

    scatter(Js, maxvals, dpi=500, linewidth=2, legend=false,  color=palette([:darkblue, :red, :darkblue], noffs), xlabel="coupling constant "*L"J_{IS}", ylabel="maximum intensity (%)")

    @. modelmax(x, p) = 100+p[1]*x^2
    for ioffs = 1:noffs
        fitmax = LsqFit.curve_fit(modelmax, Js[1:iJ], maxvals[1:iJ,ioffs], p0)
        fitcurve = [100+coef(fitmax)[1]*x^2 for x in Js[1:iJ]]
        # SStot = sum((maxvals[1:iJ,ioffs].-mean(maxvals[1:iJ,ioffs])).^2)
        # SSres = sum((maxvals[1:iJ,ioffs].-fitcurve).^2)
        # Rsq = 1-SSres/SStot
        plot!(Js[1:iJ], fitcurve, xlims=xlimits, ylims=ylimits, color=palette([:darkblue, :red, :darkblue], noffs)[ioffs], dpi=500, linewidth=2)
        # @show coef(fitmax)
    end

    display(plot!())

end


"function to write the maximum intensities of the simulation at different Js into a csv file\n
!!! intensity is divided by 2.505 (only for a linewidth of 6 Hz, dwell = 0.2 ms and npoints = 5000) !!!
maxintens: maximum intensity of center peak\n
maxintensside: intensity of largest sideband\n
filename: for saving the csv file\n
output: saved to file"
function write_maxintensities_to_csv(maxintens, maxintensside, noffs, filename)
    textside = "At J,offs1,offs2,offs3,offs4,offs5,offs6,offs7,offs8,offs9,offs10,offs11,offs12,offs13,offs14,offs15,offs16,offs17,offs18,offs19,offs20,offs21,offs22,offs23,offs24,offs25,offs26,offs27,offs28,offs29,offs30,offs31\n"
    textmax = "At J,offs1,offs2,offs3,offs4,offs5,offs6,offs7,offs8,offs9,offs10,offs11,offs12,offs13,offs14,offs15,offs16,offs17,offs18,offs19,offs20,offs21,offs22,offs23,offs24,offs25,offs26,offs27,offs28,offs29,offs30,offs31\n"

    for iJ = 1:size(maxintens[1,:])[1]
        textside = textside * string(allJs[iJ])
        textmax = textmax * string(allJs[iJ])
        for ioffs = 1:noffs
            textside = textside * string(",", maxintensside[ioffs,iJ]/2.505,)
            textmax = textmax * string(",", maxintens[ioffs,iJ]/2.505,)
        end
        textside = textside * "\n"
        textmax = textmax * "\n"
    end


    write_to_file(filename*"_sidebands.csv", textside)
    write_to_file(filename*"_maxintens.csv", textmax)

end



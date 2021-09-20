using DataFrames: DataAPI
using ArgParse
function parse_cmd()
    cfg = ArgParseSettings()
    @add_arg_table cfg begin
        "--annotate", "-a"
        help = "add annotation in the figure"
        action = :store_true
        "--pdf", "-p"
        help = "generate figure in PDF format"
        action = :store_true
    end
    return parse_args(cfg)
end


using FilePaths
using Glob
using Printf
function get_file_list(series::String, lev::Integer, ws::Integer, url::Integer)
    list = glob(string("bench/", series, "_lev", lev, "_ws", ws, "_url", Printf.@sprintf("%04d", url), "_*.time*.mean.dat"))
    # println(string("bench/", series, "_lev", lev, "_ws", ws, "_url", sprintf("%04d", url), "_*.time*.mean.dat"))
    # println(list)
    # println(length(list))
    return list
end

using DataFrames
using CSV
function read_tsv(file::String)
    tsv = DataFrame(CSV.File(file))
    calcGrav_dev = tsv[!, :"#calcGrav_dev"]
    # println(calcGrav_dev)
    # println(calcGrav_dev[begin])
    return calcGrav_dev[begin]
end


using PyPlot
include("../util/pyplot.jl")


using Statistics
function main()
    # read options
    argv = parse_cmd()
    output_pdf = argv["pdf"]
    add_memo = argv["annotate"]

    # find the latest series of simulation results
    series = "m31"
    l2lev = [0, 4, 5, 6, 7, 8]
    wsopt = [0, 1]
    Nurl = [1, 2, 4, 8, 16, 32, 64, 128]

    # initialize matplotlib
    util_pyplot.config()

    # show dependence on L2lev
    nx = length(Nurl)
    ny = length(wsopt)
    fig = util_pyplot.set_Panel(nx = nx, ny = ny)

    style_cap = Dict("color" => util_pyplot.call(fig.color, id = 0))
    style_box = Dict("color" => util_pyplot.call(fig.color, id = 0), "lw" => 0.75 * fig.lw)
    style_whi = Dict("color" => util_pyplot.call(fig.color, id = 0), "markersize" => fig.ms)
    style_med = Dict("color" => util_pyplot.call(fig.color, id = 1), "lw" => 0.5 * fig.lw)

    # generate summary CSV file
    summary = DataFrame(L2lev = [], shfl = [], unroll = [], minimum = [], median = [])

    tmin = typemax(Float64)
    tmax = typemin(Float64)
    for jj in 1:length(wsopt)
        ws = wsopt[jj]
        for ii in 1:length(Nurl)
            url = Nurl[ii]

            # allocate array for print out the measured scores
            measured_median  = Array{Float64, 1}(undef, length(l2lev))
            # measured_average = Array{Float64, 1}(undef, length(l2lev))
            measured_fastest = Array{Float64, 1}(undef, length(l2lev))
            measured_slowest = Array{Float64, 1}(undef, length(l2lev))
            # read measured results and plot the data
            for kk in 1:length(l2lev)
                lev = l2lev[kk]
                list = get_file_list(series, lev, ws, url)

                # initialize array for results
                elapse = zeros(Float32, length(list))
                # read measured time in for loop
                for ii in 1:length(list)
                    elapse[ii] = read_tsv(list[ii])
                end

                # obtain statistics (minimum, median, average, maximum, #entries)
                fastest = minimum(elapse)
                slowest = maximum(elapse)
                median  = Statistics.median(elapse)
                # average = Statistics.mean(elapse)
                # println("L2 = ", lev, ", WS = ", ws, ", Nunroll = ", url
                # , ": minimum = ", fastest
                # , ", median = ", median
                # , ", mean = ", average
                # , ", maximum = ", slowest
                # , "; N_measured = ", length(list))
                tmin = minimum([tmin, fastest])
                tmax = maximum([tmax, slowest])
                measured_median[kk]  = median
                # measured_average[kk] = average
                measured_fastest[kk] = fastest
                measured_slowest[kk] = slowest
                push!(summary, (lev, ws, url, fastest, median))

                fig.ax[ii, jj].boxplot(elapse, positions = [kk], labels = [lev], widths = [0.5]
                , capprops = style_cap
                , boxprops = style_box
                , whiskerprops = style_whi
                , medianprops = style_med
                , flierprops = Dict("color" => util_pyplot.call(fig.color, id = 0), "marker" => "x", "markersize" => 0.5 * fig.ms)
                )
            end

            # add captions
            caption = string("(", Char(97 + (ii - 1) + fig.nx * (fig.ny - jj)), ")")
            wslab = ws == 0 ? "~w/o shfl" : "~w/ shfl"
            at = fig.ax[ii, jj]
            at.text(0.03, 0.97, string(caption, wslab, L", $N_\mathrm{unroll} =$~", url), color = "black", fontsize = fig.fs, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))
            if add_memo
                for kk in 1:length(l2lev)
                    ypos = 0.95 - kk * 0.04
                    caption = string(L"$N_\mathrm{lev} =$~", l2lev[kk]
                    , L": $t_\mathrm{min} =$~", Printf.@sprintf("%.2f", measured_fastest[kk] * 1.0e+3), L"~$\si{ms}$"
                    , L", $t_\mathrm{med} =$~", Printf.@sprintf("%.2f", measured_median[kk] * 1.0e+3), L"~$\si{ms}$"
                    , L", $t_\mathrm{max} =$~", Printf.@sprintf("%.2f", measured_slowest[kk] * 1.0e+3), L"~$\si{ms}$")
                    # , L", $t_\mathrm{ave} =$~", Printf.@sprintf("%.2f", measured_average[kk] * 1.0e+3), L"~$\si{ms}$")
                    at.text(0.08, ypos, caption, color = "black", fontsize = fig.fs * 0.8, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))
                end
            end
        end
    end

    for at in fig.ax
        at.grid()
        at.set_ylim(util_pyplot.scale_axis(tmin, tmax, logPlt = false))
    end
    for jj in 1:ny
        fig.ax[begin, jj].set_ylabel(L"$\overline{t_\mathrm{walk}}~\bqty{\si{s}}$", fontsize = fig.fs)
    end
    for ii in 1:nx
        fig.ax[ii, begin].set_xlabel(L"$N_\mathrm{lev}$", fontsize = fig.fs)
    end

    # save figures
    tag = add_memo ? "_compare_note" : "_compare"
    fig.fig.savefig(string("fig/", series, tag, ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, tag, ".pdf"), format = "pdf", bbox_inches = "tight")
    end

    # write summary CSV file
    summary |> CSV.write(string(series, "_measured.csv"), delim = ',', writeheader = true)

	fig = nothing
    PyPlot.close("all")
    return nothing
end


main()

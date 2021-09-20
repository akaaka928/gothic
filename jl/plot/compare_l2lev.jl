using ArgParse
function parse_cmd()
    cfg = ArgParseSettings()
    @add_arg_table cfg begin
        "--pdf", "-p"
        help = "generate figure in PDF format"
        action = :store_true
    end
    return parse_args(cfg)
end


using FilePaths
using Glob
function get_file_list(root::String, series::String, lev::Integer)
    list = glob(string(root, "bench/", series, "_lev", lev, "_*.time*.mean.dat"))
    return list
end

using DataFrames
using CSV
function read_tsv(file::String)
    tsv = DataFrame(CSV.File(file))
    calcGrav_dev = tsv[!, :"#calcGrav_dev"]
    return calcGrav_dev[begin]
end

using PyPlot
include("../util/pyplot.jl")


function main()
    # read options
    argv = parse_cmd()
    output_pdf = argv["pdf"]

    # find the latest series of simulation results
    root_async = "../l2dep/"
    root_woasync = "../l2dep_woasync/"
    series = "m31"
    # l2lev = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    l2lev = [0, 1, 2, 3, 4, 5, 6, 7, 8]

    # initialize matplotlib
    util_pyplot.config()

    # show dependence on L2lev
    fig = util_pyplot.set_Panel(nx = 2)
    box = util_pyplot.set_Panel(nx = 2)
    # for lev in [0]
    ymin = typemax(Float64)
    ymax = typemin(Float64)
    for lev in l2lev
        # generate file list
        list_async = get_file_list(root_async, series, lev)
        list_woasync = get_file_list(root_woasync, series, lev)

        # initialize array for results
        elapse_async = zeros(Float32, length(list_async))
        elapse_woasync = zeros(Float32, length(list_woasync))
        # read measured time in for loop
        for ii in 1:length(list_async)
            elapse_async[ii] = read_tsv(list_async[ii])
        end
        for ii in 1:length(list_woasync)
            elapse_woasync[ii] = read_tsv(list_woasync[ii])
        end
        # println(elapse)
        # min = minimum([minimum(elapse_async), minimum(elapse_woasync)])
        # max = maximum([maximum(elapse_async), maximum(elapse_woasync)])
        # ymin = min(minimum(minimum[elapse_async]), minimum(minimum[elapse_woasync]))
        # println(minimum(elapse_async))
        # ymin = minimum(minimum(elapse_async))
        # println("next")
        # println(ymin)
        # ymax = max(maximum(maximum(elapse_async)), maximum(maximum(elapse_woasync)))
        # println(ymax)
        ymin = minimum([ymin, minimum(elapse_async), minimum(elapse_woasync)])
        ymax = maximum([ymax, maximum(elapse_async), maximum(elapse_woasync)])

        # show box-and-whisker plot
        style_cap = Dict("color" => util_pyplot.call(fig.color, id = 0))
        style_box = Dict("color" => util_pyplot.call(fig.color, id = 0), "lw" => 0.75 * fig.lw)
        style_whi = Dict("color" => util_pyplot.call(fig.color, id = 0), "markersize" => fig.ms)
        style_med = Dict("color" => util_pyplot.call(fig.color, id = 1), "lw" => 0.5 * fig.lw)
        fig.ax[1].boxplot(elapse_async, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , flierprops = Dict("color" => util_pyplot.call(fig.color, id = 0), "marker" => "x", "markersize" => 0.5 * fig.ms)
        )
        fig.ax[2].boxplot(elapse_woasync, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , flierprops = Dict("color" => util_pyplot.call(fig.color, id = 0), "marker" => "x", "markersize" => 0.5 * fig.ms)
        )
        box.ax[1].boxplot(elapse_async, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , showfliers = false
        )
        box.ax[2].boxplot(elapse_woasync, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , showfliers = false
        )
        elapse_async = nothing
        elapse_woasync = nothing
    end
    fig.ax[begin].set_ylabel(L"$\overline{t_\mathrm{walk}}~\bqty{\si{s}}$", fontsize = fig.fs)
    box.ax[begin].set_ylabel(L"$\overline{t_\mathrm{walk}}~\bqty{\si{s}}$", fontsize = box.fs)
    for ii in 1:2
        fig.ax[ii].set_xlabel(L"$N_\mathrm{lev}$", fontsize = fig.fs)
        box.ax[ii].set_xlabel(L"$N_\mathrm{lev}$", fontsize = box.fs)
        fig.ax[ii].grid()
        box.ax[ii].grid()
        fig.ax[ii].set_ylim(util_pyplot.scale_axis(ymin, ymax, logPlt = false))
        box.ax[ii].set_ylim(util_pyplot.scale_axis(ymin, ymax, logPlt = false))
    end

    # add caption
    caption = string("(", Char(97 + (1 - 1)), ")")
    at = fig.ax[1]
    at.text(0.03, 0.97, string(caption, "~w/~", L"\texttt{memcpy\_async()}"), color = "black", fontsize = fig.fs, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))
    at = box.ax[1]
    at.text(0.03, 0.97, string(caption, "~w/~", L"\texttt{memcpy\_async()}"), color = "black", fontsize = box.fs, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))
    caption = string("(", Char(97 + (2 - 1)), ")")
    at = fig.ax[2]
    at.text(0.03, 0.97, string(caption, "~w/o~", L"\texttt{memcpy\_async()}"), color = "black", fontsize = fig.fs, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))
    at = box.ax[2]
    at.text(0.03, 0.97, string(caption, "~w/o~", L"\texttt{memcpy\_async()}"), color = "black", fontsize = box.fs, horizontalalignment = "left", verticalalignment = "top", transform = at.transAxes, bbox = Dict("facecolor" => "white", "edgecolor" => "None", "alpha" => 0.75))

    # save figures
    fig.fig.savefig(string("fig/", series, "_l2dep_compare", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    box.fig.savefig(string("fig/", series, "_l2box_compare", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, "_l2dep_compare", ".pdf"), format = "pdf", bbox_inches = "tight")
        box.fig.savefig(string("fig/", series, "_l2box_compare", ".pdf"), format = "pdf", bbox_inches = "tight")
    end

	fig = nothing
    PyPlot.close("all")
    return nothing
end


main()

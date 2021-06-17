using ArgParse
function parse_cmd()
    cfg = ArgParseSettings()
    @add_arg_table cfg begin
        # "--target", "-t"
        # help = "filename of the target files"
        # arg_type = String
        # default = ""
        "--pdf", "-p"
        help = "generate figure in PDF format"
        action = :store_true
    end
    return parse_args(cfg)
end

using FilePaths
using Glob
function get_file_list(folder::String, series::String, tag::String)
    list = glob(string(folder, series, "*acceleration.block.", tag, ".8388608.time.log"))
    # println(list)
    # println(length(list))
    return list
end

using DataFrames
using CSV
function read_tsv(file::String)
    tsv = CSV.read(file, DataFrame; delim = "\t", header = ["dacc", "sec", "step", "sec_per_step"], comment = "#")
    return tsv
end

# using StatsBase
# function read_dacc(whole::Vector{Float64})
#     # println(typeof(whole))
#     tmp = StatsBase.countmap(whole)
#     println(tmp)
#     return [key for (key, val) in tmp if val == 1]
# end


function get_measured_results(series::String, path::String, gen::String)
    # read measured results
    list = get_file_list(path, series, gen)
    result = read_tsv(list[begin])
    for ii in 2:length(list)
        result = vcat(result, read_tsv(list[ii]))
    end
    # println(result)
    # println(result[!, "dacc"])
    # println(result[result.dacc .== 0.5, :])
    # println(result[result.dacc .== 0.5, :].sec_per_step)
    dacc = sort(unique(result.dacc), rev = true)
    # println(dacc)

    time = Array{Float64, 1}(undef, length(dacc))
    for ii in 1:length(dacc)
        data = result[result.dacc .== dacc[ii], :].sec_per_step
        time[ii] = minimum(data)
    end
    # println(time)
    return dacc, time
end


using PyPlot
include("../util/pyplot.jl")


function main()
    # read options
    argv = parse_cmd()
    output_pdf = argv["pdf"]

    # find the latest series of simulation results
    # series = argv["target"] == "" ? util_find.latest_series() : argv["target"]
    series = "m31"

    # initialize matplotlib
    util_pyplot.config()
    fig = util_pyplot.set_Panel()

    p100_dacc, p100_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc60")
    # v100v_dacc, v100v_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc70")
    v100p_dacc, v100p_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc70_60")
    a100p_dacc, a100p_time = get_measured_results(series, "../../../a100/201225gothic/adep/1410/log/", "cc80")
    a100s_dacc, a100s_time = get_measured_results(series, "log/", "cc80")
    # println(a100_dacc)
    # println(a100_time)
    fig.ax[begin].plot(p100_dacc, p100_time, util_pyplot.call(fig.point, id = 2), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 3), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 0), label = "P100 (SXM)")
    # fig.ax[begin].plot(v100v_dacc, v100v_time, util_pyplot.call(fig.point, id = 1), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 1), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 2), label = L"V100 (SXM2, compute\_70)", markerfacecolor = "none")
    fig.ax[begin].plot(v100p_dacc, v100p_time, util_pyplot.call(fig.point, id = 1), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 2), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 2), label = "V100 (SXM2)")
    fig.ax[begin].plot(a100p_dacc, a100p_time, util_pyplot.call(fig.point, id = 0), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 1), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 1), label = "A100 (PCIe)", markerfacecolor = "none")
    fig.ax[begin].plot(a100s_dacc, a100s_time, util_pyplot.call(fig.point, id = 0), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 0), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 1), label = "A100 (SXM4)")

    fig.ax[begin].set_xlabel(L"$\varDelta_\mathrm{acc}$", fontsize = fig.fs)
    fig.ax[begin].set_ylabel(L"$t_\mathrm{step}~\bqty{\si{s}}$", fontsize = fig.fs)
    fig.ax[begin].grid()
    fig.ax[begin].loglog()

    # add legend
    handles, labels = fig.ax[begin].get_legend_handles_labels()
	fig.ax[begin].legend(handles, labels, numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig.fs)

    # save figures
    fig.fig.savefig(string("fig/", series, "_adep", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, "_adep", ".pdf"), format = "pdf", bbox_inches = "tight")
    end

	fig = nothing
    PyPlot.close("all")
    return nothing

end


main()

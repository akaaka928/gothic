using ArgParse
function parse_cmd()
    cfg = ArgParseSettings()
    @add_arg_table cfg begin
        # "--target", "-t"
        # help = "filename of the target files"
        # arg_type = String
        # default = ""
        "--speedup", "-s"
        help = "show speed-up ratio"
        action = :store_true
        "--pdf", "-p"
        help = "generate figure in PDF format"
        action = :store_true
    end
    return parse_args(cfg)
end

using FilePaths
using Glob
function get_file_list(folder::String, series::String, tag::String)
    # list = glob(string(folder, series, "*acceleration.block.", tag, ".8388608.time.log"))
    list = glob(string(folder, series, "*acceleration.block.", tag, ".*.time.log"))
    # println(list)
    # println(glob(string(folder, "*acceleration*")))
    # println(length(list))
    return list
end

using DataFrames
using CSV
function read_tsv(file::String)
    tsv = CSV.read(file, DataFrame; delim = "\t", header = ["dacc", "sec", "step", "sec_per_step"], comment = "#")
    return tsv
end

function read_csv(file::String)
    csv = CSV.read(file, DataFrame; delim = ",", comment = "#")
    # println(csv)
    # println(csv.Ntot)
    # println(csv.elapsed_step)
    return csv.Ntot, csv.elapsed_step
end


function get_measured_results(series::String, path::String, gen::String)
    # read measured results
    list = get_file_list(path, series, gen)
    result = read_tsv(list[begin])
    for ii in 2:length(list)
        result = vcat(result, read_tsv(list[ii]))
    end
    # # println(result)
    # # println(result[!, "dacc"])
    # # println(result[result.dacc .== 0.5, :])
    # # println(result[result.dacc .== 0.5, :].sec_per_step)
    dacc = sort(unique(result.dacc), rev = true)
    # # println(dacc)

    time = Array{Float64, 1}(undef, length(dacc))
    for ii in 1:length(dacc)
        data = result[result.dacc .== dacc[ii], :].sec_per_step
        time[ii] = minimum(data)
    end
    # println(time)
    # return dacc, time
    return time
end


using PyPlot
include("../util/pyplot.jl")


function main()
    # read options
    argv = parse_cmd()
    output_pdf = argv["pdf"]
    show_ratio = argv["speedup"]

    # find the latest series of simulation results
    # series = argv["target"] == "" ? util_find.latest_series() : argv["target"]
    series = "m31"

    # initialize matplotlib
    util_pyplot.config()
    fig = util_pyplot.set_Panel()

    # p100_dacc, p100_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc60")
    # # v100v_dacc, v100v_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc70")
    # v100p_dacc, v100p_time = get_measured_results(series, "../../../a100/201225gothic/adep/past/", "cc70_60")
    # a100p_dacc, a100p_time = get_measured_results(series, "../../../a100/201225gothic/adep/1410/log/", "cc80")
    # a100s_dacc, a100s_time = get_measured_results(series, "log/", "cc80")
    # # println(a100_dacc)
    # # println(a100_time)

    num = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216]
    a100s = get_measured_results(string(series, "_N001k"), "log/", "cc80")
    a100s = vcat(a100s, get_measured_results(string(series, "_N002k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N004k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N008k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N016k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N032k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N064k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N128k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N256k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N512k"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N001M"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N002M"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N004M"), "log/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series), "lnk/", "cc80"))
    a100s = vcat(a100s, get_measured_results(string(series, "_N016M"), "log/", "cc80"))

    a100p = get_measured_results(string(series, "_001k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80")
    a100p = vcat(a100p, get_measured_results(string(series, "_002k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_004k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_008k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_016k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_032k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_064k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_128k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_256k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_512k"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_001M"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_002M"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_004M"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_008M"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))
    a100p = vcat(a100p, get_measured_results(string(series, "_016M"), "../../../a100/201225gothic/ndep/1410/log/", "cc80"))

    v100p_Ntot, v100p_time = read_csv("v100_sxm2_ref.csv")


    # fig.ax[begin].plot(p100_dacc, p100_time, util_pyplot.call(fig.point, id = 2), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 3), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 0), label = "P100 (SXM)")
    # # fig.ax[begin].plot(v100v_dacc, v100v_time, util_pyplot.call(fig.point, id = 1), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 1), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 2), label = L"V100 (SXM2, compute\_70)", markerfacecolor = "none")
    fig.ax[begin].plot(v100p_Ntot, v100p_time, util_pyplot.call(fig.point, id = 1), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 2), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 0), label = "V100 (SXM2)")
    fig.ax[begin].plot(num, a100p, util_pyplot.call(fig.point, id = 0), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 1), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 2), label = "A100 (PCIe)", markerfacecolor = "none")
    fig.ax[begin].plot(num, a100s, util_pyplot.call(fig.point, id = 0), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 0), linewidth = fig.lw, color = util_pyplot.call(fig.color, id = 1), label = "A100 (SXM4)")

    fig.ax[begin].set_xlabel(L"$N_\mathrm{tot}$", fontsize = fig.fs)
    fig.ax[begin].set_ylabel(L"$t_\mathrm{step}~\bqty{\si{s}}$", fontsize = fig.fs)
    fig.ax[begin].grid()
    fig.ax[begin].loglog()

    # add legend
    handles, labels = fig.ax[begin].get_legend_handles_labels()
	fig.ax[begin].legend(reverse(handles), reverse(labels), numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig.fs)

    # save figures
    fig.fig.savefig(string("fig/", series, "_ndep", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, "_ndep", ".pdf"), format = "pdf", bbox_inches = "tight")
    end

    if show_ratio
        ratio = util_pyplot.set_Panel()
        form = util_pyplot.set_Panel()

        ratio.ax[begin].axhline(1555.0 / 900.0, linestyle = util_pyplot.call(ratio.line, id = 1), linewidth = ratio.lw, color = util_pyplot.call(ratio.color, id = 0), label = "Memory bandwidth")
        ratio.ax[begin].axhline(((64 * 108) * 1.410 * 2) / ((64 * 80) * 1.530 * 2), linestyle = util_pyplot.call(ratio.line, id = 3), linewidth = ratio.lw, color = util_pyplot.call(ratio.color, id = 0), label = "Peak performance")
        ratio.ax[begin].plot(num, v100p_time ./ a100p, util_pyplot.call(ratio.point, id = 1), markersize = ratio.ms, linestyle = util_pyplot.call(ratio.line, id = 2), linewidth = ratio.lw, color = util_pyplot.call(ratio.color, id = 2), label = "A100 (PCIe)")
        ratio.ax[begin].plot(num, v100p_time ./ a100s, util_pyplot.call(ratio.point, id = 0), markersize = ratio.ms, linestyle = util_pyplot.call(ratio.line, id = 0), linewidth = ratio.lw, color = util_pyplot.call(ratio.color, id = 1), label = "A100 (SXM4)")
        form.ax[begin].plot(num, a100p ./ a100s, util_pyplot.call(form.point, id = 0), markersize = form.ms, linestyle = util_pyplot.call(form.line, id = 0), linewidth = form.lw, color = util_pyplot.call(form.color, id = 0))

        ratio.ax[begin].set_xlabel(L"$N_\mathrm{tot}$", fontsize = ratio.fs)
        ratio.ax[begin].set_ylabel("Speed up", fontsize = ratio.fs)
        ratio.ax[begin].grid()
        ratio.ax[begin].semilogx()
        form.ax[begin].set_xlabel(L"$N_\mathrm{tot}$", fontsize = form.fs)
        form.ax[begin].set_ylabel(L"$t_\mathrm{A100 (PCIe)} / t_\mathrm{A100 (SXM4)}$", fontsize = form.fs)
        form.ax[begin].grid()
        form.ax[begin].semilogx()

        # add legend
        handles, labels = ratio.ax[begin].get_legend_handles_labels()
        ratio.ax[begin].legend(reverse(handles), reverse(labels), numpoints = 1, handlelength = 2.0, loc = "best", fontsize = ratio.fs)

        ratio.fig.savefig(string("fig/", series, "_ndep", "_speedup", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
        form.fig.savefig(string("fig/", series, "_ndep", "_form", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
        if output_pdf
            ratio.fig.savefig(string("fig/", series, "_ndep", "_speedup", ".pdf"), format = "pdf", bbox_inches = "tight")
            form.fig.savefig(string("fig/", series, "_ndep", "_form", ".pdf"), format = "pdf", bbox_inches = "tight")
        end

        ratio = nothing
        form = nothing
    end


	fig = nothing
    PyPlot.close("all")
    return nothing

end


main()

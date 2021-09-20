using ArgParse
function parse_cmd()
    cfg = ArgParseSettings()
    @add_arg_table cfg begin
        "--target", "-t"
        help = "filename of the target files"
        arg_type = String
        default = "fiducial/bench/m31_11.time00068797.bare.dat" # 2.748108e-02 is 57th line (out of 64 lines) in fiducial/log/m31_11.acceleration.block.cc80.8388608.time.log
        "--pdf", "-p"
        help = "generate figure in PDF format"
        action = :store_true
    end
    return parse_args(cfg)
end

using DataFrames
using CSV
function read_tsv(file::String)
    # step calcGrav_dev calcMome_dev makeTree setTimeStep_dev sortPHcurve cpBody_dev2hst cpBody_hst2dev prediction_dev correction_dev setLaneTime_dev adjustTime_dev examineNeighbor_dev searchNeighbor_dev
    # tsv = CSV.read(file, DataFrame; delim = "\t", header = ["dacc", "sec", "step", "sec_per_step"], comment = "#")
    tsv = CSV.read(file, DataFrame; delim = "\t", header = ["step", "calcGrav_dev", "calcMome_dev", "makeTree", "setTimeStep_dev", "sortPHcurve", "cpBody_dev2hst", "cpBody_hst2dev", "prediction_dev", "correction_dev", "setLaneTime_dev", "adjustTime_dev", "examineNeighbor_dev", "searchNeighbor_dev"], comment = "#")
    # tsv = CSV.read(file, DataFrame; delim = "\t")
    # println(tsv)
    return tsv
end


using PyPlot
include("../util/pyplot.jl")


using FilePaths
using Glob
function main()
    # read options
    argv = parse_cmd()
    target = argv["target"]
    output_pdf = argv["pdf"]
    series = filename(Path(target))
    # println(series)

    # initialize matplotlib
    util_pyplot.config()
    fig = util_pyplot.set_Panel(xscale = 2.5f0)

    dat = read_tsv(target)
    fig.ax[begin].plot(dat.step[2:end], dat.correction_dev[2:end], "+", markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 0), linewidth = 0.5 * fig.lw, color = util_pyplot.call(fig.color, id = 0), label = "corrector")
    fig.ax[begin].plot(dat.step, dat.prediction_dev, util_pyplot.call(fig.point, id = 4), markersize = fig.ms, color = util_pyplot.call(fig.color, id = 0), label = "predictor")
    fig.ax[begin].plot(dat.step, dat.calcMome_dev, util_pyplot.call(fig.point, id = 3), markersize = fig.ms, color = util_pyplot.call(fig.color, id = 4), label = "calc MAC")
    fig.ax[begin].plot(dat.step, dat.sortPHcurve, util_pyplot.call(fig.point, id = 2), markersize = fig.ms, color = util_pyplot.call(fig.color, id = 2), label = "PH-key")
    fig.ax[begin].plot(dat.step, dat.makeTree, util_pyplot.call(fig.point, id = 1), markersize = fig.ms, color = util_pyplot.call(fig.color, id = 2), label = "make tree")
    fig.ax[begin].plot(dat.step, dat.calcGrav_dev, util_pyplot.call(fig.point, id = 0), markersize = fig.ms, linestyle = util_pyplot.call(fig.line, id = 0), linewidth = 0.5 * fig.lw, color = util_pyplot.call(fig.color, id = 1), label = "walk tree")

    fig.ax[begin].set_xlim([-0.7, 100.7])

    fig.ax[begin].set_xlabel("Step", fontsize = fig.fs)
    fig.ax[begin].set_ylabel(L"Execution time~$\bqty{\si{s}}$", fontsize = fig.fs)
    fig.ax[begin].grid()
    fig.ax[begin].semilogy()

    # add legend
    handles, labels = fig.ax[begin].get_legend_handles_labels()
	# fig.ax[begin].legend(reverse(handles), reverse(labels), numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig.fs)
	fig.ax[begin].legend(reverse(handles), reverse(labels), numpoints = 1, handlelength = 2.0, bbox_to_anchor = (1.0, 1.0), loc = "upper left", fontsize = fig.fs)

    # save figures
    fig.fig.savefig(string("fig/", series, "_breakdown", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, "_breakdown", ".pdf"), format = "pdf", bbox_inches = "tight")
    end


	fig = nothing
    PyPlot.close("all")
    return nothing

end


main()

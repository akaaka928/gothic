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

# using Parameters
# @with_kw mutable struct Conservatives
# 	dt::Array{Real, 1}
# 	err_Etot::Array{Real, 1}
# 	worst_err_Etot::Array{Real, 1}
# 	err_px::Array{Real, 1}
# 	err_py::Array{Real, 1}
# 	err_pz::Array{Real, 1}
# 	worst_err_px::Array{Real, 1}
# 	worst_err_py::Array{Real, 1}
# 	worst_err_pz::Array{Real, 1}
# 	err_Lx::Array{Real, 1}
# 	err_Ly::Array{Real, 1}
# 	err_Lz::Array{Real, 1}
# 	worst_err_Lx::Array{Real, 1}
# 	worst_err_Ly::Array{Real, 1}
# 	worst_err_Lz::Array{Real, 1}
# end


using FilePaths
using Glob
function get_file_list(series::String, lev::Integer)
    list = glob(string("bench/", series, "_lev", lev, "_*.time*.mean.dat"))
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


# include("../util/find.jl")

using PyPlot
include("../util/pyplot.jl")


function main()
    # read options
    argv = parse_cmd()
    output_pdf = argv["pdf"]

    # find the latest series of simulation results
    # series = argv["target"] == "" ? util_find.latest_series() : argv["target"]
    series = "m31"
    l2lev = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    # initialize matplotlib
    # util_pyplot.config(pkg = "\\usepackage{physics,siunitx,amsmath}")
    util_pyplot.config()

    # show dependence on L2lev
    fig = util_pyplot.set_Panel()
    box = util_pyplot.set_Panel()
    # for lev in [0]
    for lev in l2lev
        # generate file list
        list = get_file_list(series, lev)

        # initialize array for results
        elapse = zeros(Float32, length(list))
        # read measured time in for loop
        for ii in 1:length(list)
            elapse[ii] = read_tsv(list[ii])
        end
        # println(elapse)

        # show box-and-whisker plot
        # fig.ax[begin].plot(elapse, elapse, util_pyplot.call(fig.point, id = 0), color = util_pyplot.call(fig.color, id = 0), markersize = fig.ms, label = L"$l = 0$")
        style_cap = Dict("color" => util_pyplot.call(fig.color, id = 0))
        style_box = Dict("color" => util_pyplot.call(fig.color, id = 0), "lw" => 0.75 * fig.lw)
        style_whi = Dict("color" => util_pyplot.call(fig.color, id = 0), "markersize" => fig.ms)
        style_med = Dict("color" => util_pyplot.call(fig.color, id = 1), "lw" => 0.5 * fig.lw)
        fig.ax[begin].boxplot(elapse, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , flierprops = Dict("color" => util_pyplot.call(fig.color, id = 0), "marker" => "x", "markersize" => 0.5 * fig.ms)
        )
        box.ax[begin].boxplot(elapse, positions = [lev], widths = [0.5]
        , capprops = style_cap
        , boxprops = style_box
        , whiskerprops = style_whi
        , medianprops = style_med
        , showfliers = false
        )
        elapse = nothing
    end
    fig.ax[begin].set_xlabel(L"$N_\mathrm{lev}$", fontsize = fig.fs)
    box.ax[begin].set_xlabel(L"$N_\mathrm{lev}$", fontsize = box.fs)
    fig.ax[begin].set_ylabel(L"$\overline{t_\mathrm{walk}}~\bqty{\si{s}}$", fontsize = fig.fs)
    box.ax[begin].set_ylabel(L"$\overline{t_\mathrm{walk}}~\bqty{\si{s}}$", fontsize = box.fs)
    fig.ax[begin].grid()
    box.ax[begin].grid()



    # read TSV file
    # dat = read_csv(string("dat/", series, "_nbody.csv"))
    # read "bench/m31_lev?_?.time*.mean.dat": TSV file, first line is header, 1st value is calcGrav_dev

    # # show error scaling of the energy conservation
    # fig_ene = util_pyplot.set_Panel()
    # fig_ene.ax[begin].plot(dat.dt, abs.(dat.worst_err_Etot), util_pyplot.call(fig_ene.point, id = 1), color = util_pyplot.call(fig_ene.color, id = 0), markersize = fig_ene.ms, label = L"$t = t_\mathrm{worst}$")
    # fig_ene.ax[begin].plot(dat.dt, abs.(dat.      err_Etot), util_pyplot.call(fig_ene.point, id = 0), color = util_pyplot.call(fig_ene.color, id = 1), markersize = fig_ene.ms, label = L"$t = t_\mathrm{final}$")
    # fig_ene.ax[begin].set_xlabel(L"$\varDelta t$", fontsize = fig_ene.fs)
    # fig_ene.ax[begin].set_ylabel(L"$\abs{E(t) / E(t = 0) - 1}$", fontsize = fig_ene.fs)
    # fig_ene.ax[begin].loglog()
    # handles, labels = fig_ene.ax[begin].get_legend_handles_labels()
    # fig_ene.ax[begin].legend(handles[end:-1:begin], labels[end:-1:begin], numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig_ene.fs)

    # # show error scaling of the linear-momentum conservation
    # fig_mom = util_pyplot.set_Panel()
    # fig_mom.ax[begin].plot(dat.dt, dat.      err_px, util_pyplot.call(fig_mom.point, id = 0), color = util_pyplot.call(fig_mom.color, id = 0), markersize = fig_mom.ms, label = L"$p_x (t = t_\mathrm{final})$")
    # fig_mom.ax[begin].plot(dat.dt, dat.      err_py, util_pyplot.call(fig_mom.point, id = 1), color = util_pyplot.call(fig_mom.color, id = 1), markersize = fig_mom.ms, label = L"$p_y (t = t_\mathrm{final})$")
    # fig_mom.ax[begin].plot(dat.dt, dat.      err_pz, util_pyplot.call(fig_mom.point, id = 2), color = util_pyplot.call(fig_mom.color, id = 2), markersize = fig_mom.ms, label = L"$p_z (t = t_\mathrm{final})$")
    # fig_mom.ax[begin].plot(dat.dt, dat.worst_err_px, util_pyplot.call(fig_mom.point, id = 0), color = util_pyplot.call(fig_mom.color, id = 0), markersize = fig_mom.ms, label = L"$p_x (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_mom.ax[begin].plot(dat.dt, dat.worst_err_py, util_pyplot.call(fig_mom.point, id = 1), color = util_pyplot.call(fig_mom.color, id = 1), markersize = fig_mom.ms, label = L"$p_y (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_mom.ax[begin].plot(dat.dt, dat.worst_err_pz, util_pyplot.call(fig_mom.point, id = 2), color = util_pyplot.call(fig_mom.color, id = 2), markersize = fig_mom.ms, label = L"$p_z (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_mom.ax[begin].set_xlabel(L"$\varDelta t$", fontsize = fig_mom.fs)
    # fig_mom.ax[begin].set_ylabel(L"$p(t) - p(t = 0)$", fontsize = fig_mom.fs)
    # fig_mom.ax[begin].yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(util_pyplot.scientific))
    # fig_mom.ax[begin].semilogx()
    # fig_mom.ax[begin].grid()
    # handles, labels = fig_mom.ax[begin].get_legend_handles_labels()
    # fig_mom.ax[begin].legend(handles, labels, numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig_mom.fs)

    # # show error scaling of the angular-momentum conservation
    # fig_spn = util_pyplot.set_Panel()
    # fig_spn.ax[begin].plot(dat.dt, dat.      err_Lx, util_pyplot.call(fig_spn.point, id = 0), color = util_pyplot.call(fig_spn.color, id = 0), markersize = fig_spn.ms, label = L"$L_x (t = t_\mathrm{final})$")
    # fig_spn.ax[begin].plot(dat.dt, dat.      err_Ly, util_pyplot.call(fig_spn.point, id = 1), color = util_pyplot.call(fig_spn.color, id = 1), markersize = fig_spn.ms, label = L"$L_y (t = t_\mathrm{final})$")
    # fig_spn.ax[begin].plot(dat.dt, dat.      err_Lz, util_pyplot.call(fig_spn.point, id = 2), color = util_pyplot.call(fig_spn.color, id = 2), markersize = fig_spn.ms, label = L"$L_z (t = t_\mathrm{final})$")
    # fig_spn.ax[begin].plot(dat.dt, dat.worst_err_Lx, util_pyplot.call(fig_spn.point, id = 0), color = util_pyplot.call(fig_spn.color, id = 0), markersize = fig_spn.ms, label = L"$L_x (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_spn.ax[begin].plot(dat.dt, dat.worst_err_Ly, util_pyplot.call(fig_spn.point, id = 1), color = util_pyplot.call(fig_spn.color, id = 1), markersize = fig_spn.ms, label = L"$L_y (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_spn.ax[begin].plot(dat.dt, dat.worst_err_Lz, util_pyplot.call(fig_spn.point, id = 2), color = util_pyplot.call(fig_spn.color, id = 2), markersize = fig_spn.ms, label = L"$L_z (t = t_\mathrm{worst})$", markerfacecolor = "none")
    # fig_spn.ax[begin].set_xlabel(L"$\varDelta t$", fontsize = fig_spn.fs)
    # fig_spn.ax[begin].set_ylabel(L"$L(t) - L(t = 0)$", fontsize = fig_spn.fs)
    # fig_spn.ax[begin].yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(util_pyplot.scientific))
    # fig_spn.ax[begin].semilogx()
    # fig_spn.ax[begin].grid()
    # handles, labels = fig_spn.ax[begin].get_legend_handles_labels()
    # fig_spn.ax[begin].legend(handles, labels, numpoints = 1, handlelength = 2.0, loc = "best", fontsize = fig_spn.fs)

    # save figures
    fig.fig.savefig(string("fig/", series, "_l2dep", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    box.fig.savefig(string("fig/", series, "_l2box", ".png"), format = "png", dpi = 100, bbox_inches = "tight")
    if output_pdf
        fig.fig.savefig(string("fig/", series, "_l2dep", ".pdf"), format = "pdf", bbox_inches = "tight")
        box.fig.savefig(string("fig/", series, "_l2box", ".pdf"), format = "pdf", bbox_inches = "tight")
    end

	fig = nothing
    PyPlot.close("all")
    return nothing
end


main()




# f0705 = pandas.read_csv('0705/log/' + particle + series + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
# f0705_a00 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 0]]
# f0705_a01 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 1]]
# f0705_a02 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 2]]
# f0705_a03 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 3]]
# f0705_a04 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 4]]
# f0705_a05 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 5]]
# f0705_a06 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 6]]
# f0705_a07 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 7]]
# f0705_a08 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 8]]
# f0705_a09 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 9]]
# f0705_a10 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[10]]
# f0705_a11 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[11]]
# f0705_a12 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[12]]
# f0705_a13 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[13]]
# f0705_a14 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[14]]
# f0705_a15 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[15]]
# f0705_a16 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[16]]
# f0705_a17 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[17]]
# f0705_a18 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[18]]
# f0705_a19 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[19]]

# style_cap = dict(color = col[5])
# style_box = dict(color = col[5], lw = 0.75 * lw_box)
# style_whi = dict(color = col[5], markersize = ms_box)
# style_med = dict(color = col[5], lw = 0.5 * lw_box)
# style_fli = dict(markeredgecolor = col[5], marker = 'x', markersize = 0.5 * ms_box)
# ax_box[0].boxplot([f0705_a00, f0705_a01, f0705_a02, f0705_a03, f0705_a04, f0705_a05, f0705_a06, f0705_a07, f0705_a08, f0705_a09, f0705_a10, f0705_a11, f0705_a12, f0705_a13, f0705_a14, f0705_a15, f0705_a16, f0705_a17, f0705_a18, f0705_a19], positions = acc, widths = numpy.array(acc) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

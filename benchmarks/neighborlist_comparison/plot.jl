# Usage: julia plot.jl  (reads results/*.json, writes docs/src/assets/b_nl_cd.png and b_nl_cv.png)
using Plots, JSON3, LaTeXStrings

const here = @__DIR__
const cfg = JSON3.read(read(joinpath(here, "common_sizes.json"), String))

res = Dict{Int,Dict{String,Any}}()
for t in (1, 8)
    res[t] = Dict{String,Any}()
    for lang in ("julia", "python")
        f = joinpath(here, "results", "$(lang)_$(t)t.json")
        for (name, r) in JSON3.read(read(f, String), Dict{String,Any})
            res[t][name] = (; lang, version=r["version"], cd=r["cd"], cv=r["cv"])
        end
    end
end

# fixed color per package (identity), line style per language
const colors = [
    "CellListMap" => "#2a78d6", "NearestNeighbors" => "#eb6834", "PointNeighbors" => "#1baf7a",
    "scipy" => "#eda100", "vesin" => "#e87ba4", "freud" => "#008300", "matscipy" => "#6250d6",
]

label(name, r) = "$name$(r.lang == "julia" ? ".jl" : "") $(r.version)"
fmt(x) = isinteger(x) ? string(Int(x)) : string(x)

# GR's LaTeX renderer drops plain spaces, so they are escaped
tex(s) = latexstring(replace(s, " " => "\\ "))
const titles = Dict(
    "cd" => tex("\\mathrm{Constant density (}\\rho = $(cfg.density) \\mathrm{atoms/\\AA}^3\\mathrm{, water), cutoff = $(fmt(cfg.cutoff)) \\AA}"),
    "cv" => tex("\\mathrm{Constant volume (}V = $(fmt(cfg.cv_side))^3 \\mathrm{\\AA}^3\\mathrm{), cutoff = $(fmt(cfg.cutoff)) \\AA}"),
)

default(; fontfamily="Computer Modern", grid=true, gridcolor="#e4e4e2", gridalpha=1, gridlinewidth=0.8,
    minorgrid=false, framestyle=:box, foreground_color_axis="#888", foreground_color_border="#888",
    tickfontcolor="#444", guidefontcolor="#222", tickfontsize=8, guidefontsize=10)

for case in ("cd", "cv")
    panels = map((1, 8)) do t
        p = plot(; xscale=:log10, yscale=:log10, legend=false, title=t == 1 ? "Serial" : "8 threads",
            titlefontsize=11, xlabel="Number of particles", ylabel=t == 1 ? "Time to build the neighbor list / s" : "")
        for (name, color) in colors
            # scipy has no parallel pair search (only full lists through query_ball_point)
            (haskey(res[t], name) && !(name == "scipy" && t > 1)) || continue
            r = res[t][name]
            data = getproperty(r, Symbol(case))
            plot!(p, first.(data), getindex.(data, 2); color, linewidth=name == "CellListMap" ? 2.6 : 1.6,
                linestyle=r.lang == "julia" ? :solid : :dash, marker=:circle, markersize=3.5,
                markerstrokewidth=0, markercolor=color)
        end
        t > 1 && annotate!(p, (0.98, 0.04), text("scipy: no parallel pair search", 7, "#666", :right, "Computer Modern"))
        p
    end
    # shared y range, as in sharey
    ylims = extrema(vcat(collect.(Plots.ylims.(panels))...))
    decades = floor(Int, log10(ylims[1])):ceil(Int, log10(ylims[2]))
    foreach(p -> plot!(p; ylims, yticks=10.0 .^ decades), panels)
    # legend-only panel built from the serial results, which have every package
    leg = plot(; framestyle=:none, legend=:top, legend_column=4, foreground_color_legend=nothing,
        legendfontsize=8, grid=false, xlims=(0, 1), ylims=(0, 1))
    for (name, color) in colors
        haskey(res[1], name) || continue
        r = res[1][name]
        plot!(leg, [NaN], [NaN]; color, linewidth=name == "CellListMap" ? 2.6 : 1.6,
            linestyle=r.lang == "julia" ? :solid : :dash, marker=:circle, markersize=3.5,
            markerstrokewidth=0, markercolor=color, label=label(name, r))
    end
    fig = plot(panels..., leg; layout=@layout([a b; c{0.1h}]), size=(1100, 480),
        plot_title=titles[case], plot_titlefontsize=12, plot_titlevspan=0.08, left_margin=5Plots.mm, bottom_margin=3Plots.mm, dpi=130)
    out = joinpath(here, "..", "..", "docs", "src", "assets", "b_nl_$case.png")
    savefig(fig, out)
    println("wrote ", basename(out))
end

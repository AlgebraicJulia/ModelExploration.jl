using DataFrames
using CSV
using Plots; plotlyjs()
using AlgebraicPetri
theme(:ggplot2)

f = "Computer Modern"

sample_df = CSV.read("sample_data.csv", DataFrame)

function make_traj_plot(sample_df, fn, plt_name)
    sample_times = sample_df.times
    S_samples = sample_df.S_samples
    I_samples = sample_df.I_samples
    R_samples = sample_df.R_samples
    D_samples = sample_df.D_samples
    sample_data = hcat(S_samples, I_samples, R_samples, D_samples)

    p = palette(:default)
    colors = [p[1] p[14] p[2] p[4]]

    plt = scatter(sample_times, sample_data, 
        #thickness_scaling=3,
        label="",
        title=plt_name,
        titlefontsize=30,
        xlabel="Time",
        ylabel="People",
        yticks=0:200:600,
        fontfamily=f,
        markercolor=colors
    )

    df = CSV.read(fn, DataFrame)
    times = df.times
    S_vals = df.S_vals
    I_vals = df.I_vals
    R_vals = df.R_vals
    D_vals = df.D_vals
    vals = hcat(S_vals, I_vals, R_vals, D_vals)

    plot!(times, vals, label="", linewidth=2,
        linecolor=colors
    )
    return plt
end


fnames = [
    "SIR_traj.csv",
    "SIR2_traj.csv",
    "SIRS_traj.csv",
    "SIRS2_traj.csv",
    "SIRD_traj.csv",
    "SIRD2_traj.csv",
    "SIRSD_traj.csv",
    "SIRSD2_traj.csv"
]

mnames = [
    "SIR", "SIR 2-City",
    "SIRS", "SIRS 2-City",
    "SIRD", "SIRD 2-City",
    "SIRSD", "SIRSD 2-City"
]

plts = map(x->make_traj_plot(sample_df, x...), zip(fnames, mnames))
scalefontsizes(2.5)
display.(plts)



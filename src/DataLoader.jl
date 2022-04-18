using CSV
using DataFrames
using Plots

fl_df = CSV.read("covid_florida.csv", DataFrame)

omicron_start_index = findfirst(x -> x == "Dec  1 2021", fl_df[:, :Date])
omicron_end_index = findfirst(x -> x == "Mar 20 2022", fl_df[:, :Date])

plt1 = scatter(fl_df[omicron_start_index:-1:omicron_end_index, Symbol("7-Day Moving Avg")])

fl_inf_data = fl_df[omicron_start_index:-1:omicron_end_index, Symbol("7-Day Moving Avg")]

const FL_POP = 20e6

fl_susc_data = Int64[]

prev_val = FL_POP
for inf_val in fl_inf_data
    push!(fl_susc_data, prev_val - inf_val)
    global prev_val -= inf_val
end

plt2 = scatter(fl_susc_data)
l = @layout [a; b]
plot(plt1, plt2, layout=l)

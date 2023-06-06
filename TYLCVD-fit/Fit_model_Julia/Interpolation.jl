using CSV, DataFrames
using DataInterpolations, Plots
N_p = 1000
df_pscl_4 = CSV.read("data_PSCL-4.csv", DataFrame)
#interpolation:
t_end = df_pscl_4.t[end]
timeLine = 0:t_end
scaled_i_p = N_p * df_pscl_4.y
obs_t = df_pscl_4.t
itp = CubicSpline(scaled_i_p, obs_t)
scatter(obs_t, scaled_i_p)
plot!(itp)

interpolated_infected = [itp(t) for t in timeLine]
interpolated_df = DataFrame(
    time = timeLine,
    interpolated_ip = interpolated_infected
)

CSV.write("interpolated_data_PSCL-4.csv", interpolated_df)
plot(interpolated_df.time, interpolated_df.interpolated_ip)
scatter!(interpolated_df.time, interpolated_df.interpolated_ip)

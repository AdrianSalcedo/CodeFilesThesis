using Plots; plotly()
using CSV
using IterableTables, DataFrames, DataTables
using Statistics

path1 = "/home/gabrielsalcedo/Documentos/Julia_code_for_tomate_SDE_paper-main-master/"
path2 = "Persistence_Rs0_noise_additing/"
path0 = path1 * path2
Data = CSV.read(path0 * "Data_noise_incidence_Tyking_additing.csv", DataFrame)

Data_Mean = DataFrame()

for i in 27:70
    Sub = Data[Data.t .== i, :]
    Mean_t = mean(Sub.t)
    # Mean_S_p = mean(Sub.Sp)
    # Mean_I_p = mean(Sub.Ip)
    # Mean_I_v = mean(Sub.Iv)
     Mean_Y_p = mean(Sub.Yp)
    # Quartile_I_p = quantile!(Sub.Ip, [0.05, 0.25, 0.50, 0.75, 0.95])
    # Quartile_S_p = quantile!(Sub.Sp, [0.05, 0.25, 0.50, 0.75, 0.95])
    # Quartile_I_v = quantile!(Sub.Iv, [0.05, 0.25, 0.50, 0.75, 0.95])
    Quartile_Y_p = quantile!(Sub.Yp, [0.025, 0.25, 0.50, 0.75, 0.975])
    Data_aux = DataFrame(
    t = Mean_t,
    # Q05_S_p = Quartile_S_p[1],
    # Q25_S_p = Quartile_S_p[2],
    # Q50_S_p = Quartile_S_p[3],
    # Mean_S_p = Mean_S_p,
    # Q75_S_p = Quartile_S_p[4],
    # Q95_S_p = Quartile_S_p[5],
    # Q05_I_p = Quartile_I_p[1],
    # Q25_I_p = Quartile_I_p[2],
    # Q50_I_p = Quartile_I_p[3],
    # Mean_I_p = Mean_I_p,
    # Q75_I_p = Quartile_I_p[4],
    # Q95_I_p = Quartile_I_p[5],
    # Q05_I_v = Quartile_I_v[1],
    # Q25_I_v = Quartile_I_v[2],
    # Q50_I_v = Quartile_I_v[3],
    # Mean_I_v = Mean_I_v,
    # Q75_I_v = Quartile_I_v[4],
    # Q95_I_v = Quartile_I_v[5],
    Q05_Y_p = Quartile_Y_p[1],
    Q25_Y_p = Quartile_Y_p[2],
    Q50_Y_p = Quartile_Y_p[3],
    Mean_Y_p = Mean_Y_p,
    Q75_Y_p = Quartile_Y_p[4],
    Q95_Y_p = Quartile_Y_p[5],
      )
    Data_Mean = append!(Data_Mean,Data_aux)
end
CSV.write(path0 * "Data_mean_persistence_incidence_resistance_Tyking_add.csv",Data_Mean)

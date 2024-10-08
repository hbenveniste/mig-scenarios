using CSV, DataFrames, DelimitedFiles
using Plots, VegaLite, FileIO, VegaDatasets, FilePaths
using Statistics, Query, Distributions, StatsPlots


# Read data on calculated SSP population and GDP with and without migration
sspall = CSV.File(joinpath(@__DIR__, "../results/sspall_6_update.csv")) |> DataFrame


################################################# Compute Gini coefficients between countries with and without migration ############################################
# Treat separately version with and without migration
ginibtw_mig = sspall[:,[:scen, :country, :period, :pop_mig, :gdp_mig, :ypc_mig]]
ginibtw_nomig = sspall[:,[:scen, :country, :period, :pop_nomig, :gdp_nomig, :ypc_nomig]]

periods = unique(ginibtw_mig[!,:period]) ; ssps = unique(ginibtw_mig[!,:scen]) ; countries = unique(ginibtw_mig[!,:country])

# Order countries per growing ypc 
sort!(ginibtw_mig, [:scen, :period, :ypc_mig])
sort!(ginibtw_nomig, [:scen, :period, :ypc_nomig])

# Compute cumulative share of people
pop_world_mig = rename(combine(groupby(ginibtw_mig, [:scen, :period]), :pop_mig => sum),:pop_mig_sum => :pop_world_mig)
cum_pop=[]
for i in 0:length(periods)*length(ssps)-1
    cp = []
    for c in 1:length(countries)
        s = sum(ginibtw_mig[!,:pop_mig][i*length(countries)+1:i*length(countries)+c])/pop_world_mig[!,:pop_world_mig][i+1]
        append!(cp, s)
    end
    append!(cum_pop, cp)
end
ginibtw_mig[!,:cum_pop] = cum_pop

pop_world_nomig = rename(combine(groupby(ginibtw_nomig, [:scen, :period]), :pop_nomig => sum),:pop_nomig_sum => :pop_world_nomig)
cum_pop=[]
for i in 0:length(periods)*length(ssps)-1
    cp = []
    for c in 1:length(countries)
        s = sum(ginibtw_nomig[!,:pop_nomig][i*length(countries)+1:i*length(countries)+c])/pop_world_nomig[!,:pop_world_nomig][i+1]
        append!(cp, s)
    end
    append!(cum_pop, cp)
end
ginibtw_nomig[!,:cum_pop] = cum_pop

# Compute cumulative share of income
gdp_world_mig = rename(combine(groupby(ginibtw_mig, [:scen, :period]), :gdp_mig => sum), :gdp_mig_sum => :gdp_world_mig)
cum_gdp=[]
for i in 0:length(periods)*length(ssps)-1
    cp = []
    for c in 1:length(countries)
        s = sum(ginibtw_mig[!,:gdp_mig][i*length(countries)+1:i*length(countries)+c])/gdp_world_mig[!,:gdp_world_mig][i+1]
        append!(cp, s)
    end
    append!(cum_gdp, cp)
end
ginibtw_mig[!,:cum_gdp] = cum_gdp

gdp_world_nomig = rename(combine(groupby(ginibtw_nomig, [:scen, :period]), :gdp_nomig => sum),:gdp_nomig_sum => :gdp_world_nomig)
cum_gdp=[]
for i in 0:length(periods)*length(ssps)-1
    cp = []
    for c in 1:length(countries)
        s = sum(ginibtw_nomig[!,:gdp_nomig][i*length(countries)+1:i*length(countries)+c])/gdp_world_nomig[!,:gdp_world_nomig][i+1]
        append!(cp, s)
    end
    append!(cum_gdp, cp)
end
ginibtw_nomig[!,:cum_gdp] = cum_gdp

# Compute Gini coefficients based on area under Lorenz curve
gini_world_mig = []
for i in 0:length(periods)*length(ssps)-1
    g_prel = 0.0
    for c in 1:length(countries)-1
        prel = 0.5 * (ginibtw_mig[!,:cum_pop][i*length(countries)+c+1] - ginibtw_mig[!,:cum_pop][i*length(countries)+c]) * (ginibtw_mig[!,:cum_gdp][i*length(countries)+c+1] + ginibtw_mig[!,:cum_gdp][i*length(countries)+c])
        g_prel += prel
    end
    append!(gini_world_mig, 1 - 2 * g_prel)
end

gini_world_nomig = []
for i in 0:length(periods)*length(ssps)-1
    g_prel = 0.0
    for c in 1:length(countries)-1
        prel = 0.5 * (ginibtw_nomig[!,:cum_pop][i*length(countries)+c+1] - ginibtw_nomig[!,:cum_pop][i*length(countries)+c]) * (ginibtw_nomig[!,:cum_gdp][i*length(countries)+c+1] + ginibtw_nomig[!,:cum_gdp][i*length(countries)+c])
        g_prel += prel
    end
    append!(gini_world_nomig, 1 - 2 * g_prel)
end

# Gather calculated between-countries Gini with and without migration
gini_world = DataFrame(scen = gdp_world_mig[!,:scen], period = gdp_world_mig[!,:period])
gini_world[!,:gini_world_mig] = gini_world_mig
gini_world[!,:gini_world_nomig] = gini_world_nomig


################################################# Plot results for world #####################################
gini_w = stack(gini_world, [:gini_world_mig, :gini_world_nomig], [:scen, :period])
rename!(gini_w, :variable => :gini_type, :value => :gini)
gini_w |> @vlplot(
    width=300, height=250,
    mark={:point, size=60}, x = {"period:o", axis={labelFontSize=16}, title=nothing}, y = {"gini:q", title="Between countries Gini", axis={labelFontSize=16,titleFontSize=16}}, 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    shape = {"gini_type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/gini_btw/", "gini_world_6_update.png"))

# Plot relative differences in Gini levels without migration compared to with
gini_world[!,:reldif] = gini_world[!,:gini_world_mig] ./ gini_world[!,:gini_world_nomig] .- 1

gini_world |> @vlplot(
    width=300, height=250,
    mark = {:point, filled=true, size=80}, x = {"period:o", axis={labelFontSize=16}, title=nothing}, y = {"reldif:q", title="Relative change with migration", axis={labelFontSize=16, titleFontSize=16}}, 
    color = {"scen:n", scale={scheme=:category10}, legend=nothing}
) |> save(joinpath(@__DIR__, "../results/gini_btw/", "ginirel_world_6_update.png"))


################################################# Plot Lorenz curves ########################################
ginibtw_mig[!,:type] = repeat(["mig"], size(ginibtw_mig,1))
ginibtw_nomig[!,:type] = repeat(["nomig"], size(ginibtw_nomig,1))
lorenz = vcat(ginibtw_mig[:,[:scen, :country, :period, :cum_gdp, :cum_pop, :type]], ginibtw_nomig[:,[:scen, :country, :period, :cum_gdp, :cum_pop, :type]])
lorenz[(lorenz[!,:period].==2065),:] |> @vlplot(
    :line, x = {"cum_pop:q", axis={labelFontSize=16}}, y = {"cum_gdp:q", axis={labelFontSize=16}}, 
    title = "World Lorenz curve for 2065", 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    shape = {"type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/gini_btw/", "lorenz_2065_world_6_update.png"))


####################################### Write output files with results ############################################################################
CSV.write(joinpath(@__DIR__, "../results/gini_btw_6_update.csv"), gini_world)
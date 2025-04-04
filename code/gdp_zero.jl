using CSV, DataFrames, DelimitedFiles, XLSX
using Plots, VegaLite, FileIO, VegaDatasets, FilePaths
using Statistics, Query


################## Prepare population data: original SSP and no-migration version ####################
# Original version:
# Source:  Wittgenstein Center (WIC) Population and Human Capital Projections, version v.14 (February 2025). 
# https://zenodo.org/records/14718294
ssp1 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/PROJresult_AGE_SSP1_V14.csv", DataFrame)
ssp2 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/PROJresult_AGE_SSP2_V14.csv", DataFrame)
ssp3 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/PROJresult_AGE_SSP3_V14.csv", DataFrame)
ssp4 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/PROJresult_AGE_SSP4_V14.csv", DataFrame)
ssp5 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/PROJresult_AGE_SSP5_V14.csv", DataFrame)

ssp1.scen = repeat(["SSP1"], size(ssp1,1))
ssp2.scen = repeat(["SSP2"], size(ssp2,1))
ssp3.scen = repeat(["SSP3"], size(ssp3,1))
ssp4.scen = repeat(["SSP4"], size(ssp4,1))
ssp5.scen = repeat(["SSP5"], size(ssp5,1))

ssp = vcat(ssp1, ssp2, ssp3, ssp4, ssp5)

# Create net migration variable
ssp[!,:mig] = ssp[!,:emi] .- ssp[!,:imm]

# Version for zero migration
# Provided by Samir KC in March 2025
mig0_ssp1 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/res_x05ssp1_emort_efert_0mig.csv", DataFrame)
mig0_ssp2 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/res_x01ssp2_emort_efert_0mig.csv", DataFrame)
mig0_ssp3 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/res_x07ssp3_emort_efert_0mig.csv", DataFrame)
mig0_ssp4 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/res_x08ssp4_emort_efert_0mig.csv", DataFrame)
mig0_ssp5 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/res_x06ssp5_emort_efert_0mig.csv", DataFrame)

mig0_ssp1.scen = repeat(["SSP1"], size(mig0_ssp1,1))
mig0_ssp2.scen = repeat(["SSP2"], size(mig0_ssp2,1))
mig0_ssp3.scen = repeat(["SSP3"], size(mig0_ssp3,1))
mig0_ssp4.scen = repeat(["SSP4"], size(mig0_ssp4,1))
mig0_ssp5.scen = repeat(["SSP5"], size(mig0_ssp5,1))

mig0_ssp = vcat(mig0_ssp1, mig0_ssp2, mig0_ssp3, mig0_ssp4, mig0_ssp5)
select!(mig0_ssp, [:region,:Time,:sex,:edu,:agest,:pop,:births,:emi,:imm,:deaths,:scen])

# Sum projections for all sexes, ages, and education levels: population + net migration per country and time period
ssp_cy = combine(
    groupby(
        ssp, 
        [:region, :Time, :scen]
    ), 
    d -> (
        pop_mig = sum(d.pop), 
        inmigsum = sum(d.emi), 
        outmigsum = sum(d.imm)
    )
)

mig0_ssp_cy = combine(
    groupby(
        mig0_ssp, 
        [:region, :Time, :scen]
    ), 
    d -> sum(d.pop)
)
rename!(mig0_ssp_cy, :x1 => :pop_nomig)

# Join population datasets
sspall = innerjoin(ssp_cy, mig0_ssp_cy, on = [:region, :Time, :scen])
sspall[!,:region] = map(x -> parse(Int, SubString(x, 4)), sspall[!,:region])

# Create variable for differences in population projections between mig and nomig
sspall[!,:diffmig] = sspall[!,:pop_mig] .- sspall[!,:pop_nomig]   


#################### Prepare GDP data ################################
# Using the GDP versions of Dellink et al. 2017
# Unit: billion US$ 2005 / year PPP
gdp_oecd = XLSX.readdata(joinpath(@__DIR__, "../data/gdp_sspdb/gdp_oecd.xlsx"), "data!A1:X921") 
gdp_oecd = DataFrame(gdp_oecd, :auto)
rename!(gdp_oecd, Symbol.(Vector(gdp_oecd[1,:])))
deleteat!(gdp_oecd,1)
select!(gdp_oecd, Not([:Model, :Variable, :Unit]))
gdps_oecd = stack(gdp_oecd, 3:length(gdp_oecd[1,:]))
rename!(gdps_oecd, :variable => :year, :value => :gdp)
gdps_oecd[!,:year] = map( x -> parse(Int, String(x)), gdps_oecd[!,:year])

# Convert into iso3c country codes
iso3c_isonum = CSV.read(joinpath(@__DIR__,"../data/iso3c_isonum.csv"), DataFrame)
sspall = leftjoin(sspall, rename(iso3c_isonum, :iso3c => :country, :isonum => :region), on = :region)

# Micronesia (Federated states of) ISO numerical code is not 954, but 583. Its ISO 3 code is FSM
[if sspall[i,:region] == 583 ; sspall[i,:country] = "FSM" end for i in 1:length(sspall[!,1])]

sspall = innerjoin(sspall, rename(gdps_oecd, :year => :Time, :Region => :country, :Scenario => :scen), on = [:Time, :country, :scen])         
rename!(sspall, :gdp => :gdp_mig)


###################### Project migrant flows ###################################################
years = unique(sspall[!,:Time]) ; nbyears = length(years)
countries = unique(sspall[!,:country]) ; nbcountries = length(countries)
ssps = unique(sspall[!,:scen]) ; nbssps = length(ssps)
migflow = DataFrame(
    scen = repeat(ssps, inner = nbcountries * nbcountries * nbyears), 
    origin = repeat(countries, inner = nbyears * nbcountries, outer = nbssps), 
    destination = repeat(countries, inner = nbyears , outer = nbcountries * nbssps), 
    Time = repeat(years, outer = nbcountries * nbcountries * nbssps)
)
rename!(migflow, :origin => :country)
migflow = innerjoin(migflow, sspall, on = [:country, :Time, :scen])
select!(migflow, Not([:region, :pop_nomig, :inmigsum, :diffmig]))
rename!(migflow, :country => :origin, :pop_mig => :pop_orig, :outmigsum => :outmig_orig, :gdp_mig => :gdp_orig, :destination => :country)
migflow = innerjoin(migflow, sspall, on = [:country, :Time, :scen])
select!(migflow, Not([:region, :pop_nomig, :outmigsum, :diffmig]))
rename!(migflow, :country => :destination, :pop_mig => :pop_dest, :inmigsum => :inmig_dest, :gdp_mig => :gdp_dest)

distance = CSV.File(joinpath(@__DIR__, "../data/distance_unpop/distances.csv")) |> DataFrame
rename!(distance, :orig => :origin, :dest => :destination)

countries_migflow = unique(migflow[!,:origin])
countries_distance = unique(distance[!,:origin])
setdiff(countries_migflow, countries_distance)

migflow = leftjoin(migflow, distance, on = [:origin, :destination])

# Making units consistent 
migflow[!,:pop_orig] .*= 1000        # pop is in thousands
migflow[!,:pop_dest] .*= 1000
migflow[!,:gdp_orig] .*= 10^9        # gdp is in billion $
migflow[!,:gdp_dest] .*= 10^9

# Creating gdp per capita variables
migflow[!,:ypc_orig] = migflow[!,:gdp_orig] ./ migflow[!,:pop_orig]
migflow[!,:ypc_dest] = migflow[!,:gdp_dest] ./ migflow[!,:pop_dest]
migflow[!,:ypc_ratio] = migflow[!,:ypc_dest] ./ migflow[!,:ypc_orig]

# Joining data on remittances
rho = CSV.File(joinpath(@__DIR__, "../data/rem_wb/rho.csv")) |> DataFrame
phi = CSV.File(joinpath(@__DIR__,"../data/rem_wb/phi.csv")) |> DataFrame
rename!(rho, :rho => :remshare)
rename!(phi, :phi => :remcost)
migflow = leftjoin(migflow, rho, on = [:origin, :destination])
migflow = leftjoin(migflow, phi, on = [:origin, :destination])
# For corridors with no remittance data, attribute 0
migflow[!,:remshare] = Missings.coalesce.(migflow[!,:remshare], 0)
# For corridors with no cost data, assume that the cost of sending migflow is the mean of all available corridors
for i in eachindex(migflow[:,1])
    if ismissing(migflow[i,:remcost])
        migflow[i,:remcost] = (migflow[i,:origin] == migflow[i,:destination] ? 0.0 : mean(phi[!,:remcost]))
    end
end

# Joining data on common official languages
comofflang = CSV.File(joinpath(@__DIR__, "../data/common_official_language/comofflang.csv")) |> DataFrame
rename!(comofflang, :orig => :origin, :dest => :destination)
migflow = leftjoin(migflow, comofflang, on = [:origin,:destination])
# Country not covered in comofflang: TWN, common language with CHN, HKG, SGP, MAC 
for i in eachindex(migflow[:,1])
    if ismissing(migflow[i,:comofflang])
        cmis = ["CHN", "TWN","HKG","SGP","MAC"]
        migflow[i,:comofflang] = (in(cmis, migflow[i,:origin]) && in(cmis, migflow[i,:destination])) ? 1 : 0
    end
end

# Calibration results from gravcalib.jl
# This is calibration based on Azose and Raftery's data for 1990-2015, compiled in Abel and Cohen (2019) 
# Main specification: year fixed effects (reg_ar_yfe). All resulting files and graphs indexed _6.
# Robustness runs: origin/destination/year fixed effects (reg_ar_odyfe). All resulting files and graphs indexed _7.

beta = CSV.File(joinpath(@__DIR__,"../data/gravity_calib/beta.csv")) |> DataFrame
ireg = findfirst(beta[!,:regtype].=="reg_ar_yfe")

migflow[!,:move_prelim] = 
    exp(beta[ireg,:b0]) .* 
    migflow[:,:pop_orig].^beta[ireg,:b1] .* 
    migflow[:,:pop_dest].^beta[ireg,:b2] .* 
    migflow[:,:ypc_orig].^beta[ireg,:b4] .* 
    migflow[:,:ypc_dest].^beta[ireg,:b5] .* 
    migflow[:,:distance].^beta[ireg,:b7] .* 
    map(x -> exp(beta[ireg,:b8] * x), migflow[:,:remshare]) .* 
    map(x -> exp(beta[ireg,:b9] * x), migflow[:,:remcost]) .* 
    map(x -> exp(beta[ireg,:b10] * x), migflow[:,:comofflang])

for i in eachindex(migflow[:,1])
    if migflow[i,:origin] == migflow[i,:destination]
        migflow[i,:move_prelim] = 0
    end
end

# Rescale so that sum(move_prelim) = inmig_dest/5, with move_prelim from the gravity equation and inmig_dest from Samir KC's data
movein = combine(
    groupby(
        migflow, 
        [:destination, :Time, :scen]
    ), 
    d -> (move_in_prelim = sum(d.move_prelim), inmig_dest_ = mean(d.inmig_dest)), 
)
select!(movein, Not(:inmig_dest_))
migflow = innerjoin(migflow, movein, on = [:destination, :Time, :scen])
migflow[!,:move_in_frac] = migflow[!,:move_prelim] ./ migflow[!,:move_in_prelim]
migflow[!,:move_inmigbase] = migflow[!,:move_in_frac] .* migflow[!,:inmig_dest] .* 1000 ./ 5         # data on inmig and outmig flows are for a 5-year period


###################### Prepare data for remittances ###################################################
# Computing age of migrants at time of migration, using Samir KC's projection data
# No migrants numbers for the updated SSPs; we assume that migrants age distributions remain the same for the update
agegroup = combine(groupby(ssp, [:agest, :region, :Time, :scen]), d -> (inmig = sum(d.emi), pop = sum(d.pop)))

iso3c_isonum = CSV.File(joinpath(@__DIR__,"../data/iso3c_isonum.csv")) |> DataFrame
rename!(iso3c_isonum, :iso3c => :country, :isonum => :region)
agegroup[!,:region] = map(x -> parse(Int, SubString(x, 4)), agegroup[!,:region])
agegroup = leftjoin(agegroup, iso3c_isonum, on = :region)
# Micronesia (Federated states of) ISO numerical code is not 954, but 583. Its ISO 3 code is FSM
[if agegroup[i,:region] == 583 ; agegroup[i,:country] = "FSM" end for i in 1:length(agegroup[!,1])]
filter!(
    row -> (row.region < 900),
    agegroup
) 
select!(agegroup, Not(:region))

# Computing shares of migrants by age
ageshare = combine(groupby(agegroup, [:scen, :country, :Time]), d -> sum(d.inmig))
rename!(ageshare, :x1 => :inmig_sum)
agegroup = innerjoin(agegroup, ageshare, on = [:scen, :country, :Time])
agegroup[!,:share] = agegroup[!,:inmig] ./ agegroup[!,:inmig_sum]
for i in eachindex(agegroup[:,1]) ; if agegroup[i,:inmig] == 0 && agegroup[i,:inmig_sum] == 0 ; agegroup[i,:share] = 0 end end

# Introducing life expectancy as projected in the SSP. Data from the Wittgenstein Centre
lifeexp = CSV.File(joinpath(@__DIR__, "../data/lifeexp_sspwc/lifeexp.csv");header=9) |> DataFrame
select!(lifeexp, Not(:Area))
lifeexp[!,:Period] = map( x -> parse(Int, SubString(x, 1:4)), lifeexp[!,:Period])
# Compute life expectancy for overall population as average between male and female
lifeexp = combine(groupby(lifeexp, [:Scenario, :Period, :ISOCode]), d -> mean(d.Years))   
rename!(lifeexp, :x1 => :lifeexp)

rename!(iso3c_isonum, :region => :ISOCode)
lifeexp = leftjoin(lifeexp, iso3c_isonum, on = :ISOCode)
# Micronesia (Federated states of) ISO numerical code is not 954, but 583. Its ISO 3 code is FSM
[if lifeexp[i,:ISOCode] == 583 ; lifeexp[i,:country] = "FSM" end for i in 1:length(lifeexp[!,1])]
# The Channels Islands do not have a proper ISO code, instead Jersey (832, JEY) and Guernsey (831, GGY) do. 
channelsind = findall(lifeexp[!,:ISOCode] .== 830)
for i in channelsind
    push!(lifeexp, [lifeexp[i,:Scenario], lifeexp[i,:Period], 831, lifeexp[i,:lifeexp], "GGY"])
    push!(lifeexp, [lifeexp[i,:Scenario], lifeexp[i,:Period], 832, lifeexp[i,:lifeexp], "JEY"])
end
delete!(lifeexp, channelsind)
delete!(lifeexp, findall(lifeexp[!,:ISOCode] .== 900))        # Delete data for world
select!(lifeexp, Not(:ISOCode))
sort!(lifeexp, [:Scenario, :country])

# We only have data for SSP1,2,3, but based on mortality assumptions mentioned in KC and Lutz (2017), we assume that SSP4 ~ SSP2, and SSP5 ~ SSP1
lifeexp_ssp4 = lifeexp[(lifeexp[!,:Scenario] .== "SSP2"),:]
lifeexp_ssp4[!,:Scenario] = repeat(["SSP4"], size(lifeexp_ssp4,1))
lifeexp_ssp5 = lifeexp[(lifeexp[!,:Scenario] .== "SSP1"),:]
lifeexp_ssp5[!,:Scenario] = repeat(["SSP5"], size(lifeexp_ssp5,1))
lifeexp = vcat(lifeexp, lifeexp_ssp4)
lifeexp = vcat(lifeexp, lifeexp_ssp5)

# Join data
rename!(lifeexp, :Scenario => :scen, :country => :destination, :Period => :Time)
migflow = innerjoin(migflow, lifeexp, on = [:scen, :destination, :Time])


######################## Adding a stock variable indicating how many immigrants from a region are in another region ##############################
# Computing the part of the migrants stock that comes from people who migrate during the projections, i.e. starting 2020
rename!(agegroup, :country => :destination)
migflow_all = innerjoin(
    migflow[:,[:scen, :origin, :destination, :Time, :inmig_dest, :move_in_frac, :move_inmigbase, :lifeexp]], 
    agegroup[:,[:scen, :destination, :Time, :agest, :pop, :share]], 
    on = [:scen, :destination, :Time]
)
l_or = length(unique(migflow_all[!,:origin]))
l_dest = length(unique(migflow_all[!,:destination]))
periods = unique(migflow_all[!,:Time])
l_period = length(periods)
l_age = length(unique(migflow_all[!,:agest]))
l_ssp = length(unique(migflow_all[!,:scen]))
migflow_all[!,:duration_prep] = (migflow_all[!,:Time]) .- (migflow_all[!,:lifeexp] .- migflow_all[!,:agest])

# First, people who migrate year after year
# Compute stocks of migrants for 5-year period, per age group
migflow_all[!,:stock0] = [migflow_all[i, :move_inmigbase] * migflow_all[i,:share] * 5 for i in eachindex(migflow_all[:,1])] 
# For simplification's sake, we do not take into account the migrants of age above life expectancy: they only represent 0.06% of migrants in SSP2
for i in eachindex(migflow_all[:,1])
    if migflow_all[i,:lifeexp] < migflow_all[i,:agest]
        migflow_all[i,:stock0] = 0.0
    end
end

# Second, how each period's migrants stock declines over time
for p in eachindex(periods)
    migflow_all[!,Symbol(string("stock",periods[p]))] = zeros(size(migflow_all,1))
    for i in 0:(l_or*l_dest*l_ssp-1)
        for j in (p-1)*l_age+1:l_age*l_period
            ind = i * l_period * l_age + j
            value = (migflow_all[ind,:duration_prep] - periods[p] < 0) ? migflow_all[i * l_period * l_age + (p-1) * l_age + (rem(j, l_age) == 0 ? l_age : rem(j, l_age)),:stock0] : 0.0
            migflow_all[ind,Symbol(string("stock",periods[p]))] = value
        end
    end
end

# Third, add all periods' migrants stocks
migflow_all[!,:stock_dyn] = migflow_all[!,:stock2020] .+ migflow_all[!,:stock2025] .+ migflow_all[!,:stock2030] .+ 
    migflow_all[!,:stock2035] .+ migflow_all[!,:stock2040] .+ migflow_all[!,:stock2045] .+ migflow_all[!,:stock2050] .+ 
    migflow_all[!,:stock2055] .+ migflow_all[!,:stock2060] .+ migflow_all[!,:stock2065] .+ migflow_all[!,:stock2070] .+ 
    migflow_all[!,:stock2075] .+ migflow_all[!,:stock2080] .+ migflow_all[!,:stock2085] .+ migflow_all[!,:stock2090] .+ migflow_all[!,:stock2095]
# Making room in migflow_all
select!(migflow_all, Not([:stock2020, :stock2025, :stock2030, :stock2035, :stock2040, :stock2045, :stock2050, :stock2055, :stock2060, :stock2065, :stock2070, :stock2075, :stock2080, :stock2085, :stock2090, :stock2095]))

# Computing the initial stock of migrants, and how it declines over time
# We use bilateral migration stocks from 2017 from the World Bank
# In order to get its age distribution, we assume that it is the average of two age distributions in the destination country: 
# the one of migrants at time of migration in the period 2020-2025 (computed from the SSP as "share")
# and the one of the overall destination population in the period 2020-2025

# Reading bilateral migrant stocks from 2017
migstock_matrix = XLSX.readdata(joinpath(@__DIR__, "../data/rem_wb/WB_Bilateral_Estimates_Migrant_Stocks_2017.xlsx"), "Bilateral_Migration_2017!A2:HJ219")
migstock_matrix = DataFrame(migstock_matrix, :auto)
header = 2
countries = migstock_matrix[(header):(length(migstock_matrix[:,1]) - 3), 1]
migstock = DataFrame(
    origin = repeat(countries, inner = length(countries)), 
    destination = repeat(countries, outer = length(countries))
)
stock = []
for o in (header):(length(countries)+header-1)
    ostock = migstock_matrix[o, 2:(end - 3)]
    append!(stock, ostock)
end
migstock.migrantstock = stock
indmissing = findall([typeof(migstock[i,:migrantstock]) != Int64 for i in eachindex(migstock[:,1])])
for i in indmissing ; migstock[i,:migrantstock] = 0.0 end

# Converting into country codes
ccode = XLSX.readdata(joinpath(@__DIR__,"../data/rem_wb/GDPpercap2017.xlsx"), "Data!A1:E218")
ccode = DataFrame(ccode, :auto)
rename!(ccode, Symbol.(Vector(ccode[1,:])))
deleteat!(ccode,1)
select!(ccode, Not([Symbol("Series Code"), Symbol("Series Name"), Symbol("2017 [YR2017]")]))
rename!(ccode, Symbol("Country Name") => :country, Symbol("Country Code") => :country_code)
rename!(ccode, :country => :destination)
indnkorea = findfirst(x -> x == "Korea, Dem. People’s Rep.", ccode[!,:destination])
ccode[indnkorea,:destination] = "Korea, Dem. Rep."
migstock = innerjoin(migstock, ccode, on = :destination)
rename!(migstock, :country_code => :dest_code)
rename!(ccode, :destination => :origin)
migstock = innerjoin(migstock, ccode, on = :origin)
rename!(migstock, :country_code => :orig_code)
select!(migstock, Not([:origin, :destination]))

# Getting age distributions
rename!(agegroup, :destination => :country)
agedist = @from i in agegroup begin
    @where i.Time == 2020
    @select {i.scen, i.country, i.agest, i.pop, migshare = i.share}
    @collect DataFrame
end
popshare = combine(groupby(agedist, [:scen, :country]), d -> sum(d.:pop))
rename!(popshare, :x1 => :pop_sum)
agedist = innerjoin(agedist, popshare, on = [:scen, :country])
agedist[!,:popshare] = agedist[!,:pop] ./ agedist[!,:pop_sum]
select!(agedist, Not([:pop, :pop_sum]))
agedist[!,:meanshare] = (agedist[!,:popshare] .+ agedist[!,:migshare]) ./ 2

# Join data
rename!(migstock, :dest_code => :destination, :orig_code => :origin)
rename!(agedist, :country => :destination)
migstock = innerjoin(migstock, agedist, on = :destination)
migstock[!,:stock_by_age] = migstock[!,:migrantstock] .* migstock[!,:meanshare]
migstock[!,:Time] = repeat([2020], size(migstock,1))

# Making the initial migrant stock decrease over time
stock_by_age = migstock[:,[:scen, :origin, :destination, :agest, :stock_by_age, :Time]]

countries_ms_o = setdiff(unique(migstock[!,:origin]),unique(migflow_all[!,:origin]))
countries_ms_d = setdiff(unique(migstock[!,:destination]),unique(migflow_all[!,:destination]))
indrem = []
for i in eachindex(stock_by_age[:,1])
    if in(stock_by_age[i,:origin], countries_ms_o) == true || in(stock_by_age[i,:destination], countries_ms_d) == true
        append!(indrem, i)
    end
end
delete!(stock_by_age, indrem)

countries_missing = setdiff(unique(migflow_all[!,:origin]),unique(migstock[!,:origin]))     # SWZ, TWN
co_un = vcat(unique(stock_by_age[!,:origin]), ["SWZ", "TWN"])
sort!(co_un)
sba_swz_o = DataFrame(scen = repeat(ssps, inner = l_age*l_or), origin = repeat(["SWZ"], l_ssp*l_age*l_or), destination = repeat(co_un, outer = l_ssp, inner = l_age), agest = repeat(unique(migflow_all[!,:agest]), outer = l_ssp*l_or), stock_by_age = repeat([0], l_ssp*l_age*l_or), Time = repeat(["2020"], l_ssp*l_age*l_or))
stock_by_age = vcat(stock_by_age, sba_swz_o)
sba_twn_o = DataFrame(scen = repeat(ssps, inner = l_age*l_or), origin = repeat(["TWN"], l_ssp*l_age*l_or), destination = repeat(co_un, outer = l_ssp, inner = l_age), agest = repeat(unique(migflow_all[!,:agest]), outer = l_ssp*l_or), stock_by_age = repeat([0], l_ssp*l_age*l_or), Time = repeat(["2020"], l_ssp*l_age*l_or))
stock_by_age = vcat(stock_by_age, sba_twn_o)
sba_swz_d = DataFrame(scen = repeat(ssps, inner = l_age*l_or), origin = repeat(co_un, outer = l_ssp, inner = l_age), destination = repeat(["SWZ"], l_ssp*l_age*l_or), agest = repeat(unique(migflow_all[!,:agest]), outer = l_ssp*l_or), stock_by_age = repeat([0], l_ssp*l_age*l_or), Time = repeat(["2020"], l_ssp*l_age*l_or))
ind_u = []
for i in eachindex(sba_swz_d[:,1])
    if (sba_swz_d[i,:origin] == "SWZ" && sba_swz_d[i,:destination] == "SWZ") || (sba_swz_d[i,:origin] == "TWN" && sba_swz_d[i,:destination] == "SWZ")
        append!(ind_u, i)
    end
end
delete!(sba_swz_d, ind_u)       # delete doubled rows
stock_by_age = vcat(stock_by_age, sba_swz_d)
sba_twn_d = DataFrame(scen = repeat(ssps, inner = l_age*l_or), origin = repeat(co_un, outer = l_ssp, inner = l_age), destination = repeat(["TWN"], l_ssp*l_age*l_or), agest = repeat(unique(migflow_all[!,:agest]), outer = l_ssp*l_or), stock_by_age = repeat([0], l_ssp*l_age*l_or), Time = repeat(["2020"], l_ssp*l_age*l_or))
ind_u = []
for i in eachindex(sba_twn_d[:,1])
    if (sba_twn_d[i,:origin] == "SWZ" && sba_twn_d[i,:destination] == "TWN") || (sba_twn_d[i,:origin] == "TWN" && sba_twn_d[i,:destination] == "TWN")
        append!(ind_u, i)
    end
end
delete!(sba_twn_d, ind_u)       # delete doubled rows
stock_by_age = vcat(stock_by_age, sba_twn_d)

sort!(migflow_all, [:Time, :scen, :origin, :destination, :agest])
sort!(stock_by_age, [:scen, :origin, :destination, :agest])
migflow_all[!,:stock_by_age] = repeat(stock_by_age[!,:stock_by_age], l_period)
sort!(migflow_all, [:scen, :origin, :destination, :Time, :agest])       # IMPORTANT: need to be ordered by [:scen, :origin, :destination, :Time, :agest] for next step

# Pursue calculation of stock_by_age over time 
for i in 0:l_or*l_dest*l_ssp-1
    for j in (1+l_age):l_period*l_age
        migflow_all[i*l_period*l_age + j,:stock_by_age] = (migflow_all[i*l_period*l_age + j,:duration_prep] - 2020 < 0) ? migflow_all[i*l_period*l_age + (rem(j,l_age) == 0 ? l_age : rem(j,l_age)),:stock_by_age] : 0
    end
end

# Summing initial and dynamic migrants stocks
migflow_all[!,:stock_total] = migflow_all[!,:stock_dyn] .+ migflow_all[!,:stock_by_age]

# Summing by age and finalizing migrants stock data for remittances
migrantstock = combine(groupby(migflow_all, [:scen, :origin, :destination, :Time]), d -> sum(d.stock_total))
rename!(migrantstock, :x1 => :migrantstock)
sort!(migrantstock, [:scen, :origin, :destination, :Time])


################################### Calculating remittances sent by migrants to their origin communities ###############################
# Calculate share of income sent as remittance (remshare)
remittances = migflow[:,[:scen, :origin, :destination, :Time, :ypc_orig, :ypc_dest]]
remittances = innerjoin(remittances, migrantstock, on = [:scen, :origin, :destination, :Time])

rho = CSV.File(joinpath(@__DIR__, "../data/rem_wb/rho.csv")) |> DataFrame
phi = CSV.File(joinpath(@__DIR__,"../data/rem_wb/phi.csv")) |> DataFrame
remittances = leftjoin(remittances, rho, on = [:origin, :destination])
remittances = leftjoin(remittances, phi, on = [:origin, :destination])
# For corridors with no remittance data (with Taiwan and Swaziland), we attribute 0
remittances[!,:rho] = Missings.coalesce.(remittances[!,:rho], 0)
# For corridors with no cost data, we assume that the cost of sending remittances is the mean of all available corridors
for i in eachindex(remittances[:,1])
    if ismissing(remittances[i,:phi])
        remittances[i,:phi] = (remittances[i,:origin] == remittances[i,:destination] ? 0.0 : mean(phi[!,:phi]))
    end
end

# Calculate remittances
# We assume that the remittances share of income stays constant over the lifetime of the migrant at remshare
remittances[!,:remittances] = remittances[!,:migrantstock] .* remittances[!,:ypc_dest] .* remittances[!,:rho] .* (1 .- remittances[!,:phi])

received = combine(groupby(remittances, [:scen, :origin, :Time]), d -> sum(d.remittances))
rename!(received, :x1 => :received, :origin => :country)
sent = combine(groupby(remittances, [:scen, :destination, :Time]), d -> sum(d.remittances))
rename!(sent, :x1 => :sent, :destination => :country)
remittances_country = innerjoin(received, sent, on = [:scen, :country, :Time])
remittances_country[!,:rem] = (remittances_country[!,:received] .- remittances_country[!,:sent]) ./ 10^9

# We account for remittances at t-1 for gdp_mig at t
remittances_country[:,:Time] .+= 5


####################################### Calculate new GDP projections ###########################################
sspall = leftjoin(sspall, remittances_country[:,[:scen, :country, :Time, :rem]], on = [:scen, :country, :Time])

# GDP without migration rescaled through changes in population size and remittances
sspall[!,:gdp_nomig] = (sspall[!,:gdp_mig] .- sspall[!,:rem]) .* sspall[!,:pop_nomig] ./ sspall[!,:pop_mig]

# For 2020, assign same value of GDP to mig and nomig
for i in eachindex(sspall[:,1])
    if ismissing(sspall[i,:gdp_nomig])
        sspall[i,:gdp_nomig] = sspall[i,:gdp_mig]
    end
end

sspall[!,:ypc_mig] = sspall[!,:gdp_mig] .* 10^9 ./ (sspall[!,:pop_mig] .* 10^3)
sspall[!,:ypc_nomig] = sspall[!,:gdp_nomig] .* 10^9 ./ (sspall[!,:pop_nomig] .* 10^3)


####################################### Correct results for specific countries ###############################################
# For 5 countries, |rem| > GDP at some periods: Sri Lanka, Netherlands, Serbia, Jordan, Tajikistan
sspall[!,:remgdpshare] = sspall[!,:rem] ./ sspall[!,:gdp_mig]
# We force |rem| < GDP by keeping constant |rem|/GDP ratio once issues arise
sort!(sspall, [:scen, :country, :Time])
for i in eachindex(sspall[:,1])
    if !ismissing(sspall[i,:remgdpshare]) && (sspall[i,:remgdpshare] > 1 || sspall[i,:remgdpshare] < -1)
        ind = findlast(x -> (!ismissing(x) && x > -1 && x < 1), sspall[1:i,:remgdpshare])
        sspall[i,:gdp_nomig] = sspall[i,:gdp_mig] * (1 - sspall[ind,:remgdpshare]) * sspall[i,:pop_nomig] / (sspall[i,:pop_mig] + sspall[i,:diffmig])
    end
end

# We recompute gdp per capita
sspall[!,:ypc_mig] = sspall[!,:gdp_mig] .* 10^9 ./ (sspall[!,:pop_mig] .* 10^3)
sspall[!,:ypc_nomig] = sspall[!,:gdp_nomig] .* 10^9 ./ (sspall[!,:pop_nomig] .* 10^3)


####################################### Plot results for world #########################################################
gdp_world = combine(groupby(sspall, [:Time, :scen]), d -> (gdp_mig_sum = sum(d.gdp_mig), gdp_nomig_sum = sum(d.gdp_nomig), pop_mig_sum = sum(d.pop_mig), pop_nomig_sum = sum(d.pop_nomig), immigrants_sum = sum(d.inmigsum)))
gdp_world[!,:immig_share] = gdp_world[!,:immigrants_sum] ./ gdp_world[!,:pop_mig_sum] .* 100
gdp_world[!,:ypc_mig] = gdp_world[!,:gdp_mig_sum] .* 10^9 ./ (gdp_world[!,:pop_mig_sum] .* 10^3)
gdp_world[!,:ypc_nomig] = gdp_world[!,:gdp_nomig_sum] .* 10^9 ./ (gdp_world[!,:pop_nomig_sum] .* 10^3)

gdp_world |> @filter(_.Time == 2020 || _.Time == 2050 || _.Time == 2075 || _.Time == 2095) |> @vlplot(
    width=130, height=250,
    mark = :bar, 
    column = {"Time:o", axis={labelFontSize=16, titleFontSize=20}}, 
    x = {"scen:n", scale={rangeStep=8}, axis={labelFontSize=16}, title=nothing},
    y = {"immig_share:q", title="Immigrants in World Population, %", axis={labelFontSize=16, titleFontSize=16}}, 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}},
    spacing=10,
    config={view={stroke=:transparent},axis={domainWidth=1}}
) |> save(joinpath(@__DIR__, "../results/pop_all/", "immigshare_world_bar_6_update.png"))

gdp_w = stack(gdp_world, [:gdp_mig_sum, :gdp_nomig_sum], [:scen, :Time])
rename!(gdp_w, :variable => :gdp_type, :value => :gdp)
gdp_w |> @vlplot(
    width=300, height=250,
    mark={:point, size=60}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"gdp:q", title=nothing, axis={labelFontSize=16}}, 
    title = "World GDP, in billion USD2005", 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    shape = {"gdp_type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/gdp_all/", "gdp_world_6_update.png"))

ypc_w = stack(gdp_world, [:ypc_mig, :ypc_nomig], [:scen, :Time])
rename!(ypc_w, :variable => :ypc_type, :value => :ypc)
ypc_w |> @vlplot(
    width=300, height=250,
    mark={:point, size=60}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"ypc:q", title="World GDP/cap, USD2005/cap/yr", axis={labelFontSize=16,titleFontSize=16}}, 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    shape = {"ypc_type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/ypc_all/", "ypc_world_6_update.png"))

gdp_world[!,:reldif] = gdp_world[!,:gdp_nomig_sum] ./ gdp_world[!,:gdp_mig_sum] .- 1
gdp_world[!,:reldif_ypc] = gdp_world[!,:ypc_mig] ./ gdp_world[!,:ypc_nomig] .- 1

# Plot relative differences in ypc levels without migration compared to with
gdp_world |> @vlplot(
    width=300, height=250,
    mark = {:point, filled=true, size=80}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"reldif_ypc:q", title="Relative change with migration", axis={labelFontSize=16,titleFontSize=16}}, 
    #title = "Relative changes with vs without migration", 
    color = {"scen:n", scale={scheme=:category10}, legend=nothing}
) |> save(joinpath(@__DIR__, "../results/ypc_rel/", "ypcrel_world_6_update.png"))


####################################### Gather some results as illustration #######################################################
# Change in population when zero migration
sspall[!,:popdiff] = sspall[!,:pop_mig] ./ sspall[!,:pop_nomig] .- 1
popdiff = combine(groupby(sspall, :country), d -> popdiff = maximum(abs.(d.popdiff)))
rename!(popdiff, :x1 => :popdiff)

pop_all = stack(sspall, [:pop_mig, :pop_nomig], [:scen, :country, :Time])
rename!(pop_all, :variable => :pop_type, :value => :pop)
# Plot Mexico as an example
pop_all |> @filter(_.country == "MEX") |> @vlplot(
    width=300, height=250,
    mark={:point, size=60}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"pop:q", title=nothing, axis={labelFontSize=16}}, 
    title = "MEX, Population, in thousands", 
    color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    shape = {"pop_type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/pop_all/", "pop_MEX_6_update.png"))

# Remittance flows over time
remdiff = combine(groupby(remittances_country, [:scen, :country]), d -> (minrem = minimum(d.rem), maxrem = maximum(d.rem)))
remittances_country = innerjoin(remittances_country, remdiff, on = [:scen, :country])

# Change in gdp when zero migration
sspall[!,:gdpdiff] = sspall[!,:gdp_mig] ./ sspall[!,:gdp_nomig] .- 1

# Change in gdp per capita when zero migration
sspall[!,:ypcdiff] = (sspall[!,:ypc_mig] .- sspall[!,:ypc_nomig]) ./ abs.(sspall[!,:ypc_nomig])


###################################### Plot geographical maps #####################################
world110m = dataset("world-110m")

countries = unique(sspall[!,:country])
country_code = DataFrame(country = countries, code = unique(sspall[!,:region]))

# Change in population when zero migration
for s in ssps
    @vlplot(width=800, height=600) + @vlplot(mark={:geoshape, stroke = :lightgray}, 
        data={values=world110m, format={type=:topojson, feature=:countries}}, 
        transform = [{lookup=:id, from={data=filter(row -> row[:scen] == s && row[:Time] == 2100, sspall), key=:region, fields=["popdiff"]}}],
        projection={type=:naturalEarth1}, title = {text=string("Relative changes by 2100 with vs without migration, ", s),fontSize=20}, 
        color = {"popdiff:q", scale={domain=[-1,1], scheme=:redblue}, legend={title=nothing, symbolSize=40, labelFontSize=16}}
    ) |> save(joinpath(@__DIR__, "../results/world_maps/", string("popdiff_", s, "_6_update.png")))
end


# Change in gdp when zero migration
for s in ssps
    @vlplot(width=800, height=600) + @vlplot(mark={:geoshape, stroke = :lightgray}, 
        data={values=world110m, format={type=:topojson, feature=:countries}}, 
        transform = [{lookup=:id, from={data=filter(row -> row[:scen] == s && row[:Time] == 2100, sspall), key=:region, fields=["gdpdiff"]}}],
        projection={type=:naturalEarth1}, 
        color = {"gdpdiff:q", scale={domain=[-1,1], scheme=:redblue}, legend={title="Relative change", symbolSize=40, labelFontSize=16}}
    ) |> save(joinpath(@__DIR__, "../results/world_maps/", string("gdpdiff_corr_", s, "_6_update.png")))
end

# Change in gdp per capita when zero migration
for s in ssps
    @vlplot(width=800, height=600) + @vlplot(mark={:geoshape, stroke = :lightgray}, 
        data={values=world110m, format={type=:topojson, feature=:countries}}, 
        transform = [{lookup=:id, from={data=filter(row -> row[:scen] == s && row[:Time] == 2100, sspall), key=:region, fields=["ypcdiff"]}}],
        projection={type=:naturalEarth1}, #title = {text=string("Relative changes in GDP/cap by 2100 with vs without migration, ", s),fontSize=20}, 
        color = {"ypcdiff:q", scale={domain=[-0.19,0.19], scheme=:redblue}, legend={title="Relative change", symbolSize=40, labelFontSize=16}}
    ) |> save(joinpath(@__DIR__, "../results/world_maps/", string("ypcdiff_corr_", s, "_6_update.png")))
end


###################################### Plot results at continent level ############################
# Define continents based on KC and Lutz (2017): 
# Northern America (NOA), Latin America and the Caribbean (LAC), Europe (EUR), Africa (AFR), Asia (ASIA), Oceania (OCE)
regions = DataFrame(
    country = sort!(countries), 
    worldregion = [
        "LAC","ASIA","AFR","EUR","ASIA","LAC","EUR","OCE","EUR","ASIA","AFR","EUR","AFR","AFR","ASIA","EUR","ASIA","LAC","EUR","EUR","LAC","LAC","LAC",
        "LAC","ASIA","ASIA","AFR","AFR","NOA","EUR","LAC","ASIA","AFR","AFR","AFR","AFR","LAC","AFR","AFR","LAC","LAC","EUR","EUR","EUR","AFR","EUR",
        "LAC","AFR","LAC","AFR","AFR","EUR","EUR","AFR","EUR","OCE","EUR","AFR","EUR","EUR","AFR","AFR","AFR","AFR","AFR","EUR","LAC","LAC","ASIA",
        "LAC","EUR","LAC","EUR","ASIA","ASIA","EUR","ASIA","ASIA","EUR","ASIA","EUR","LAC","ASIA","ASIA","ASIA","AFR","ASIA","ASIA","ASIA","ASIA","ASIA","ASIA",
        "AFR","AFR","LAC","ASIA","AFR","EUR","EUR","EUR","ASIA","AFR","EUR","AFR","ASIA","LAC","EUR","AFR","EUR","ASIA","EUR","ASIA","AFR","AFR","AFR",
        "AFR","ASIA","AFR","OCE","AFR","AFR","LAC","EUR","EUR","ASIA","OCE","ASIA","ASIA","LAC","LAC","ASIA","OCE","EUR","LAC","EUR","LAC","ASIA","OCE",
        "ASIA","EUR","ASIA","AFR","ASIA","AFR","AFR","ASIA","OCE","AFR","LAC","AFR","EUR","AFR","LAC","EUR","EUR","EUR","AFR","ASIA","AFR","AFR","ASIA",
        "ASIA","ASIA","ASIA","OCE","LAC","AFR","ASIA","ASIA","AFR","AFR","EUR","LAC","NOA","ASIA","LAC","LAC","ASIA","OCE","OCE","ASIA","AFR","AFR","AFR"
    ]
)

# Plot ypc levels per world region with and without migration
sspall = innerjoin(sspall, regions, on = :country)
sspall_reg = combine(groupby(sspall, [:scen, :Time, :worldregion]), d -> (pop_mig_reg = sum(d.pop_mig), pop_nomig_reg = sum(d.pop_nomig), gdp_mig_reg = sum(d.gdp_mig), gdp_nomig_reg = sum(d.gdp_nomig)))

sspall_reg[!,:ypc_mig_reg] = sspall_reg[!,:gdp_mig_reg] .* 10^9 ./ (sspall_reg[!,:pop_mig_reg] .* 10^3)
sspall_reg[!,:ypc_nomig_reg] = sspall_reg[!,:gdp_nomig_reg] .* 10^9 ./ (sspall_reg[!,:pop_nomig_reg] .* 10^3)
ypc_reg = stack(sspall_reg, [:ypc_mig_reg, :ypc_nomig_reg], [:scen, :worldregion, :Time])
rename!(ypc_reg, :variable => :ypc_type, :value => :ypc)
for r in unique(regions[!,:worldregion])
    ypc_reg |> @filter(_.worldregion == r) |> @vlplot(
        width=300, height=250,
        mark={:point, size=60}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"ypc:q", title="GDP/cap, USD2005/cap/yr", axis={labelFontSize=16,titleFontSize=16}}, 
        title = string(r), 
        color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
        shape = {"ypc_type:o", scale={range=["circle","triangle-up"]}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}
    ) |> save(joinpath(@__DIR__, "../results/ypc_reg/", string("ypc_", r, "_6_update.png")))
end

# Plot relative differences in ypc levels without migration compared to with
sspall_reg[!,:popdif_reg] = sspall_reg[!,:pop_mig_reg] ./ sspall_reg[!,:pop_nomig_reg] .- 1
sspall_reg[!,:reldif_reg] = sspall_reg[!,:ypc_mig_reg] ./ sspall_reg[!,:ypc_nomig_reg] .- 1
for r in unique(regions[!,:worldregion])
    sspall_reg |> @filter(_.worldregion == r) |> @vlplot(
        width=300, height=250,
        mark = {:point, filled=true, size=80}, x = {"Time:o", axis={labelFontSize=16}, title=nothing}, y = {"reldif_reg:q", title="Relative change with migration", axis={labelFontSize=16,titleFontSize=16}}, 
        title = string(r), 
        color = {"scen:n", scale={scheme=:category10}, legend=nothing}
    ) |> save(joinpath(@__DIR__, "../results/ypc_reg/", string("ypcrel_", r, "_6_update.png")))
end

# Plot heat tables and geographical maps of migrant and remittances flows in 2050 and end of century
rename!(regions, :country => :origin, :worldregion => :or_reg)
migflow = innerjoin(migflow, regions, on = :origin)
remittances = innerjoin(remittances, regions, on = :origin)
rename!(regions, :origin => :destination, :or_reg => :dest_reg)
migflow = innerjoin(migflow, regions, on = :destination)
remittances = innerjoin(remittances, regions, on = :destination)
migflow_reg = combine(groupby(migflow, [:scen, :Time, :or_reg, :dest_reg]), d -> (move_reg = sum(d.move_inmigbase), share_or_reg = sum(d.move_inmigbase) / sum(d.pop_orig), share_dest_reg = sum(d.move_inmigbase) / sum(d.pop_dest)))
remittances_reg = combine(groupby(remittances, [:scen, :Time, :or_reg, :dest_reg]), d -> (rem_reg = sum(d.remittances)))

rename!(remittances_reg, :x1 => :rem_reg)
rename!(sspall_reg, :worldregion => :or_reg, :gdp_mig_reg => :gdp_or_reg)
remittances_reg = innerjoin(remittances_reg, sspall_reg[:,[:scen, :Time, :or_reg, :gdp_or_reg]], on = [:scen, :Time, :or_reg])
rename!(sspall_reg, :or_reg => :dest_reg, :gdp_or_reg => :gdp_dest_reg)
remittances_reg = innerjoin(remittances_reg, sspall_reg[:,[:scen, :Time, :dest_reg, :gdp_dest_reg]], on = [:scen, :Time, :dest_reg])
remittances_reg[!,:remshare_or_reg] = remittances_reg[!,:rem_reg] ./ (remittances_reg[!,:gdp_or_reg] .* 10^9)
remittances_reg[!,:remshare_dest_reg] = remittances_reg[!,:rem_reg] ./ (remittances_reg[!,:gdp_dest_reg] .* 10^9)

for d in [2050, 2095]
    for s in ssps
        migflow_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"move_reg:q", scale={domain=[0,10^6], scheme=:goldred}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/pop_reg/", string("migflow_reg_", d,"_", s, "_6_update.png")))
        migflow_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"share_or_reg:q", scale={domain=[0,0.0001], scheme=:goldred}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/pop_reg/", string("migflow_reg_share_or_", d,"_", s, "_6_update.png")))
        migflow_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"share_dest_reg:q", scale={domain=[0,0.0001], scheme=:goldred}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/pop_reg/", string("migflow_reg_share_dest_", d,"_", s, "_6_update.png")))
    end
end

for d in [2050, 2095]
    for s in ssps
        remittances_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"rem_reg:q", scale={domain=[0,10^12]}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/ypc_reg/", string("remittances_reg_", d,"_", s, "_6_update.png")))
        remittances_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"remshare_or_reg:q", scale={domain=[0,0.01]}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/ypc_reg/", string("remittances_reg_share_or_", d,"_", s, "_6_update.png")))
        remittances_reg |> @filter(_.Time == d && _.scen == s) |> @vlplot(
            :rect, y=:or_reg, x=:dest_reg, color={"remshare_dest_reg:q", scale={domain=[0,0.01]}},title = string(s, ", ", d)
        ) |> save(joinpath(@__DIR__, "../results/ypc_reg/", string("remittances_reg_share_dest_", d,"_", s, "_6_update.png")))
    end
end

# For illustration, plot world regions used for displaying results
@vlplot(width=800, height=600) + @vlplot(mark={:geoshape, stroke = :lightgray}, 
    data={values=world110m, format={type=:topojson, feature=:countries}
    #, transform = [{typ = filter, expr = (:id in country_code[:code])}]
    }, 
    transform = [{lookup=:id, from={data=filter(row -> row[:scen] == "SSP2" && row[:Time] == 2100, sspall), key=:region, fields=["worldregion"]}}],
    projection={type=:naturalEarth1}, title = {text=string("World regions considered"),fontSize=20}, 
    color = {"worldregion:n", scale={scheme=:set3}, legend={title=nothing, symbolSize=40, labelFontSize=16}}
) |> save(joinpath(@__DIR__, "../results/pop_reg/", string("regions.png")))


####################################### Write output files with results ############################################################################
CSV.write(joinpath(@__DIR__, "../results/sspall_6_update.csv"), sspall)
CSV.write("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/results_large/migflow_6_update.csv", migflow)
CSV.write("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/results_large/remittances_6_update.csv", remittances)

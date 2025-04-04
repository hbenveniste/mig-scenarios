using CSV, DataFrames, DelimitedFiles, FileIO, FilePaths, XLSX, Statistics
using FixedEffectModels, RegressionTables
using GLM
using Query
using Distances


########################## Prepare migration flows data ########################################
# Using Abel (2018). Note: too large to be stored on Github; available from https://guyabel.com/publication/global-migration-estimates-by-gender/ 
migflow_allstockdemo = CSV.File("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Abel_data/gf_imr.csv") |> DataFrame
# Using Abel and Cohen (2019)
migflow_alldata = CSV.File(joinpath(@__DIR__, "../data/migflow_all/ac19.csv")) |> DataFrame

# From Abel 2018, we choose demographic data from the UN's WPP2015, and migrant stock data from the World Bank for 1960-1990 and the UN for 1990-2015.
migflow = @from i in migflow_allstockdemo begin
    @where i.demo == "wpp2015" && i.sex == "b" && ((i.stock == "wb11" && i.year0 < 1990) || (i.stock == "un15" && i.year0 >= 1990))
    @select {i.stock, i.demo, i.sex, i.year0, i.interval, i.orig, i.dest, i.orig_code, i.dest_code, i.flow}
    @collect DataFrame
end

# From Abel and Cohen 2019, we choose Azose and Raftery's data for 1990-2015, based a demographic accounting, pseudo-Bayesian method, which performs the best against most correlation measures with actual data
migflow_ar = migflow_alldata[:,[:year0, :orig, :dest, :da_pb_closed]]


########################## Prepare population data from the Wittgenstein Centre, based on historical data from the WPP 2019 ##################################
pop_allvariants = CSV.File("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Pop_hist_data/WPP2019.csv") |> DataFrame
# We use the Medium variant, the most commonly used. Unit: thousands
pop = @from i in pop_allvariants begin
    @where i.Variant == "Medium" && i.Time < 2016 
    @select {i.LocID, i.Location, i.Time, i.PopTotal}
    @collect DataFrame
end


########################## Prepare GDP data from the World Bank's WDI as available at the IIASA SSP database #####################
# Unit: billion US$ 2005 / year PPP
gdp_unstacked = XLSX.openxlsx(joinpath(@__DIR__, "../data/gdphist_wdisspdb/gdphist.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["data"];first_row=2)...)
end
select!(gdp_unstacked, Not([:Model, Symbol("Scenario (History)"), :Variable, :Unit]))
gdp = stack(gdp_unstacked, 2:size(gdp_unstacked, 2))
rename!(gdp, :variable => :year0, :value => :gdp)
gdp[!,:year0] = map( x -> parse(Int, String(x)), gdp[!,:year0])


########################## Prepare distance data based on UN POP data on capital cities coordinates ##############################
loc = load(joinpath(@__DIR__, "../data/distance_unpop/WUP2014-F13-Capital_Cities.xls"), "DATA!A17:J257") |> DataFrame
# choose one capital per country
inddupl = []
for c in ["Yamoussoukro", "Pretoria", "Bloemfontein", "Porto-Novo", "Sucre", "Sri Jayewardenepura Kotte", "s-Gravenhage (The Hague)", "St. Peter Port"]
    ind = findfirst(x -> x == c, loc[!,Symbol("Capital City")])
    append!(inddupl, ind)
end
sort!(inddupl)
delete!(loc, inddupl)
select!(loc, Not([:Index, Symbol("Capital City"), :Note, Symbol("Capital Type"), Symbol("City code"), Symbol("Population (thousands)")]))
loc[!,Symbol("Country code")] = map( x -> trunc(Int, x), loc[!,Symbol("Country code")])
rename!(loc, Symbol("Country code") => :country_code, Symbol("Country or area") => :country)
push!(loc, ["Taiwan", 158, 25.0330, 121.5654])

earthradius = 6372.8        # in km

# Calculating distances between countries as distances between their capital cities, using the Haversine formula of the Distances package
dist = DataFrame(orig_code = repeat(loc[!,:country_code], inner = size(loc,1)), lat_or = repeat(loc[!,:Latitude], inner = size(loc,1)), lon_or = repeat(loc[!,:Longitude], inner = size(loc,1)), dest_code = repeat(loc[!,:country_code], outer = size(loc,1)), lat_dest = repeat(loc[!,:Latitude], outer = size(loc,1)), lon_dest = repeat(loc[!,:Longitude], outer = size(loc,1)))
dist[!,:loc_or] = [tuple(dist[!,:lat_or][i], dist[!,:lon_or][i]) for i in 1:size(dist,1)]
dist[!,:loc_dest] = [tuple(dist[!,:lat_dest][i], dist[!,:lon_dest][i]) for i in 1:size(dist,1)]
dist[!,:distance] = [haversine(dist[!,:loc_or][i], dist[!,:loc_dest][i], earthradius) for i in 1:size(dist, 1)]

# Add country codes
iso3c_isonum = CSV.File(joinpath(@__DIR__,"../data/iso3c_isonum.csv")) |> DataFrame
rename!(iso3c_isonum, :isonum => :orig_code)
distances = innerjoin(dist, iso3c_isonum, on = :orig_code)
rename!(distances, :iso3c => :orig)
rename!(iso3c_isonum, :orig_code => :dest_code)
distances = innerjoin(distances, iso3c_isonum, on = :dest_code)
rename!(distances, :iso3c => :dest)
select!(distances, Not(1:8))


########################################### Prepare remittance data based on World Bank data ##############################
rho = CSV.File(joinpath(@__DIR__, "../data/rem_wb/rho.csv")) |> DataFrame
phi = CSV.File(joinpath(@__DIR__,"../data/rem_wb/phi.csv")) |> DataFrame
remittances = leftjoin(rho, phi, on = [:origin, :destination])
# For corridors with no cost data, we assume that the cost of sending remittances is the mean of all available corridors
for i in eachindex(remittances[:,1])
    if ismissing(remittances[i,:phi])
        remittances[i,:phi] = (remittances[i,:origin] == remittances[i,:destination] ? 0.0 : mean(phi[!,:phi]))
    end
end


############################################ Prepare land surface data from WDI (World Bank) #######################
land = XLSX.openxlsx(joinpath(@__DIR__, "../data/landarea_wb/land_area.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["Data"])...)
end
rename!(land, Symbol("Country Code") => :country, Symbol("2017 [YR2017]") => :area)
select!(land, [:country, :area])

# Deal with Sudan and South Sudan: Sudan surface is 1,886,068 km^2
ind_ssd = findfirst(land[!,:country].=="SSD")
delete!(land,ind_ssd)
ind_sdn = findfirst(land[!,:country].=="SDN")
land[ind_sdn,:area] = 1886068 


############################################# Prepare data on common official languages #############################
countries = unique(migflow_ar[:,:orig])
comol = DataFrame(
    orig = repeat(countries, inner = length(countries)), 
    dest = repeat(countries, outer=length(countries)), 
    comofflang = zeros(length(countries)^2)
)

# Assign common official languages
offlang = Dict(
    "Chinese" => ["CHN", "TWN","HKG","SGP","MAC"],
    "Spanish" => ["MEX", "COL","ARG","ESP","VEN","PER","CHL","ECU","CUB","GTM","DOM","HND","BOL","SLV","NIC","CRI","PRY","URY","PAN","PRI","GNQ"],
    "English" => ["USA","GBR","CAN","AUS","ZAF","IRL","GHA","NZL","SGP","LBR","PAN","ZWE","ZMB","HKG","JEY","VGB","GUM","BMU","CYM","BWA","LCA","GIB","MLT","BLZ","MNP","VUT","SYC","COK","FLK","PLW","ASM","WSM","NFK","NIU","TKL","NRU","ATG","BHS","BRB","BDI","CMR","DMA","SWZ","FJI","GMB","GRD","GUY","IND","JAM","KEN","KIR","LSO","MWI","MHL","FSM","NAM","NGA","PAK","PNG","PHL","RWA","KNA","VCT","SLE","SLB","SSD","SDN","TZA","TON","TTO","TUV","UGA"],
    "Arabic" => ["EGY","DZA","SAU","IRQ","YEM","MAR","SDN","SYR","JOR","TUN","LBY","LBN","ARE","OMN","KWT","TCD","ISR","QAT","BHR","ESH","DJI","COM","ERI","MRT","PSE","SOM","TZA"],
    "Portuguese" => ["BRA","AGO","PRT","MOZ","GNB","STP","CPV","MAC","GNQ","TLS"],
    "Russian" => ["RUS","BLR","KGZ","KAZ"],
    "German" => ["DEU","AUT","CHE","BEL","LIE","LUX"],
    "Tamil" => ["SGP","LKA","IND"],
    "French" => ["FRA","CAN","BEL","CHE","MLI","PYF","NCL","MYT","BDI","LUX","MCO","RWA","WLF","VUT","SYC","COD","CMR","MDG","CIV","NER","BFA","SEN","TCD","GIN","BEN","HTI","TGO","CAF","COG","GAB","GNQ","DJI"],
    "Korean" => ["PRK","KOR"],
    "Turkish" => ["TUR","CYP"],
    "Italian" => ["ITA","CHE","SMR"],
    "Malay" => ["IDN","MYS","SGP","BRN"],
    "Hindustani" => ["FJI","IND","PAK"],
    "Sotho" => ["ZAF","LSO","ZWE"],
    "Quechua" => ["PER","BOL","ECU"],
    "Persian" => ["IRN","AFG","TJK"],
    "Dutch" => ["NLD","BEL","SUR"],
    "Yoruba" => ["NGA","BEN","TGO","BRA"],
    "Swahili" => ["TZA","KEN","UGA","RWA"],
    "Hausa" => ["NGA","NER","CMR","GHA","BEN","CIV","TGO","SDN"],
    "Aymara" => ["PER","BOL"],
    "Bengali" => ["BGD","IND"],
    "Berber" => ["DZA","MAR"],
    "Greek" => ["GRC","CYP"],
    "Guarani" => ["BOL","PRY"],
    "Romanian" => ["ROU","MDA"],
    "Rundi" => ["BDI","RWA"],
    "Swati" => ["ZAF","SWZ"],
    "Swedish" => ["SWE","FIN"],
    "Tswana" => ["ZAF","BWA"]
)

for i in eachindex(comol[:,1]) ; if comol[i,:orig] == comol[i,:dest] ; comol[i,:comofflang] = 1 end end
for l in offlang
    for c1 in l[2]
        for c2 in l[2]
            ind = intersect(findall(comol[!,:orig].==c1),findall(comol[!,:dest].==c2))
            if !isempty(ind)
                comol[ind[1],:comofflang] = 1
            end
        end
    end
end


########################## Join datasets #########################################
# Joining the datasets
data = innerjoin(migflow, distances, on = [:orig, :dest])
select!(data, Not([:stock, :demo, :sex]))
data[!,:flow] = float(data[!,:flow])
rename!(data, :flow => :flow_Abel)
data_ar = innerjoin(migflow_ar, distances, on = [:orig, :dest])
rename!(data_ar, :da_pb_closed => :flow_AzoseRaftery)

rename!(pop, :LocID => :orig_code, :Time => :year0)
data = innerjoin(data, pop, on = [:year0, :orig_code])
rename!(data, :Location => :orig_name, :PopTotal => :pop_orig)
rename!(pop, :orig_code => :dest_code)
data = innerjoin(data, pop, on = [:year0, :dest_code])
rename!(data, :Location => :dest_name, :PopTotal => :pop_dest)
iso3c_isonum = CSV.File(joinpath(@__DIR__, "../data/iso3c_isonum.csv")) |> DataFrame
rename!(pop, :dest_code => :isonum)
pop = innerjoin(pop, iso3c_isonum, on = :isonum)
rename!(pop, :iso3c => :orig, :PopTotal => :pop_orig)
data_ar = innerjoin(data_ar, pop[:,[:year0, :pop_orig, :orig]], on = [:year0, :orig])
rename!(pop, :orig => :dest, :pop_orig => :pop_dest)
data_ar = innerjoin(data_ar, pop[:,[:year0, :pop_dest, :dest]], on = [:year0, :dest])

rename!(gdp, :Region => :orig)
data = innerjoin(data, gdp, on = [:year0, :orig])
data_ar = innerjoin(data_ar, gdp, on = [:year0, :orig])
rename!(data, :gdp => :gdp_orig)
rename!(data_ar, :gdp => :gdp_orig)
rename!(gdp, :orig => :dest)
data = innerjoin(data, gdp, on = [:year0, :dest])
data_ar = innerjoin(data_ar, gdp, on = [:year0, :dest])
rename!(data, :gdp => :gdp_dest)
rename!(data_ar, :gdp => :gdp_dest)

rename!(remittances, :origin => :orig, :destination => :dest, :rho => :remshare, :phi => :remcost)
data = innerjoin(data, remittances, on = [:orig,:dest])
data_ar = innerjoin(data_ar, remittances, on =[:orig, :dest])

rename!(land, :country => :orig, :area => :area_orig)
data = innerjoin(data, land, on = :orig)
data_ar = innerjoin(data_ar, land, on = :orig)
rename!(land, :orig => :dest, :area_orig => :area_dest)
data = innerjoin(data, land, on = :dest)
data_ar = innerjoin(data_ar, land, on = :dest)

data = innerjoin(data, comol, on = [:orig,:dest])
data_ar = innerjoin(data_ar, comol, on = [:orig,:dest])

# Making units consistent 
# Flows are for a multiple-year period. We compute annual migrant flows as constant over said period
data[!,:flow_Abel] ./= data[!,:interval] 
data[!,:pop_orig] .*= 1000        # pop is in thousands
data[!,:pop_dest] .*= 1000
data[!,:gdp_orig] .*= 10^9        # gdp is in billion $
data[!,:gdp_dest] .*= 10^9
data_ar[!,:flow_AzoseRaftery] ./= 5
data_ar[!,:pop_orig] .*= 1000        # pop is in thousands
data_ar[!,:pop_dest] .*= 1000
data_ar[!,:gdp_orig] .*= 10^9        # gdp is in billion $
data_ar[!,:gdp_dest] .*= 10^9

# Creating gdp per capita variables
data[!,:ypc_orig] = data[!,:gdp_orig] ./ data[!,:pop_orig]
data[!,:ypc_dest] = data[!,:gdp_dest] ./ data[!,:pop_dest]
data[!,:ypc_ratio] = data[!,:ypc_dest] ./ data[!,:ypc_orig]
data_ar[!,:ypc_orig] = data_ar[!,:gdp_orig] ./ data_ar[!,:pop_orig]
data_ar[!,:ypc_dest] = data_ar[!,:gdp_dest] ./ data_ar[!,:pop_dest]
data_ar[!,:ypc_ratio] = data_ar[!,:ypc_dest] ./ data_ar[!,:ypc_orig]

# Create ratios of move_od / stay_o variables
emig = combine(groupby(data, [:year0,:orig]), d->sum(d.flow_Abel))
rename!(emig, :x1 => :emig)
data = innerjoin(data, emig, on=[:year0,:orig])
data[!,:stay_orig] = data[!,:pop_orig] .- data[!,:emig]
data[!,:mig_ratio] = data[!,:flow_Abel] ./ data[!,:stay_orig]
data[!,:mig_ratio_tot] = data[!,:flow_Abel] ./ data[!,:pop_orig]
emig_ar = combine(groupby(data_ar, [:year0,:orig]), d->sum(d.flow_AzoseRaftery))
rename!(emig_ar, :x1 => :emig)
data_ar = innerjoin(data_ar, emig_ar, on=[:year0,:orig])
data_ar[!,:stay_orig] = data_ar[!,:pop_orig] .- data_ar[!,:emig]
data_ar[!,:mig_ratio] = data_ar[!,:flow_AzoseRaftery] ./ data_ar[!,:stay_orig]
data_ar[!,:mig_ratio_tot] = data_ar[!,:flow_AzoseRaftery] ./ data_ar[!,:pop_orig]

# Create density variables
data[!,:density_orig] = data[!,:pop_orig] ./ data[!,:area_orig]
data[!,:density_dest] = data[!,:pop_dest] ./ data[!,:area_dest]
data[!,:density_ratio] = data[!,:density_dest] ./ data[!,:density_orig]
data_ar[!,:density_orig] = data_ar[!,:pop_orig] ./ data_ar[!,:area_orig]
data_ar[!,:density_dest] = data_ar[!,:pop_dest] ./ data_ar[!,:area_dest]
data_ar[!,:density_ratio] = data_ar[!,:density_dest] ./ data_ar[!,:density_orig]


############################# Calibrate gravity equation ##################################
# log transformation
logdata = DataFrame(
    year = data[!,:year0], 
    orig = data[!,:orig], 
    dest = data[!,:dest],
    remshare = data[!,:remshare], 
    remcost = data[!,:remcost],
    comofflang = data[!,:comofflang]
)
for name in [:flow_Abel, :mig_ratio, :mig_ratio_tot, :pop_orig, :pop_dest, :area_orig, :area_dest, :density_orig, :density_dest, :density_ratio, :gdp_orig, :gdp_dest, :ypc_orig, :ypc_dest, :ypc_ratio, :distance]
    logdata[!,name] = [log(data[!,name][i]) for i in 1:size(logdata, 1)]
end
logdata_ar = DataFrame(
    year = data_ar[!,:year0], 
    orig = data_ar[!,:orig], 
    dest = data_ar[!,:dest], 
    remshare = data_ar[!,:remshare], 
    remcost = data_ar[!,:remcost],
    comofflang = data_ar[!,:comofflang]
)
for name in [:flow_AzoseRaftery, :mig_ratio, :mig_ratio_tot, :pop_orig, :pop_dest, :area_orig, :area_dest, :density_orig, :density_dest, :density_ratio, :gdp_orig, :gdp_dest, :ypc_orig, :ypc_dest, :ypc_ratio, :distance]
    logdata_ar[!,name] = [log(data_ar[!,name][i]) for i in 1:size(logdata_ar, 1)]
end

# Remove rows with distance = 0 or flow = 0
gravity = @from i in logdata begin
    @where i.distance != -Inf && i.flow_Abel != -Inf && i.mig_ratio_tot != -Inf && i.mig_ratio != -Inf
    @select {i.year, i.orig, i.dest, i.flow_Abel, i.mig_ratio, i.mig_ratio_tot, i.pop_orig, i.pop_dest, i.area_orig, i.area_dest, i.density_orig, i.density_dest, i.density_ratio, i.gdp_orig, i.gdp_dest, i.ypc_orig, i.ypc_dest, i.ypc_ratio, i.distance, i.remshare, i.remcost, i.comofflang}
    @collect DataFrame
end
dropmissing!(gravity)       # remove rows with missing values in ypc_orig or ypc_dest
gravity_ar = @from i in logdata_ar begin
    @where i.distance != -Inf && i.mig_ratio_tot != -Inf && i.mig_ratio != -Inf && i.flow_AzoseRaftery != -Inf
    @select {i.year, i.orig, i.dest, i.flow_AzoseRaftery, i.mig_ratio, i.mig_ratio_tot, i.pop_orig, i.pop_dest, i.area_orig, i.area_dest, i.density_orig, i.density_dest, i.density_ratio, i.gdp_orig, i.gdp_dest, i.ypc_orig, i.ypc_dest, i.ypc_ratio, i.distance, i.remshare, i.remcost, i.comofflang}
    @collect DataFrame
end
dropmissing!(gravity_ar)       # remove rows with missing values in ypc_orig or ypc_dest

# Compute linear regression with FixedEffectModels package
rr1 = reg(gravity_ar, @formula(flow_AzoseRaftery ~ pop_orig + pop_dest + ypc_orig + ypc_dest + distance + remshare + remcost + comofflang + fe(year)), Vcov.cluster(:orig, :dest), save=true)
rr2 = reg(gravity_ar, @formula(flow_AzoseRaftery ~ pop_orig + pop_dest +  ypc_orig + ypc_dest + distance + remshare + remcost + comofflang + fe(orig) + fe(dest) + fe(year)), Vcov.cluster(:orig, :dest), save=true)
rr3 = reg(gravity_ar, @formula(mig_ratio_tot ~ density_ratio + ypc_ratio + distance + remshare + remcost + comofflang + fe(year)), Vcov.cluster(:orig, :dest), save=true)
rr4 = reg(gravity_ar, @formula(mig_ratio_tot ~ density_ratio + ypc_ratio + distance + remshare + remcost + comofflang + fe(orig) + fe(dest) + fe(year)), Vcov.cluster(:orig, :dest), save=true)
rr5 = reg(gravity, @formula(flow_Abel ~ pop_orig + pop_dest + ypc_orig + ypc_dest + distance + remshare + remcost + comofflang + fe(year)), Vcov.cluster(:orig, :dest), save=true)

regtable(rr1, rr2, rr3, rr4, rr5; regression_statistics=[:nobs, :r2, :r2_within])     

beta = DataFrame(
    regtype = ["reg_ar_yfe","reg_ar_odyfe","reg_abel_yfe"],
    b1 = [0.686,0.746,0.575],       # pop_orig
    b2 = [0.681,-0.748,0.603],       # pop_dest
    b4 = [0.428,0.155,0.105],       # ypc_orig
    b5 = [0.838,-0.004,0.791],       # ypc_dest
    b7 = [-1.291,-1.472,-1.038],       # distance
    b8 = [0.155,0.190,0.090],       # remshare
    b9 = [-8.382,-11.706,-14.203],       # remcost
    b10 = [1.727,1.606,1.439]     # comofflang
)
beta_ratio = DataFrame(
    regtype = ["reg_ratio_ar_yfe","reg_ratio_ar_odyfe"],
    b3 = [0.037,-0.249],       # density_ratio
    b6 = [0.145,-0.088],       # ypc_ratio
    b7 = [-1.166,-1.472],       # distance
    b8 = [-0.053,0.190],       # remshare
    b9 = [-8.860,-11.703],       # remcost
    b10 = [1.217,1.607]       # comofflang
)

# Compute constant including year fixed effect as average of beta0 + yearFE
cst_ar_yfe = hcat(gravity_ar[:,Not([:orig,:dest,:year])], fe(rr1))
cst_ar_yfe[!,:constant] = cst_ar_yfe[!,:flow_AzoseRaftery] .- beta[1,:b1] .* cst_ar_yfe[!,:pop_orig] .- beta[1,:b2] .* cst_ar_yfe[!,:pop_dest] .- beta[1,:b4] .* cst_ar_yfe[!,:ypc_orig] .- beta[1,:b5] .* cst_ar_yfe[!,:ypc_dest] .- beta[1,:b7] .* cst_ar_yfe[!,:distance] .- beta[1,:b8] .* cst_ar_yfe[!,:remshare] .- beta[1,:b9] .* cst_ar_yfe[!,:remcost] .- beta[1,:b10] .* cst_ar_yfe[!,:comofflang]
constant_ar_yfe = mean(cst_ar_yfe[!,:constant])

cst_ar_odyfe = hcat(gravity_ar[:,Not([:orig,:dest,:year])], fe(rr2))
cst_ar_odyfe[!,:constant] = cst_ar_odyfe[!,:flow_AzoseRaftery] .- beta[2,:b1] .* cst_ar_odyfe[!,:pop_orig] .- beta[2,:b2] .* cst_ar_odyfe[!,:pop_dest] .- beta[2,:b4] .* cst_ar_odyfe[!,:ypc_orig] .- beta[2,:b5] .* cst_ar_odyfe[!,:ypc_dest] .- beta[2,:b7] .* cst_ar_odyfe[!,:distance] .- beta[2,:b8] .* cst_ar_odyfe[!,:remshare] .- beta[2,:b9] .* cst_ar_odyfe[!,:remcost] .- beta[2,:b10] .* cst_ar_odyfe[!,:comofflang] .- cst_ar_odyfe[!,:fe_orig] .- cst_ar_odyfe[!,:fe_dest]
constant_ar_odyfe = mean(cst_ar_odyfe[!,:constant])

cst_ratio_ar_yfe = hcat(gravity_ar[:,Not([:orig,:dest,:year])], fe(rr3))
cst_ratio_ar_yfe[!,:constant] = cst_ratio_ar_yfe[!,:flow_AzoseRaftery] .- beta_ratio[1,:b3] .* cst_ratio_ar_yfe[!,:density_ratio] .- beta_ratio[1,:b6] .* cst_ratio_ar_yfe[!,:ypc_ratio] .- beta_ratio[1,:b7] .* cst_ratio_ar_yfe[!,:distance] .- beta_ratio[1,:b8] .* cst_ratio_ar_yfe[!,:remshare] .- beta_ratio[1,:b9] .* cst_ratio_ar_yfe[!,:remcost] .- beta_ratio[1,:b10] .* cst_ratio_ar_yfe[!,:comofflang]
constant_ratio_ar_yfe = mean(cst_ratio_ar_yfe[!,:constant])

cst_ratio_ar_odyfe = hcat(gravity_ar[:,Not([:orig,:dest,:year])], fe(rr4))
cst_ratio_ar_odyfe[!,:constant] = cst_ratio_ar_odyfe[!,:flow_AzoseRaftery] .- beta_ratio[2,:b3] .* cst_ratio_ar_odyfe[!,:density_ratio] .- beta_ratio[2,:b6] .* cst_ratio_ar_odyfe[!,:ypc_ratio] .- beta_ratio[2,:b7] .* cst_ratio_ar_odyfe[!,:distance] .- beta_ratio[2,:b8] .* cst_ratio_ar_odyfe[!,:remshare] .- beta_ratio[2,:b9] .* cst_ratio_ar_odyfe[!,:remcost] .- beta_ratio[2,:b10] .* cst_ratio_ar_odyfe[!,:comofflang] .- cst_ratio_ar_odyfe[!,:fe_orig] .- cst_ratio_ar_odyfe[!,:fe_dest]
constant_ratio_ar_odyfe = mean(cst_ratio_ar_odyfe[!,:constant])

cst_abel_yfe = hcat(gravity[:,Not([:orig,:dest,:year])], fe(rr5))
cst_abel_yfe[!,:constant] = cst_abel_yfe[!,:flow_Abel] .- beta[3,:b1] .* cst_abel_yfe[!,:pop_orig] .- beta[3,:b2] .* cst_abel_yfe[!,:pop_dest] .- beta[3,:b4] .* cst_abel_yfe[!,:ypc_orig] .- beta[3,:b5] .* cst_abel_yfe[!,:ypc_dest] .- beta[3,:b7] .* cst_abel_yfe[!,:distance] .- beta[3,:b8] .* cst_abel_yfe[!,:remshare] .- beta[3,:b9] .* cst_abel_yfe[!,:remcost] .- beta[3,:b10] .* cst_abel_yfe[!,:comofflang]
constant_abel_yfe = mean(cst_abel_yfe[!,:constant])

beta[!,:b0] = [constant_ar_yfe,constant_ar_odyfe,constant_abel_yfe] 
beta_ratio[!,:b0] = [constant_ratio_ar_yfe,constant_ratio_ar_odyfe] 

# Gather FE values
fe_ar_yfe = hcat(gravity_ar[:,[:year,:orig,:dest]], fe(rr1))
fe_ar_odyfe = hcat(gravity_ar[:,[:year,:orig,:dest]], fe(rr2))
fe_ratio_ar_yfe = hcat(gravity_ar[:,[:year,:orig,:dest]], fe(rr3))
fe_ratio_ar_odyfe = hcat(gravity_ar[:,[:year,:orig,:dest]], fe(rr4))
fe_abel_yfe = hcat(gravity[:,[:year,:orig,:dest]], fe(rr5))


CSV.write(joinpath(@__DIR__,"../data/gravity_calib/beta.csv"), beta)
CSV.write(joinpath(@__DIR__,"../data/gravity_calib/beta_ratio.csv"), beta_ratio)

CSV.write(joinpath(@__DIR__,"../data/gravity_calib/fe_ar_yfe.csv"), fe_ar_yfe)
CSV.write(joinpath(@__DIR__,"../data/gravity_calib/fe_ar_odyfe.csv"), fe_ar_odyfe)
CSV.write(joinpath(@__DIR__,"../data/gravity_calib/fe_ratio_ar_yfe.csv"), fe_ratio_ar_yfe)
CSV.write(joinpath(@__DIR__,"../data/gravity_calib/fe_ratio_ar_odyfe.csv"), fe_ratio_ar_odyfe)
CSV.write(joinpath(@__DIR__,"../data/gravity_calib/fe_abel_yfe.csv"), fe_abel_yfe)

CSV.write("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/results_large/gravity.csv", gravity)
CSV.write("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/results_large/gravity_ar.csv", gravity_ar)
CSV.write("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/results_large/data_ar.csv", data_ar)

# Main specification: year fixed effects (reg_ar_yfe). All resulting files and graphs indexed _6.
# Robustness runs: origin/destination/year fixed effects (reg_ar_odyfe). All resulting files and graphs indexed _7.

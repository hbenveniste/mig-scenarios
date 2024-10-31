using CSV, DataFrames, DelimitedFiles
using Plots, VegaLite, FileIO, VegaDatasets, FilePaths
using Statistics, Query, Distributions, StatsPlots

# Approach: use Rao et al. (2018) and redo their Gini projections, assuming only education profiles and spending are affected by migration 
# (Hence no changes for TFP, non-resource trade, labor share of income, political orientation, public spending on health and social services)


################## Read data for original SSP Gini projections #####################################
gini = CSV.File(joinpath(@__DIR__, "../data/gini_rao/ssp_ginis.csv")) |> DataFrame


################## Prepare population data: original SSP and no-migration version ####################
# Working version provided by Samir KC in 2019
mig0_ssp1 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/mig0_SSP1.csv", DataFrame)
ssp1 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP1.csv", DataFrame)
mig0_ssp2 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/mig0_SSP2.csv", DataFrame)
ssp2 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP2.csv", DataFrame)
mig0_ssp3 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/mig0_SSP3.csv", DataFrame)
ssp3 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP3.csv", DataFrame)
mig0_ssp4 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/mig0_SSP4.csv", DataFrame)
ssp4 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP4.csv", DataFrame)
mig0_ssp5 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/mig0_SSP5.csv", DataFrame)
ssp5 = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP5.csv", DataFrame)

select!(mig0_ssp1, Not(:Column1))
select!(ssp1, Not(:Column1))
select!(mig0_ssp2, Not(:Column1))
select!(ssp2, Not(:Column1))
select!(mig0_ssp3, Not(:Column1))
select!(ssp3, Not(:Column1))
select!(mig0_ssp4, Not([:Column1, :pattern, :nSx, :pop1, :deaths, :pop2, :pop3, :asfr, Symbol("pop3.shift"), :births, Symbol("deaths.nb"), :pop4, :area]))
select!(ssp4, Not([:Column1, :pattern, :nSx, :pop1, :deaths, :pop2, :pop3, :asfr, Symbol("pop3.shift"), :births, Symbol("deaths.nb"), :pop4, :area]))
select!(mig0_ssp5, Not([:Column1, :pattern, :nSx, :pop1, :deaths, :pop2, :pop3, :asfr, Symbol("pop3.shift"), :births, Symbol("deaths.nb"), :pop4, :area]))
select!(ssp5, Not([:Column1, :pattern, :nSx, :pop1, :deaths, :pop2, :pop3, :asfr, Symbol("pop3.shift"), :births, Symbol("deaths.nb"), :pop4, :area]))

mig0_ssp4[!,:scen] = repeat(["mig0_SSP4"], size(mig0_ssp4,1))
ssp4[!,:scen] = repeat(["SSP4"], size(ssp4,1))
mig0_ssp5[!,:scen] = repeat(["mig0_SSP5"], size(mig0_ssp5,1))
ssp5[!,:scen] = repeat(["SSP5"], size(ssp5,1))

select!(mig0_ssp4, [2,1,3,4,5,6,7,8,9])
select!(ssp4, [2,1,3,4,5,6,7,8,9])
select!(mig0_ssp5, [2,1,3,4,5,6,7,8,9])
select!(ssp5, [2,1,3,4,5,6,7,8,9])

# Join ssp datasets
ssp = vcat(ssp1, ssp2, ssp3, ssp4, ssp5)
mig0 = vcat(mig0_ssp1, mig0_ssp2, mig0_ssp3, mig0_ssp4, mig0_ssp5)
mig0[!,:scen] = map(x -> SubString(x, 6:9), mig0[!,:scen])

# Replace missing values by zeros
for name in [:inmig, :outmig]
    for i in 1:length(ssp[!,name])
        if ssp[i,name] == "NA" 
            ssp[i,name] = "0"
        end
        if mig0[i,name] == "NA" 
            mig0[i,name] = "0"
        end
    end
    ssp[!,name] = map(x -> parse(Float64, x), ssp[!,name])
    mig0[!,name] = map(x -> parse(Float64, x), mig0[!,name])
end

# Sum projections for all sexes and ages: population per country and time period
ssp_cye = combine(groupby(ssp, [:region, :period, :scen, :edu]), d -> sum(d.pop))
rename!(ssp_cye, :x1 => :pop_mig)
mig0_cye = combine(groupby(mig0, [:region, :period, :scen, :edu]), d -> sum(d.pop))
rename!(mig0_cye, :x1 => :pop_nomig)

# Join population datasets
sspedu_6 = innerjoin(ssp_cye, mig0_cye, on = [:region, :period, :scen, :edu])


####################### Update SSP population scenarios #################################################
# Source:  K.C., S., Lutz, W. , Potančoková, M. , Abel, G. , Barakat, B., Eder, J., Goujon, A. , Jurasszovich, S., et al. (2020). 
# Global population and human capital projections for Shared Socioeconomic Pathways – 2015 to 2100, Revision-2018. 
# https://pure.iiasa.ac.at/id/eprint/17550/
ssp1_update = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP1_2018update.csv", DataFrame)
ssp2_update = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP2_2018update.csv", DataFrame)
ssp3_update = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP3_2018update.csv", DataFrame)
ssp4_update = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP4_2018update.csv", DataFrame)
ssp5_update = CSV.read("C:/Users/hmrb/Stanford_Benveniste Dropbox/Hélène Benveniste/YSSP-IIASA/Samir_data/SSP5_2018update.csv", DataFrame)

ssp1_update.scen = repeat(["SSP1"], size(ssp1_update,1))
ssp2_update.scen = repeat(["SSP2"], size(ssp2_update,1))
ssp3_update.scen = repeat(["SSP3"], size(ssp3_update,1))
ssp4_update.scen = repeat(["SSP4"], size(ssp4_update,1))
ssp5_update.scen = repeat(["SSP5"], size(ssp5_update,1))

ssp_update = vcat(ssp1_update, ssp2_update, ssp3_update, ssp4_update, ssp5_update)

# Keep only population numbers for all age groups together (ageno_0), both sexes together (sexno = 0) and disaggregated education levels (eduno != 0), and distinct countries (isono < 900)
filter!(
    row -> (row.sexno == 0 && row.eduno !=0 && row.isono < 900),
    ssp_update
)
select!(ssp_update, [:year,:isono,:eduno,:scen,:ageno_0])

# Full update: use population sizes and education distribution of updated scenarios
# Convert 10 education levels (Under 15, No Education, Incomplete Primary, Primary, Lower Secondary, Upper Secondary, Post Secondary, Short Post Secondary, Bachelor, Master and higher)
# to 6 education levels (no education, some primary, primary completed, lower secondary completed, upper secondary completed, post secondary completed)
ssp_update[!,:edu_6] = [(ssp_update[i,:eduno] == 1 || ssp_update[i,:eduno] == 2) ? "e1" : (ssp_update[i,:eduno] == 3 ? "e2" : (ssp_update[i,:eduno] == 4 ? "e3" : (ssp_update[i,:eduno] == 5 ? "e4" : (ssp_update[i,:eduno] == 6 ? "e5" : "e6")))) for i in eachindex(ssp_update[:,1])]

sspedu_6[!,:region] = map(x -> parse(Int, SubString(x, 3)), sspedu_6[!,:region])

sspedu_6 = leftjoin(
    sspedu_6,
    rename(
        combine(
            groupby(
                ssp_update, 
                [:year,:isono,:scen,:edu_6]
            ), 
            :ageno_0 => sum
        ),
        :year => :period, :isono => :region, :edu_6 => :edu, :ageno_0_sum => :pop_mig_update
    ),
    on = [:scen, :period, :region, :edu]
)

# We assume that for each country, education level, and SSP scenario, the ratio of pop_mig/pop_nomig remains the same for the update 
sspedu_6.pop_nomig_update = sspedu_6.pop_mig_update .* sspedu_6.pop_nomig ./ sspedu_6.pop_mig

# Partial update: use population sizes of updated scenarios and education distribution of earlier scenarios
sspedu_6 = innerjoin(
    sspedu_6,
    combine(
        groupby(
            sspedu_6,
            [:region,:period,:scen]
        ),
        :pop_mig => sum, :pop_nomig => sum 
    ),
    on = [:region,:period,:scen]
)
sspedu_6.popshare_mig = sspedu_6.pop_mig ./ sspedu_6.pop_mig_sum
sspedu_6.popshare_nomig = sspedu_6.pop_nomig ./ sspedu_6.pop_nomig_sum

# We assume that for each country and SSP scenario, the ratio of pop_mig_sum/pop_nomig_sum remains the same for the update 
# and that education distribution remains the same for the update 
sspedu_6 = leftjoin(
    sspedu_6,
    rename(
        combine(
            groupby(
                ssp_update, 
                [:year,:isono,:scen]
            ), 
            :ageno_0 => sum
        ),
        :year => :period, :isono => :region, :ageno_0_sum => :pop_mig_sum_update
    ),
    on = [:scen, :period, :region]
)

sspedu_6.pop_nomig_sum_update = sspedu_6.pop_mig_sum_update .* sspedu_6.pop_nomig_sum ./ sspedu_6.pop_mig_sum
sspedu_6.pop_mig_updatepart = sspedu_6.pop_mig_sum_update .* sspedu_6.popshare_mig
sspedu_6.pop_nomig_updatepart = sspedu_6.pop_nomig_sum_update .* sspedu_6.popshare_nomig


#################################### Compute changes in Gini due to migration-related changes in education composition of the population #####################
# Convert 6 education levels (no education, some primary, primary completed, lower secondary completed, upper secondary completed, post secondary completed)
# to 4 education levels (no education, primary, secondary, tertiary) following KC and Lutz (2017)
sspedu_6[!,:edu_4] = [(sspedu_6[i,:edu] == "e1") ? "noed" : ((sspedu_6[i,:edu] == "e2" || sspedu_6[i,:edu] == "e3") ? "prim" : ((sspedu_6[i,:edu] == "e4" || sspedu_6[i,:edu] == "e5") ? "sec" : "ter")) for i in eachindex(sspedu_6[:,1])]

# Compute changes in education shares of population related to migration
sspedu = combine(groupby(sspedu_6, [:region, :period, :scen, :edu_4]), d -> (pop_mig=sum(d.pop_mig_update), pop_nomig=sum(d.pop_nomig_update)))
sspedu_all = combine(groupby(sspedu, [:region, :period, :scen]), d -> (pop_mig_sum=sum(d.pop_mig), pop_nomig_sum=sum(d.pop_nomig)))
sspedu = innerjoin(sspedu, sspedu_all, on=[:region, :period, :scen])
sspedu[!,:edushare_mig] = sspedu[!,:pop_mig] ./ sspedu[!,:pop_mig_sum] 
sspedu[!,:edushare_nomig] = sspedu[!,:pop_nomig] ./ sspedu[!,:pop_nomig_sum]
sspedu[!,:edushare_diff] = sspedu[!,:edushare_nomig] .- sspedu[!,:edushare_mig]

# Compute resulting change to Gini according to Rao et al. (2018): beta2 coefficient in eq. (1) estimated in Table 2
beta2_prim = -0.27 ; beta2_sec = -0.26 ; beta2_ter = -0.61
sspedu[!,:beta2] = [(sspedu[i,:edu_4] == "prim") ? beta2_prim : ((sspedu[i,:edu_4] == "sec") ? beta2_sec : (sspedu[i,:edu_4] == "ter" ? beta2_ter : 0.0)) for i in eachindex(sspedu[:,1])]
dgini_edu = combine(groupby(sspedu, [:region, :period, :scen]), d -> sum(d.edushare_diff .* d.beta2))
rename!(dgini_edu, :x1 => :dgini_edu)


#################################### Compute changes in Gini due to migration-related changes in education spending ##############################
# Compute changes in education spending related to migration based on Rao et al. supplementary material
# Compute population index: educated population as a multiple of the educated population in 2010
# Calculate educated population in a given year with and without migration
sspeduspend = @from i in sspedu begin
    @where i.edu_4 != "noed"
    @select {i.region, i.period, i.scen, i.edu_4, i.pop_mig, i.pop_nomig}
    @collect DataFrame
end
educsum = combine(groupby(sspeduspend, [:region, :period, :scen]), d -> (pop_mig_sum=sum(d.pop_mig),pop_nomig_sum=sum(d.pop_nomig)))

# Calculate educated population in 2010 using data from the Wittgenstein Center 
educ2010_all = CSV.File(joinpath(@__DIR__, "../data/pop2010_wic/pop2010.csv");header=9) |> DataFrame
educ2010 = @from i in educ2010_all begin
    @where i.Education != "Total" && i.Education != "Under 15" && i.Education != "No Education"
    @select {i.Area, i.Year, i.Education, i.ISOCode, i.Population}
    @collect DataFrame
end
educsum_2010 = combine(groupby(educ2010, :ISOCode), d -> sum(d.Population))
rename!(educsum_2010, :x1 => :pop_2010, :ISOCode => :region)

# Calculate population index
educsum = innerjoin(educsum, educsum_2010, on=:region)
educsum[!,:popindex_mig] = educsum[!,:pop_mig_sum] ./ educsum[!,:pop_2010]
educsum[!,:popindex_nomig] = educsum[!,:pop_nomig_sum] ./ educsum[!,:pop_2010]

# Compute level index: costs of education for tertiary, secondary and primary, weighted by the population shares in each
# Consider average costs of primary, secondary and tertiary education per student for OECD countries (p.204 in https://www.oecd.org/edu/Education-at-a-Glance-2014.pdf)
cost_prim = 8296 ; cost_sec = 9280 ; cost_ter = 13958 ; cost_av = 9487
sspeduspend[!,:cost] = [(sspeduspend[i,:edu_4] == "prim") ? cost_prim : ((sspeduspend[i,:edu_4] == "sec") ? cost_sec : (sspeduspend[i,:edu_4] == "ter" ? cost_ter : 0.0)) for i in eachindex(sspeduspend[:,1])]
sspeduspend = innerjoin(sspeduspend, select(educsum, Not([:pop_2010, :popindex_mig, :popindex_nomig])), on=[:region, :period, :scen])
sspeduspend[!,:costweighted_mig] = sspeduspend[!,:pop_mig] ./ sspeduspend[!,:pop_mig_sum] .* sspeduspend[!,:cost] ./ cost_av
sspeduspend[!,:costweighted_nomig] = sspeduspend[!,:pop_nomig] ./ sspeduspend[!,:pop_nomig_sum] .* sspeduspend[!,:cost] ./ cost_av
educsum = innerjoin(educsum, combine(groupby(sspeduspend, [:region, :period, :scen]), d -> (levelindex_mig = sum(d.costweighted_mig), levelindex_nomig = sum(d.costweighted_nomig))), on=[:region, :period, :scen])

# Compute resource index: difference in priority accorded to education in public spending; only SSP1 and SSP5 up to 2050, and only to those countries whose average education spending is below the global average in 2010
# Use data on per capita education spending from SPEED (Statistics of Public Expenditure for Economic Development) from the International Food Policy Research Institute
speed = load(joinpath(@__DIR__, "../data/eduspending_ifpri/speed.xls"), "gdpeducation_ppp!A1:AK148") |> DataFrame
select!(speed, [:region, :ISO, Symbol("2010")])
rename!(speed, Symbol("2010") => :educspend2010)

# Complete data: add countries in Gini missing in SPEED + when missing data on education spending, use regional average
# Regions used in SPEED: East Asia and Pacific (EAP), Europe Central Asia (ECA), Latin America and the Caribbean (LAC), MENA, SOUTH ASIA, SSA, EURO ZONE, HIGH INCOME 
c_miss = setdiff(unique(gini[:,:iso]), unique(speed[:,:ISO]))
s_miss = DataFrame(
    region = ["LAC","ECA","ECA","EAP","SSA","LAC","SSA","SSA","SSA","SSA","LAC","EAP","LAC","LAC","MENA","EAP","EAP","MENA","LAC","EAP","ECA","ECA","SSA","EAP","LAC","LAC","EAP","MENA","EAP","SSA","SSA","LAC","SSA","ECA","ECA","EAP","EAP","ECA","EAP"],
    ISO = c_miss,
    educspend2010 = missings(length(c_miss))
)
speed = vcat(speed, s_miss)

speed_reg = load(joinpath(@__DIR__, "../data/eduspending_ifpri/speed.xls"), "gdpeducation_ppp!B151:AK160") |> DataFrame
select!(speed_reg, [:region, Symbol("2010")])
rename!(speed_reg, Symbol("2010") => :educspend2010)
for i in eachindex(speed[:,1])
    if ismissing(speed[i,:educspend2010])
        ind = findfirst(speed_reg[:,:region].==speed[i,:region])
        speed[i,:educspend2010] = speed_reg[ind,:educspend2010]
    end
end

iso3c_isonum = CSV.File(joinpath(@__DIR__,"../data/iso3c_isonum.csv")) |> DataFrame
rename!(iso3c_isonum, :iso3c => :ISO)
speed = innerjoin(speed, iso3c_isonum, on = :ISO)
educsum = leftjoin(educsum, rename(select(speed, Not([:ISO,:region])), :isonum => :region), on =:region)
av_world = speed_reg[findfirst(speed_reg[:,:region].=="World"),:educspend2010]
educsum[!,:resindex] = [((educsum[i,:scen] == "SSP1" || educsum[i, :scen] == "SSP5") && educsum[i,:period] <= 2050 && !ismissing(educsum[i,:educspend2010]) && educsum[i,:educspend2010] < av_world) ? 0.1 * (educsum[i, :period]- 2010) : 0.0 for i in eachindex(educsum[:,1])]

# Compute education spending as percentage of total public spending for each country, year and scenario, with and without migration
educsum[!,:educspend_mig] = (educsum[!,:educspend2010] .+ educsum[!,:resindex]) ./ 100 .* educsum[!, :popindex_mig] .* educsum[!, :levelindex_mig]
educsum[!,:educspend_nomig] = (educsum[!,:educspend2010] .+ educsum[!,:resindex]) ./ 100 .* educsum[!, :popindex_nomig] .* educsum[!, :levelindex_nomig]

# Compute changes in education spending related to migration
educsum[!,:educspend_diff] = educsum[!,:educspend_nomig] .- educsum[!,:educspend_mig]

# Compute resulting change to Gini according to Rao et al. (2018): beta5 coefficient in eq. (1) estimated in Table 2
beta5_edu = 0.67
dgini_eduspend = select(educsum, [:region, :period, :scen, :educspend_diff])
dgini_eduspend[!,:dgini_eduspend] = dgini_eduspend[!,:educspend_diff] .* beta5_edu


#################################### Perform new Gini projections including migration-related changes in education profiles and spending #########
# Join gini datasets
rename!(iso3c_isonum, :ISO => :iso)
gini = innerjoin(gini, iso3c_isonum, on = :iso)
rename!(gini, :year => :period, :scenario => :scen, :isonum => :region)
gini = leftjoin(gini, dgini_edu, on =[:region, :period, :scen])
gini = leftjoin(gini, dgini_eduspend, on =[:region, :period, :scen])
select!(gini, Not(:educspend_diff))

# Compute Gini projections without migration
gini[!,:gini_nomig] = gini[!,:gini] .- gini[!,:dgini_edu] .- gini[!, :dgini_eduspend]

# For 2010 (before migration starts), assign same values of Gini for mig and nomig
for i in eachindex(gini[:,1])
    if gini[i,:period] == 2010
        gini[i,:gini_nomig] = gini[i, :gini]
    end
end

sort!(gini, [:scen,:period,:iso])
rename!(gini, :iso=>:country)


################################################# Plot geographical maps #####################################
ssps = unique(gini[!,:scen]) ; countries = unique(gini[!,:country]) ; periods = unique(gini[!,:period])[1:17]
world110m = dataset("world-110m")

# Look at variations in Gini without migration
gini[!,:relgini] = gini[!,:gini_nomig] ./ gini[!,:gini] .- 1

for s in ssps
    @vlplot(width=800, height=600) + @vlplot(mark={:geoshape, stroke = :lightgray}, 
        data={values=world110m, format={type=:topojson, feature=:countries}}, 
        transform = [{lookup=:id, from={data=filter(row -> row[:scen] == s && row[:period] == 2100, gini), key=:region, fields=["relgini"]}}],
        projection={type=:naturalEarth1}, 
        color = {"relgini:q", scale={domain=[-0.25,0.25], scheme=:redblue}, legend={title="Relative change", symbolSize=40, labelFontSize=16}}
    ) |> save(joinpath(@__DIR__, "../results/world_maps/", string("ginidiff_", s, "_edu_update.png")))
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

# For each continent, show median of change in Gini for each SSP with and without migration, plus interquartile range of relative changes
gini = innerjoin(gini, regions, on =:country)
gini[!,:reldif] = gini[!,:gini] ./ gini[!,:gini_nomig] .- 1

for r in unique(regions[!,:worldregion])
    gini |> @filter(_.worldregion == r) |> @vlplot(
        width=300, height=250, x = {"period:o", axis={labelFontSize=16}, title=nothing}
    ) + @vlplot(
        mark={:errorband, extent=:iqr}, y={"gini:q", title="Gini with migration", scale={domain=[0.2,0.7]}, axis={labelFontSize=16,titleFontSize=16}},
        title = string(r), 
        color = {"scen:n", scale={scheme=:category10}, legend=nothing}, 
    ) +
    @vlplot(
        :line, y={"median(gini)", title="Gini with migration"},
        color = {"scen:n", scale={scheme=:category10}, legend=nothing}, 
    ) |> save(joinpath(@__DIR__, "../results/gini_reg/", string("ginici_mig_", r,"_edu_update.png")))
end
for r in unique(regions[!,:worldregion])
    gini |> @filter(_.worldregion == r) |> @vlplot(
        width=300, height=250, x = {"period:o", axis={labelFontSize=16}, title=nothing}
    ) + @vlplot(
        mark={:errorband, extent=:iqr}, y={"gini_nomig:q", title="Gini without migration", scale={domain=[0.2,0.7]}, axis={labelFontSize=16,titleFontSize=16}},
        title = string(r), 
        color = {"scen:n", scale={scheme=:category10}, legend=nothing}, 
    ) +
    @vlplot(
        :line, y={"median(gini_nomig)", title="Gini without migration"},
        color = {"scen:n", scale={scheme=:category10}, legend=nothing}, 
    ) |> save(joinpath(@__DIR__, "../results/gini_reg/", string("ginici_nomig_", r,"_edu_update.png")))
end
for r in unique(regions[!,:worldregion])
    gini |> @filter(_.worldregion == r) |> @vlplot(
        width=300, height=250, x = {"period:o", axis={labelFontSize=16}, title=nothing}
    ) + @vlplot(
        mark={:errorband, extent=:iqr}, y={"reldif:q", title="Relative change with migration", axis={labelFontSize=16,titleFontSize=16}},
        title = string(r), 
        color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    ) +
    @vlplot(
        :line, y={"median(reldif)", title="Relative change with migration"},
        color = {"scen:n", scale={scheme=:category10}, legend={titleFontSize=16, symbolSize=40, labelFontSize=16}}, 
    ) |> save(joinpath(@__DIR__, "../results/gini_reg/", string("ginicirel_", r,"_edu_update.png")))
end


####################################### Write output files with results ############################################################################
CSV.write(joinpath(@__DIR__, "../results/gini_edu_update.csv"), gini)
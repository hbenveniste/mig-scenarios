using CSV, DataFrames, ExcelFiles, Query, DelimitedFiles, Statistics, XLSX


##################### Calculating the share of income that migrants send to their home region as remittances ##################
# Reading remittances flows at country * country level. Data for 2017 from World Bank
remflow_matrix = XLSX.openxlsx(joinpath(@__DIR__, "../data/rem_wb/WB_Bilateral_Remittance_Estimates_2017.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["Bilateral_Remittances_2017"];first_row=2)...)
end
countriesr = remflow_matrix[1:214,1]
remflow = stack(remflow_matrix, 2:215)
select!(remflow, Not(:WORLD))
rename!(remflow, Symbol("receiving (across) / sending (down) ") => :sending, :variable => :receiving, :value => :flow)
sort!(remflow, :sending)
indworld = findall(remflow[!,:sending] .== "WORLD")
delete!(remflow, indworld)
indmissing = findall([typeof(remflow[i,:flow]) != Float64 for i in 1:size(remflow, 1)])
for i in indmissing ; remflow[i,:flow] = 0.0 end
remflow[!,:flow] = map(x -> float(x), remflow[!,:flow])
remflow[!,:receiving] = map(x -> string(x), remflow[!,:receiving])

# Reading migrant stocks at country * country level. Data for 2017 from World Bank
migstock_matrix = XLSX.openxlsx(joinpath(@__DIR__, "../data/rem_wb/WB_Bilateral_Estimates_Migrant_Stocks_2017.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["Bilateral_Migration_2017"];first_row=2)...)
end
countriesm = migstock_matrix[1:214,1]
migstock = stack(migstock_matrix, 2:215)
select!(migstock, Not([Symbol("Other North"), Symbol("Other South"), :World]))
rename!(migstock, Symbol("destination (across) / origin (down) ") => :origin, :variable => :destination, :value => :stock)
sort!(migstock, :origin)
indregion = vcat(findall(migstock[!,:origin] .== "Other North"), findall(migstock[!,:origin] .== "Other South"), findall(migstock[!,:origin] .== "World"))
delete!(migstock, indregion)
indmissing = findall([ismissing(migstock[i,:stock]) for i in 1:size(migstock, 1)])
for i in indmissing ; migstock[i,:stock] = 0.0 end
migstock[!,:stock] = map(x -> float(x), migstock[!,:stock])
migstock[!,:destination] = map(x -> string(x), migstock[!,:destination])

# Reading GDP per capita at country level. Data for 2017 from World Bank(WDI), in current USD
ypc2017 = XLSX.openxlsx(joinpath(@__DIR__, "../data/rem_wb/GDPpercap2017.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["Data"])...)
end
select!(ypc2017, Not([Symbol("Series Code"), Symbol("Series Name")]))
rename!(ypc2017, Symbol("Country Name") => :country, Symbol("Country Code") => :country_code, Symbol("2017 [YR2017]") => :ypc)
for i in eachindex(ypc2017[:,1]) ; if ypc2017[i,:ypc] == ".." ; ypc2017[i,:ypc] = missing end end      # replacing missing values by zeros

# Joining data
rename!(remflow, :sending => :destination) ; rename!(remflow, :receiving => :origin)                # remittances sending country = destination country
rho = outerjoin(remflow, migstock, on = [:origin, :destination])
sort!(rho, :origin)
rename!(rho, :flow => :remflow, :stock => :migstock)

rename!(ypc2017, :country => :destination, :country_code => :code_dest, :ypc => :ypc_dest)
indnkorea = findfirst(x -> x == "Korea, Dem. People’s Rep.", ypc2017[!,:destination])
ypc2017[!,:destination][indnkorea] = "Korea, Dem. Rep."
rho = innerjoin(rho, ypc2017, on = :destination)         
rename!(ypc2017, :destination => :origin, :code_dest => :code_or, :ypc_dest => :ypc_or)
rho = innerjoin(rho, ypc2017, on = :origin)       
rho[!,:ypc] = [max(mean([rho[i,:ypc_or],rho[i,:ypc_dest]]), rho[i,:ypc_or]) for i in eachindex(rho[:,1])]

# Calculating rho using rho * ypc * migstock = remflow
rho[!,:rho] = rho[!,:remflow] .* 1000000 ./ (rho[!,:migstock] .* rho[!,:ypc])       # Remittances are in million USD 2018
for i in eachindex(rho[:,1])
    if ismissing(rho[i,:ypc]) || rho[i,:migstock] == 0.0 || rho[i,:ypc] == 0.0
        rho[i,:rho] = 0.0
    end
end

# Sorting the data
rhofinal = rho[:,[:code_or,:code_dest,:rho]]
sort!(rhofinal, [:code_or,:code_dest])
rename!(rhofinal, :code_or => :origin, :code_dest => :destination)


CSV.write(joinpath(@__DIR__, "../data/rem_wb/rho.csv"), rhofinal)
using CSV, DataFrames, ExcelFiles, Query, Statistics, XLSX


####################### Calculating the cost of sending remittances, as share of money sent ###########################
# Reading remittances price data by World Bank: Remittance Prices Worldwide (2018)
# Aggregate data at the country * country level for 2017. 
rpw = XLSX.openxlsx(joinpath(@__DIR__, "../data/rem_wb/WB_rpw.xlsx")) do xf
    DataFrame(XLSX.gettable(xf["WB_rpw"])...)
end
select!(rpw, Not(vcat(names(rpw)[1], names(rpw)[4:8], names(rpw)[10:22], names(rpw)[27:29], names(rpw)[34:end])))
# Source is the sending country, i.e. the migrant's destination
rename!(rpw, Symbol("cc1 total cost %") => :cc1, Symbol("cc2 total cost %") => :cc2, :source_code => :destination, :destination_code => :origin)
rpw[!,:period] = map(x -> parse(Int, SubString(x,1:4)), rpw[!,:period])
indnot2017 = findall(rpw[!,:period] .!= 2017)
delete!(rpw, indnot2017)

# Drop rows for which costs are negative
rpw[!,:meancost] = [((typeof(rpw[i,:cc1]) == Float64 || typeof(rpw[i,:cc1]) == Int64) && rpw[i,:cc1]>=0) ? ((typeof(rpw[i,:cc2]) == Float64 && rpw[i,:cc2]>=0) ? mean([rpw[i,:cc1], rpw[i,:cc2]]) : max(0,rpw[i,:cc1])) : max(0,rpw[i,:cc2]) for i in eachindex(rpw[:,1])]
# Average over surveys and firms per corridor
phi = combine(groupby(rpw, [:origin, :destination]), df -> mean(df[!,:meancost]) ./ 100)      
rename!(phi, :x1 => :phi)


CSV.write("../data/rem_wb/phi.csv", phi)
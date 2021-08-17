Background details:

Please cite the International Migration Review article if you use any observations from the data. 

The data is in a `tidy` format, with 10,321,454 rows (observations) and 10 columns (variables). 

NOTE1: Due to the large number of rows the file is unlikey to be fully displayed in MS Excel.
NOTE2: There are many caveats to the estimates which are outlined in the International Migration Review article.

Column details:

stock: stock data source used for the estimated flow (un12, un13, un15 or wb11)
demo: demographic data source used for the estimated flow (wpp2010, wpp2012, wpp2015)
sex: gender of the estimated flow
year0: first year of the period of the estimated flow (ranging from 1960 to 2010, and varying in lengths depending on the stock and demo data, as described in the paper)
interval: the length of the period of the estimated flow 
orig: ISO 3166-1 alpha-3 letter code for the origin country of the estimated flow
dest: ISO 3166-1 alpha-3 letter code for the destination country of the estimated flow
orig_code: ISO 3166-1 numeric code for the origin country of the estimated flow
dest_code: ISO 3166-1 numeric code for the destination country of the estimated flow
flow: estimated flow

More details on the country codes:

https://www.iso.org/iso-3166-country-codes.html
https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3
https://en.wikipedia.org/wiki/ISO_3166-1_numeric

More details on tidy data:

http://vita.had.co.nz/papers/tidy-data.pdf
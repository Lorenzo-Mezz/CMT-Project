co2_future_mod = co2_future(11:310, :)
csvco2andtemp = table(future_dates', Y_future, co2_future_mod.CO2, bottom_pred, 'VariableNames', {'Time','Surface Temp','Atm. CO2', 'Bottom Temp'})
writetable(csvco2andtemp, 'co2_temp_for_c.csv')
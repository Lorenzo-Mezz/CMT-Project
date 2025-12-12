clear

run("surface_temp_zug.m")
pause(1)
run("bottomtempzug.m") 
pause(1)
run("chl_a_zug.m")
pause(1)
run("CO2curves.m")
pause(1)
co2_future_mod = co2_future(11:310, :);
csvco2andtemp = table(future_dates', Y_future, bottom_pred,co2_future_mod.CO2, biomass_future, 'VariableNames', {'Time','Surface Temp', 'Bottom Temp','Atm. CO2', 'Surface Biomass (micrograms/m^3)' })
writetable(csvco2andtemp, '../data/co2_temp_for_c.csv')
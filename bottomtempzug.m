sf_temp_zug2 = readtable("CTD_monitoring_Lake_Zug_by_Canton_Zugbottom")
plot(sf_temp_zug2.Time, sf_temp_zug2.BottomTemperature_degC_)

sf_apres2006 = sf_temp_zug2(sf_temp_zug2.Time > datetime(2006,03,01), :)
time = convertTo(sf_apres2006.Time, 'juliandate');
t_0 = min(time)
time = time - t_0
bot_model = fitlm(time,sf_apres2006.BottomTemperature_degC_)

future_dates = datetime(2025,11,01):calmonths(1):datetime(2050,10,31);
t_future = convertTo(future_dates, 'juliandate');
t_future = t_future - t_0;
t_future = t_future';

bottom_pred = predict(bot_model, t_future);

%plot(sf_apres2006.Time, sf_apres2020.BottomTemperature_degC_)
plot(sf_temp_zug2.Time, sf_temp_zug2.BottomTemperature_degC_)
hold on
plot(future_dates, bottom_pred)
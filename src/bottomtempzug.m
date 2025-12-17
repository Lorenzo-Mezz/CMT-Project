% read CSV file and plot to get an idea of trend
sf_temp_zug2 = readtable("../data/CTD_monitoring_Lake_Zug_by_Canton_Zugbottom");
plot(sf_temp_zug2.Time, sf_temp_zug2.BottomTemperature_degC_)

% noticed trend in graph, decided to cut off before march 2026
sf_apres2006 = sf_temp_zug2(sf_temp_zug2.Time > datetime(2006,03,01), :);
time = convertTo(sf_apres2006.Time, 'juliandate');
t_0 = min(time);
time = time - t_0;

% fit linear regression model
bot_model = fitlm(time,sf_apres2006.BottomTemperature_degC_);

% define future dates
future_dates = datetime(2025,11,01):calmonths(1):datetime(2050,10,31);

% using t_future from previous code
t_future = convertTo(future_dates, 'juliandate');
t_future = t_future - t_0;
t_future = t_future';

% predicitons of bottom temperature using linear regression model
bottom_pred = predict(bot_model, t_future);

% plot of historical data + predicitons
plot(sf_temp_zug2.Time, sf_temp_zug2.BottomTemperature_degC_)
hold on
plot(future_dates, bottom_pred)

xlabel('Time');
ylabel('Temperature');
title('Bottom Lake Temperature (195 m depth)');

saveas(gcf, "../data/bottom_temp_plot.png");

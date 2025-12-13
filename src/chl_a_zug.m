chlorophyll_zug = readtable("../data/zug_chlorophyll_satellite_summary.csv");
chl_data = split(chlorophyll_zug.Var1, ',');

dt = datetime(chl_data(:,1), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSX', 'TimeZone','local');
chl_a = str2double(chl_data(:,3))
chl_data = table(dt, chl_a)
dt = datetime(dt, 'TimeZone', 'local');

keep_index = chl_data.chl_a <= 18;
chl_data_filtered = chl_data(keep_index, :);

dt_filtered = chl_data_filtered.dt;
chl_a_filtered = chl_data_filtered.chl_a;

[dt_sorted, idx] = sort(dt_filtered);
dt_sorted.TimeZone = 'UTC';
chl_a_sorted = chl_a_filtered(idx);

num_removed = size(chl_data, 1) - size(chl_data_filtered, 1);
disp(['Total data points removed (chl_a > 20): ', num2str(num_removed)]);

moyenne_chla = mean(chl_a);
disp(['Average Chlorophyll Concentration: ', num2str(moyenne_chla)]);


t = convertTo(dt_sorted, 'juliandate');
t0 = min(t);
t = t-t0;

cpy = 1;
f = cpy ./ 365.25;
omega = 2.*pi.*f;

x = [
    t,...
    sin(omega*t), ...
    cos(omega*t), ...
    ];

returned = fitlm(x,chl_a_sorted);



figure
plot(dt_sorted, chl_a_sorted)
line(dt_sorted, returned.Fitted, 'Color', 'r')

coef = returned.Coefficients.Estimate


future_dates = datetime(2025,11,01):calmonths(1):datetime(2050,10,31);
future_dates.TimeZone = 'UTC';

t_future = convertTo(future_dates, 'juliandate');
t_future = t_future - min(t_future);
t_future = t_future';

x_future = [
    t_future,...
    sin(omega*t_future), ...
    cos(omega*t_future), ...
    ];, 'Color','b'

coef_1 = coef(1);

coef_1_matrix = repmat(coef_1, 300, 1);

chl_a_future = coef_1_matrix + x_future * coef(2:end)  ;

figure
plot(dt_sorted, chl_a_sorted, 'Color','b')
hold on
plot(future_dates, chl_a_future, 'Color','b')

figure
biomass_future = chl_a_future/0.015;
plot(future_dates', biomass_future)

xlabel('Time');
ylabel('Biomass (micrograms/m3)');
title('Biomass Prediciton at Surface');

saveas(gcf, "../data/chl_a_pred.png");

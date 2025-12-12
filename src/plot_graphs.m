growth_rate_pred = readtable("../data/biomass_growth_results.csv");

dtime = growth_rate_pred.Time;
h_0 = growth_rate_pred.Depth_0;
h_1 = growth_rate_pred.Depth_1;
h_2 = growth_rate_pred.Depth_2;
h_5 = growth_rate_pred.Depth_5;
h_10 = growth_rate_pred.Depth_10;
h_30 = growth_rate_pred.Depth_30;

figure;
plot(dtime, h_0, 'DisplayName', '0 m'); hold on;
plot(dtime, h_1, 'DisplayName', '1 m');
plot(dtime, h_2, 'DisplayName', '2 m');
plot(dtime, h_5, 'DisplayName', '5 m');
plot(dtime, h_10, 'DisplayName', '10 m');
plot(dtime, h_30, 'DisplayName', '30 m');

xlabel('Time');
ylabel('Growth');
title('Algal Growth at varying Depths');

legend;

grid on;
saveas(gcf, "../results/algal_growth_plot.png");

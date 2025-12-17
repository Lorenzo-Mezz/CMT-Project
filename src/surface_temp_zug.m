% read CSV file
sf_temp_zug = readtable("../data/CTD_monitoring_Lake_Zug_by_Canton_Zug");

% plot to get a visualisation
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_)

t = convertTo(sf_temp_zug.Time, 'juliandate');
t = t-min(t);

% seasonal frequencies : one cycle per year for the main seasonal change and a second cycle per year to better match the shape of the seasonal ups and downs
cpy = 1:2;
f = cpy ./ 365.25;
omega = 2.*pi.*f;

% building the regression design matrix:
% - linear trend
% - quadratic trend
% - annual sinusoid (1 cycle/year)
% - semi-annual sinusoid (2 cycles/year)
X = [
    t,...
    t.^2,...
    sin(omega(1)*t), ...
    cos(omega(1)*t), ...
    sin(omega(2)*t), ...
    cos(omega(2)*t),...
    ];

% observerved surface temp and fitting harmonic regression model with polynomial trend
Y = sf_temp_zug.SurfaceTemperature_degC_;
returned = fitlm(X,Y);

% plotting observed data and fitted model
figure
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_);
line(sf_temp_zug.Time, returned.Fitted, 'Color', 'r');

% extract regression coefficients
coef = returned.Coefficients.Estimate;

% define future time intervals/dates
future_dates = datetime(2025,11,01):calmonths(1):datetime(2050,10,31);
t_future = convertTo(future_dates, 'juliandate');
t_future = t_future - min(t_future);
t_future = t_future';

% building the design matrix for future predictions
X_future = [
    t_future,...
    t_future.^2,...
    sin(omega(1)*t_future), ...
    cos(omega(1)*t_future), ...
    sin(omega(2)*t_future), ...
    cos(omega(2)*t_future)
    ]; 


% extract intercept term + replicate intercept to match number of future points
coef_1 = coef(1);
coef_1_matrix = repmat(coef_1, 300, 1);

% future surface temperature predictions
Y_future = coef_1_matrix + X_future * coef(2:end);

% plot of historical and predicted surface temperature
figure
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_, 'Color','b');
hold on
plot(future_dates, Y_future, 'Color','b');
saveas(gcf, "../data/surface_temp_plot.png");



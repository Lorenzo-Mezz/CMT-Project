sf_temp_zug = readtable("CTD_monitoring_Lake_Zug_by_Canton_Zug");
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_)

t = convertTo(sf_temp_zug.Time, 'juliandate');
t = t-min(t);

cpy = 1:2;
f = cpy ./ 365.25;
omega = 2.*pi.*f;

X = [
    t,...
    t.^2,...
    sin(omega(1)*t), ...
    cos(omega(1)*t), ...
    sin(omega(2)*t), ...
    cos(omega(2)*t),...
    ];

Y = sf_temp_zug.SurfaceTemperature_degC_;
returned = fitlm(X,Y);

figure
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_);
line(sf_temp_zug.Time, returned.Fitted, 'Color', 'r');

coef = returned.Coefficients.Estimate;

future_dates = datetime(2025,11,01):calmonths(1):datetime(2050,10,31);
t_future = convertTo(future_dates, 'juliandate');
t_future = t_future - min(t_future);
t_future = t_future';

X_future = [
    t_future,...
    t_future.^2,...
    sin(omega(1)*t_future), ...
    cos(omega(1)*t_future), ...
    sin(omega(2)*t_future), ...
    cos(omega(2)*t_future)
    ]; 

coef_1 = coef(1);

coef_1_matrix = repmat(coef_1, 300, 1);

Y_future = coef_1_matrix + X_future * coef(2:end);

figure
plot(sf_temp_zug.Time, sf_temp_zug.SurfaceTemperature_degC_, 'Color','b');
hold on
plot(future_dates, Y_future, 'Color','b');




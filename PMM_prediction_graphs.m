data = load('EOP_1900tonow_MJD_xy_dxsigmadx_dysigmady.txt'); % EOPc04 series data pulled from the IERS (1900-now)
x = data(:,2); % daily final vaules of xp measured in mas
y = data(:,4); % daily final vaules of yp measured in mas

xsq = x.^2;
ysq = y.^2;
sigmaxsq = data(:,3).^2; % daily final values of sigmax measured in mas
sigmaysq = data(:,5).^2; % daily final values of sigmay measured in mas

uncert_r = sqrt((xsq.*sigmaxsq)./(xsq+ysq) + (ysq.*sigmaysq)./(xsq+ysq));

r = sqrt(x.^2 + y.^2);

t = zeros(2454,1); % create an array of zeros the length of the data given
for i = 2:length(t) % starting from t2 to then end
    t(i) = t(i-1) + .05*3.154e7; % data point taken every .05 year
    % need time in sec for PMM so 3.1547e7 sec in a year. Add time past in
    % sec from when the last one happened
end

A = [cos(2*pi*f(1)*t) sin(2*pi*f(1)*t) cos(2*pi*f(2)*t) sin(2*pi*f(2)*t)...
    cos(2*pi*f(3)*t) sin(2*pi*f(3)*t) cos(2*pi*f(4)*t) sin(2*pi*f(4)*t)...
    cos(2*pi*f(5)*t) sin(2*pi*f(5)*t)]; % Comes from Polar motion model

covA = [cos(2*pi*f(1)*t/uncertP1) sin(2*pi*f(1)*t/uncertP1) cos(2*pi*f(2)*t/uncertP2) sin(2*pi*f(2)*t/uncertP2)...
    cos(2*pi*f(3)*t/uncertP3) sin(2*pi*f(3)*t/uncertP3) cos(2*pi*f(4)*t/uncertP4) sin(2*pi*f(4)*t/uncertP4)...
    cos(2*pi*f(5)*t/uncertP5) sin(2*pi*f(5)*t)/uncertP5]; % Comes from Polar motion model

uncertPMM = sqrt(diag(cov(inv(covA.'*covA)))); % uncertainty in the PMM

% setting up to solve for constants
constants = (A.'*A)\(A.'*r);

preddata = load('2022_8_to_2022_10_MJD_x_dx_y_dy.txt');
xpred = preddata(:,2)*10^-3;
ypred = preddata(:,4)*10^-3;

rpred = sqrt(xpred.^2 + ypred.^2);

xsq = xpred.^2;
ysq = ypred.^2;
sigmaxsq = (preddata(:,3)*10^-3).^2; % daily final values of sigmax measured in mas
sigmaysq = (preddata(:,5)*10^-3).^2; % daily final values of sigmay measured in mas

uncert_r = sqrt((xsq.*sigmaxsq)./(xsq+ysq) + (ysq.*sigmaysq)./(xsq+ysq));

Npred = length(preddata);

tpred = zeros(2454+Npred,1); % create an array of zeros the length of the data given
for i = 2:length(tpred) % starting from t2 to then end
    tpred(i) = tpred(i-1) + .003*3.154e7; % data point taken every .05 year
    % need time in sec for PMM so 3.1547e7 sec in a year. Add time past in
    % sec from when the last one happened
end

tpred = tpred(2454+1:2490);

Apred = [cos(2*pi*f(1)*tpred) sin(2*pi*f(1)*tpred) cos(2*pi*f(2)*tpred) sin(2*pi*f(2)*tpred)...
    cos(2*pi*f(3)*tpred) sin(2*pi*f(3)*tpred) cos(2*pi*f(4)*tpred) sin(2*pi*f(4)*tpred)...
    cos(2*pi*f(5)*tpred) sin(2*pi*f(5)*tpred)]; % Comes from Polar motion model

covApred = [cos(2*pi*f(1)*tpred/uncertP1) sin(2*pi*f(1)*tpred/uncertP1) cos(2*pi*f(2)*tpred/uncertP2) sin(2*pi*f(2)*tpred/uncertP2)...
    cos(2*pi*f(3)*tpred/uncertP3) sin(2*pi*f(3)*tpred/uncertP3) cos(2*pi*f(4)*tpred/uncertP4) sin(2*pi*f(4)*tpred/uncertP4)...
    cos(2*pi*f(5)*tpred/uncertP5) sin(2*pi*f(5)*tpred/uncertP5)]; % Comes from Polar motpredion model

uncertPMM = sqrt(diag(cov(inv(covApred.'*covApred))));

alpha = -10^-1.9
pmm_check = Apred.'*Apred
pmm_pred = Apred.'*Apred*constants;
chir = rpred%./uncert_r
% note these are the equations without the uncertainties in PMM and r
actual = covApred.'*(chir)*alpha;
predicted = ((pmm_pred)./uncertPMM)%*alpha;

chisq = sum(actual - predicted)^2

% plots subplots of points of what PMM predicts and what r actually is
rfuture = (Apred.'*rpred)*alpha;
figure(1)
subplot(2,1,1)
plot(pmm_pred,'.r');
xlim([0,10])
ylim([-.2,.2])
title('PMM Prediction')
subplot(2,1,2)
plot(rfuture,'.b');
xlim([0,10])
ylim([-.2,.20])
title('A^{T}r_{Bulletin B}')

% combined plots of PMM and r
figure(2)
ePMM = errorbar(pmm_pred, uncertPMM);
hold on
plot(rfuture,'.r')
%ylim([-.25,.3])
title('Combined PMM Prediction and A^{T}r_{Bulletin B}')
legend({'PMM prediction';'r'},'Location','northwest')
hold off






%% FFT of EOP Data
data = load('EOP_1900tonow_MJD_xy_dxsigmadx_dysigmady.txt'); % EOPc04 series data pulled from the IERS (1900-now)
x = data(:,2); % daily final vaules of xp measured in mas
y = data(:,4); % daily final vaules of yp measured in mas
xsq = x.^2;
ysq = y.^2;
sigmaxsq = data(:,3).^2; % daily final values of sigmax measured in mas
sigmaysq = data(:,5).^2; % daily final values of sigmay measured in mas

uncert_r = sqrt((xsq.*sigmaxsq)./(xsq+ysq) + (ysq.*sigmaysq)./(xsq+ysq));

N = length(data); % Number of data points taken
fs = (N/122)/3.154e7; % sampling frequency -here one sample is taken every
% .05 years - want in Hz so 1/sec, 3.154e7 sec in a year

freq = (-N/2:N/2-1)/N*fs; % creates an 1byN array of frequencies where the
% spectrum is centered in the middle

r = sqrt(x.^2 + y.^2);

R = fft(r)/N; % scaling factor for the y axis


plot(freq/10^-8,abs(fftshift(R)));
title('Frequencies from Fast Fourier Transform of r')
xlabel('Frequencies x10^{-8}')
ylabel('|FFT of r|')
xlim([.06,3.5])
ylim([-.001,.05])
% fftshift rearranges the fourier transform by
% shifting the zero frequency component to the center of the array

%% Pick Out Frequencies for Least Squares

ypts = abs(fftshift(R));
peaksyindx = [];
n = length(ypts);

if ypts(1)> ypts(2) % sees if the first value is a peak
    peaksyindx = vertcat(peaksyindx,1);
end

for i = 2:n-1
    if ypts(i-1)<ypts(i) && ypts(i+1)<ypts(i) % finds peaks by checking if it's neighbors are less than
        peaksyindx = vertcat(peaksyindx,i);
    end
end

if ypts(n)> ypts(n-1) % sees if the last value is a peak
   peaksyindx = vertcat(peaksyindx,n);
end

ypeaks = [ypts(peaksyindx).';peaksyindx.']; % Creates a 2xn matrix that holds
% the y values in the first row and the y indices in the second

[m,n] = size(ypeaks); % picking out the length to iterate through

xindx = [];
for i = 1:n
    if ypeaks(1,i) > .012 % Sorts through only significant peaks
        xindx = vertcat(xindx, ypeaks(2,i)); % saves indices to find frequencies
    end
end

f = freq(xindx); % picks out all the frequencies
f = f(7:11); % because of symmetry only want/need positive side and convert to w for PMM
fday = (f)*(3.154e7)*365; % for chandler wobble term analysis

%% Least Squares Calculation


t = zeros(2454,1); % create an array of zeros the length of the data given
for i = 2:length(t) % starting from t2 to then end
    t(i) = t(i-1) + .05*3.154e7; % data point taken every .05 year
    % need time in sec for PMM so 3.1547e7 sec in a year. Add time past in
    % sec from when the last one happened
end

A = [cos(2*pi*f(1)*t) sin(2*pi*f(1)*t) cos(2*pi*f(2)*t) sin(2*pi*f(2)*t)...
    cos(2*pi*f(3)*t) sin(2*pi*f(3)*t) cos(2*pi*f(4)*t) sin(2*pi*f(4)*t)...
    cos(2*pi*f(5)*t) sin(2*pi*f(5)*t)]; % Comes from Polar motion model


% setting up to solve for constants
constants = (A.'*A)\(A.'*r);















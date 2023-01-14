% Hw3_r11521201
close all;clear;clc;
opengl hardware; 

% 讀入txt資料---------------------------------------------------------------
file = importdata('HW3- building seismic response data.txt');
t = file(:,1); % time vector
x = file(:,2); % base acceleration
y = file(:,3); % roof acceleration

% set parameters ----------------------------------------------------------
fs = 200;  % Hz
ts = 1/fs; % sec

%%
% (a) auto-correlation 
% calculate x autocorrelation
[auto_x, lags1] = xcorr(x, length(x)/2, 'normalize');
subplot(2,2,1);
plot(lags1, auto_x)
xlabel('Point')
ylabel('Amplitude (gal)')
title('Auto-correlation Function of Base')
set(gcf,'Position',[50,50,1000,800]);

% calculate y autocorrelation
[auto_y, lags2] = xcorr(y, length(y)/2, 'normalize');
subplot(2,2,2);
plot(lags2, auto_y)
xlabel('Point')
ylabel('Amplitude (gal)')
title('Auto-correlation Function of Roof')
set(gcf,'Position',[50,50,1000,800]);


%%
% (b) cross-correlation & identify the time lag
% xy cross-correlation
[rxy, lags3] = xcorr(x, y, length(x)/2, 'unbiased');
subplot(2,2,3);
plot(lags3, rxy)

% 計算time lag
[val1, index1] = max(rxy);
x_middle = round(length(rxy)/2, 0);
delay_time_xy = (index1 - x_middle)*ts;
annotation('textbox', [0.315, 0.06, 0.1, 0.1], 'String', "time lag is " + delay_time_xy +" sec")

xlabel('Point')
ylabel('Amplitude (gal)')
title('Cross-correlation Function of Base-Roof')
set(gcf,'Position',[50,50,1000,800]);

% -------------------------------------------------------------------------
% yx cross-correlation
[ryx, lags4] = xcorr(y, x, length(x)/2, 'unbiased');
subplot(2,2,4);
plot(lags4, ryx)

% 計算time lag
[val2, index2] = max(ryx);
x_middle = round(length(ryx)/2, 0);
delay_time_yx = (index2 - x_middle)*ts;
annotation('textbox', [0.76, 0.06, 0.1, 0.1], 'String', "time lag is " + delay_time_yx+" sec")

xlabel('Point')
ylabel('Amplitude (gal)')
title('Cross-correlation Function of Roof-Base')
set(gcf,'Position',[50,50,1000,800]);
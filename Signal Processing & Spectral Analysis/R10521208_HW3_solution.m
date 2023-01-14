% 鄭棋 R10521208 訊號處理 HW3 
clc;clear;%close all;
% import data
data = importdata('HW3- building seismic response data.txt');
t = data(:,1);
ground = data(:,2);
top = data(:,3);
% plot data
figure(1);
plot(t,ground);grid on;ylim([-80,80]);xlim([0,103]);
xlabel('Time(s)');ylabel('Acceleration(gal)');
title('Aceeleration of first floor');

figure(2);
plot(t,top);grid on;ylim([-250,300]);xlim([0,103]);
xlabel('Time(s)');ylabel('Acceleration(gal)');
title('Aceeleration of top floor');

%----------------------------Set variable---------------------------------%
T = t(2)-t(1);
fs = 1/T;                                     % sampling rate (Hz)
%----------------------------Autocorrelation------------------------------%
auto_ground = xcorr(ground,length(ground)/2); % overlap half of the data
auto_top = xcorr(top,length(ground)/2); 
point2 = (0:1:length(auto_ground)-1);
point = point2 - length(ground)/2;
%----------------------------Plot result----------------------------------%
figure(3);
subplot(2,1,1);
plot(point,auto_ground);grid on;xlim([-10240,10240]);%ylim([-80,80]);
xlabel('Point');ylabel('Amplitude');
title('Autocorrelation function of ground floor acceleration');
subplot(2,1,2);
plot(point,auto_top);grid on;xlim([-10240,10240]);%ylim([-80,80]);
xlabel('Point');ylabel('Amplitude');
title('Autocorrelation function of top floor acceleration');
%----------------------------cross correlation----------------------------%
[cross,lag] = xcorr(ground,top,length(ground)/2,'unbiased'); % overlap half of the data
figure;
plot(lag,cross);grid on;xlim([-10240,10240]);%ylim([-80,80]);
xlabel('Point');ylabel('Amplitude');
title('Autocorrelation function of ground and top floor acceleration');
[maximum,index] = max(abs(cross));
fprintf('Time lag = %3.2f s\n',T*abs(lag(index)));
% practice
% m = 0.84.^(1:10);
% n = circshift(m,3);
% figure
% [c,lags] = xcorr(m,n,'biased');
% stem(lags, c);

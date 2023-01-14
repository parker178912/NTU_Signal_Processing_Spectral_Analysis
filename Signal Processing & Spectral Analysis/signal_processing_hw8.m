clc;clear;close all;
data = importdata('Ming-Li School 1999-9-21 data.xlsx');
data1 = data(:,20);

%% plot data
SR = 200; % sampling rate
t = 1/SR*(1:length(data1))';
plot(t,data1,'linewidth',0.5);
xlim([0,103]);
grid on;
xlabel('time (s)');
ylabel('accelerarion (gal)');
title('Ch20 data');

%% Hilbert transform
data2 = hilbert(data1);
figure;
plot(t,real(data2));
hold on;
plot(t,imag(data2));
hold off;
legend('Real Part','Imaginary Part');

%% Instantaneous Frequency 
tmp1 = real(data2);
tmp2 = imag(data2);
IF_data = (tmp1.*([0;diff(tmp2)])-tmp2.*([0;diff(tmp1)]))./(tmp1.^2)+tmp2.^2;
IF_data = IF_data*SR/2/pi;
figure;
plot(t,IF_data);
xlabel('time (s)');
ylabel('Instantaneous Frequency');
title('Instantaneous Frequency by Direct Differentiation');
xlim([0,103]);
hold on;
% data smoothing
IF_data_SM = sgolayfilt(IF_data,3,201);
plot(t,IF_data_SM);
legend('original','smoothed');


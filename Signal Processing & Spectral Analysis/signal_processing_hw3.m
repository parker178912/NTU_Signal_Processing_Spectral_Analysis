clc;clear;close all;

data = importdata("HW3- building seismic response data.txt");
t = data(:,1);
ground = data(:,2);
top = data(:,3);
T = 0.005;
fs = 1/T;
set(gcf,'Position',[50,50,1000,800]);
% (a)
[autoground, lag1] = xcorr(ground, length(ground)/2);
[autotop, lag2] = xcorr(top, length(top)/2);
subplot(2,2,1);
plot(lag1, autoground)
xlabel('Point')
ylabel('Ampitude')
title('Acceleration Function of Ground Floor with Added Zeros')

subplot(2,2,2);
plot(lag2, autotop)
xlabel('Point')
ylabel('Ampitude')
title('Acceleration Function of Top Floor with Added Zeros')

% (b)
[crossgt, lag3] = xcorr(ground, top, length(top)/2);
subplot(2,2,3);
plot(lag3, crossgt)
[maximum1, maxnum1] = max(crossgt);
middlenum1 = round(length(crossgt)/2);
delay1 = (maxnum1 - middlenum1)*T;
annotation('textbox', [0.36, 0.06, 0.1, 0.1], 'String', "time lag is " + delay1 +" sec")
xlabel('Point')
ylabel('Ampitude')
title('Cross-Correlation Function of Ground-to-Top Floor(Sxy)')

[crosstg, lag4] = xcorr(top, ground, length(top)/2);
subplot(2,2,4);
plot(lag4, crosstg)
[maximum2, maxnum2] = max(crosstg);
middlenum2 = round(length(crosstg)/2);
delay2 = (maxnum2 - middlenum2)*T;
annotation('textbox', [0.8, 0.06, 0.1, 0.1], 'String', "time lag is " + delay2 +" sec")
xlabel('Point')
ylabel('Ampitude')
title('Cross-Correlation Function of Top-to-Ground Floor(Syx)')
%% initialize the environment
clc
clear
opengl hardware; % use hardware to accelerate

%% read in the data
% num = value part of data
% txt = string part of data
% raw = read in the whole data
[num,txt,raw] = xlsread('Homework-2 data set-RCF-Four Specimen Test Data.xlsx');

data_title = txt(1,[3, 7]);
t = num(:, 1);      % the total time steps
y = num(:, [3. 7]);      % respectively are RCF6_060_A04 and RCF4_110_A04
y = y - mean(y);    % initialize the data

%% define the parameters
SR = round(1/(t(2) - t(1)));                              % sampling rate, 200Hz
pt = length(t);                                           % data length
sigma = 0.4;                                              % standard deviation of Gaussian distribution, it could be 0.4 or 1
winlen = 5;                                               % in sec
tmp = (-winlen:(1/SR):winlen)';                           % discrete the window we need(which is same time interval as data)
wpt = length(tmp);                                        % data length per window
wfun = sqrt(2/(sqrt(2*pi)*sigma))*exp(-tmp.^2/sigma^2);   % window function (Gaussian distribution with mean = 0)

%% do physicl spectrum
fps = cell(2,1);
ps = cell(2,1);
for i = 1:2
    [fps{i},ps{i}] = PhysicalSpectrum(y(:, i), pt, wpt, wfun, SR);
end

%% plot the result
figure
for i = 1:2
    subplot(2,1,i);
    surf(t,fps{i},ps{i},'LineStyle','none','FaceLighting','phong');
    view(0,90),colormap('jet'),xlabel('Time (sec)'),ylabel('Frequency (Hz)');
    xlim([20,90]);ylim([0,20]);
    title(['Synthetic Time series Physical Spectrum of ', data_title(i), 'with sigma = ', sigma]);
    colorbar;
    set(gcf,'Position',[50,50,1024,768]);
end

%% define function
function [fps, ps] = PhysicalSpectrum(y, pt, wpt, wfun, SR)
fps = SR*((1:wpt)'-1)/wpt;          % what is this?
ps = zeros(wpt,pt);
for I = 1:pt
    idx1 = I-(wpt-1)/2;             % index of window's left bound
    idx2 = idx1+wpt-1;              % right bound
    if (idx1 <= 0)                  % left bound out of range
        % symmetry on the edge
        tmp = [flipud(y(1:(wpt-idx2)));y(1:idx2)];
        tmp = wfun.*tmp;
    elseif (idx1 > (pt-wpt+1))      % right bound out of range
        % symmetry on the edge
        tmp = wpt-(pt-idx1+1);
        tmp = [y(idx1:end);flipud(y((pt-tmp+1):end))];
        tmp = wfun.*tmp;
    else
        % normal procedure
        tmp = wfun.*y(idx1:idx2);
    end
    ps(:,I) = abs(fft(tmp))/wpt*2;  % do  fft to the data
end
end
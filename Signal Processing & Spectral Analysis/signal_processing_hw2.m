clc;clear;close all;
[num,txt,raw] = xlsread('Homework-2 data set-RCF-Four Specimen Test Data.xlsx');
t = num(:,1);
dtitle = txt(1,[3,5,7,9]);
y = num(:,[3,5,7,9]);
% 下面兩行可調整參數
sigma = 4.0;      
winlen = 5;
%
fps = cell(4,1);
ps = cell(4,1);
for i = 1:4
    [fps{i},ps{i}] = PhysicalSpectrum(t,y(:,i),sigma,winlen);
end
figure
for i = 1:4
    subplot(2,2,i);
    surf(t,fps{i},ps{i},'LineStyle','none','FaceLighting','phong');
    view(0,90),colormap('jet'),xlabel('Time (sec)'),ylabel('Frequency (Hz)');
    xlim([20,90]);
    ylim([0,20]);
    title(['Synthetic Time series Physical Spectrum',dtitle(i)]);
    colorbar;
    set(gcf,'Position',[50,50,1024,768]);
end
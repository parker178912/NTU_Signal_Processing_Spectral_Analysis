function [coef_save,spectrogram]=MCMW_VCF_2(tt,Acc_all,sigma,dt,f_range,str)
n = length(sigma);
fs = 1/dt;
f = f_range;
type_coef = 'energy';
spectrogram=figure;
for i = 1:3
    coef_save(i) = 0;
end
for j=1:length(sigma)
    [coefs]=MCMW_VCE(Acc_all,f,fs,type_coef,sigma(j));
    % coef_save{j} = coefs;
    subplot(n,1,j)
    imagesc(tt,f,coefs)
    h = colorbar;
    colormap jet
    ax = gca;
    ax.YDir = 'normal';
    xlabel('Time(s)')
    ylabel('Frequency(Hz)')
    ylabel(h,'Amplitude')
    title(['Spectrogram of the ',str,' with \sigma=',num2str(sigma(j))])
end
end

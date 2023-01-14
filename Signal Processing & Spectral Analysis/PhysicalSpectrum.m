function [fps,ps] = PhysicalSpectrum(t,y,sigma,winlen)
y = y-mean(y);
SR = round(1/(t(2)-t(1)));
pt = length(t);
tmp = (-winlen:(1/SR):winlen)';
wpt = length(tmp);
wfun = sqrt(2/(sqrt(2*pi)*sigma))*exp(-tmp.^2/sigma^2);
fps = SR*((1:wpt)'-1)/wpt;
ps = zeros(wpt,pt);
for I = 1:pt
    idx1 = I-(wpt-1)/2;
    idx2 = idx1+wpt-1;
    if (idx1 <= 0)
        tmp = [flipud(y(1:(wpt-idx2)));y(1:idx2)];
        tmp = wfun.*tmp;
    elseif (idx1 > (pt-wpt+1))
        tmp = wpt-(pt-idx1+1);
        tmp = [y(idx1:end);flipud(y((pt-tmp+1):end))];
        tmp = wfun.*tmp;
    else
        tmp = wfun.*y(idx1:idx2);
    end
    ps(:,I) = abs(fft(tmp))/wpt*2;
end
end
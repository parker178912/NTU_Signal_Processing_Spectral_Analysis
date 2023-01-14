function [coefs]=MCMW_VCE(Acc,f,fs,type_coef,sigma)
if length(Acc)==1
    Acc = Acc';
end
scales=1;
i=(-1)^0.5;
Acc = Acc';
for s=1:length(f)
    j=1;
    for w=(0:length(Acc)-1)*fs/length(Acc)*2*pi
        psi(j) = scales^0.5*((2*pi)^0.5*sigma*exp(-0.5*(w*scales-2*pi*f(s))^2*sigma^2));
        j = j+1;
    end
    coefs(s,:) = 1/2/pi*(ifft(fft(Acc).*conj(psi)));
end
switch type_coef
    case 'complex'
        coefs = coefs;
    case 'absolute'
        coefs = abs(coefs);
    case 'energy'
        coefs = abs(coefs);
        coefs = coefs.*coefs;
end

end 
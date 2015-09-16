function y=FMRSpectrumDerivativeUniaxial(Bvalues,gtruevalues,frequency,Banvalues,K2toK1values,linewidthvalues)

    Bvalues=reshape(Bvalues,length(Bvalues),1);
    dB = diff(Bvalues);
    Bvalues2=[Bvalues-[.5*dB ; .5*dB(end)] ; Bvalues(end) + .5*dB(end)];

    w = FMRSpectrumUniaxial(Bvalues2,gtruevalues,frequency,Banvalues,K2toK1values,linewidthvalues);

    matsize=size(w);
    matsize=matsize(2:end);
    
    B = repmat(Bvalues2,[1 matsize]);
 
    y=diff(w)./diff(B);
    
end
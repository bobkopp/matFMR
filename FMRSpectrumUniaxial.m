function y=FMRSpectrumUniaxial(Bvalues,gtruevalues,frequency,Banvalues,K2toK1values,linewidthvalues)
    PlanckConst=6.62607e-34;
    BohrMagneton=9.28476e-24;
    
% following Griscom (1974)

    Btruevalues=PlanckConst*frequency./(gtruevalues*BohrMagneton);
    
    % ok, so we're going to set up our phase space for calculation of Bres
    
    % note here that we're assuming, following Griscom, that the resonance
    % field is independent of applied field (or maybe not... I have to look
    % at the derivation more closely)
    
    dtheta = 0.01;

    [nulla,nullb,Btrue,Ban,K2toK1,theta] = ndgrid (1, 1, Btruevalues,Banvalues,K2toK1values,[0:dtheta:.5*pi]);

    sinTheta = sin(theta);
    sinThetaSq = sinTheta.^2;
    
    weightings1 = 2*pi*sinTheta*dtheta;
    clear sinTheta;    

    factorA = 4*(1-sinThetaSq)-2*sinThetaSq + K2toK1.*(16*(1-sinThetaSq).*sinThetaSq-4*sinThetaSq.*sinThetaSq);

    Bres1 = (sqrt(Btrue.^2+(Ban.*factorA).^2/16)) - .25 * Ban .* factorA;
    
    clear sinThetaSq;
    
    
    Bvalues=reshape(Bvalues,[length(Bvalues) 1]);
    linewidthvalues=reshape(linewidthvalues,[1 length(linewidthvalues)]);
    
    matsize=size(Bres1);
    matsize=matsize(3:6);
    
    B = repmat(Bvalues,[1 length(linewidthvalues) matsize]);
    Bres = repmat(Bres1,[length(Bvalues) length(linewidthvalues) 1]);
    weightings = repmat(weightings1,[length(Bvalues) length(linewidthvalues) 1]);
    linewidth = repmat(linewidthvalues,[length(Bvalues) 1 matsize]);

    contribs = weightings .* sqrt(1/(2*pi)) ./linewidth .* exp (-((B - Bres).^2)./(2*linewidth.^2));

    
    % ok, now collapse along the sixth and seventh dimension
    
    absorption=sum(contribs,6);
    y=absorption;  
    
end


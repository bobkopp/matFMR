function y=FMRSpectrumCubic(Bvalues,gtruevalues,frequency,Banvalues,K2toK1values,linewidthvalues)
    PlanckConst=6.62607e-34;
    BohrMagneton=9.28476e-24;
    

    Btruevalues=PlanckConst*frequency./(gtruevalues*BohrMagneton);
    
    % ok, so we're going to set up our phase space for calculation of Bres
        
    
    dtheta = .02;
    dphi=min(max(0.005,min(linewidthvalues)*3),0.05);
    
    [nulla,nullb,Btrue,Ban,K2toK1,theta,phi] = ndgrid (1, 1, Btruevalues,Banvalues,K2toK1values,[0:dtheta:.5*pi],[0:dphi:.5*pi]);

    sinTheta = sin(theta);
    sinThetaSq = sinTheta.^2;
    sinPhiSq = sin(phi).^2;
    
    weightings1 = 4*dphi*sinTheta*dtheta;
    clear sinTheta;    
    
    A = 4*((1 - 5*((1-sinThetaSq) .* sinThetaSq + sinThetaSq.*sinThetaSq .* sinPhiSq .* (1-sinPhiSq))) - .5 * K2toK1 .* (((1-sinThetaSq) .* sinThetaSq + sinThetaSq.*sinThetaSq .* sinPhiSq .* (1-sinPhiSq)) - 21 * (sinThetaSq.*sinThetaSq .* (1-sinThetaSq) .* sinPhiSq .* (1-sinPhiSq) )));
    Bres1 = sqrt(Btrue.^2+(Ban.*A).^2/16) - .25 * (Ban.*A);
    
    clear sinThetaSq sinPhiSq;
    
    
    Bvalues=reshape(Bvalues,[length(Bvalues) 1]);
    linewidthvalues=reshape(linewidthvalues,[1 length(linewidthvalues)]);
    
    matsize=size(Bres1);
    matsize=matsize(3:7);
    
    B = repmat(Bvalues,[1 length(linewidthvalues) matsize]);
    Bres = repmat(Bres1,[length(Bvalues) length(linewidthvalues) 1]);
    weightings = repmat(weightings1,[length(Bvalues) length(linewidthvalues) 1]);
    linewidth = repmat(linewidthvalues,[length(Bvalues) 1 matsize]);

    contribs = weightings .* sqrt(1/(2*pi)) ./linewidth .* exp (-((B - Bres).^2)./(2*linewidth.^2));

    
    % ok, now collapse along the sixth and seventh dimension
    
    absorption=sum(sum(contribs,7),6);
    y=absorption;   
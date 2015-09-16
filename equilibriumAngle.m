function theta=equilibriumAngle(Bvalues,Banvalues,K2toK1values,orientThetavalues)
%    factor=.5*(Ban./(B+1e-9)).*sin(2*orientTheta);
%    offset=(factor<=orientTheta).*factor + (factor>orientTheta).*orientTheta;
%    theta=orientTheta-offset;

%    testThetavalues=[0:.01:.5*pi];
%    testTheta=reshape(testThetavalues,[ones(1,length(size(orientTheta))) length(testThetavalues)]);
%    testTheta=repmat(testTheta,[size(orientTheta) 1]);

%    Bext=repmat(B,[ones(1,length(size(orientTheta))) length(testThetavalues)])
%    Banext=repmat(Ban,[ones(1,length(size(orientTheta))) length(testThetavalues)]);
%    K2toK1ext=repmat(K2toK1,[ones(1,length(size(orientTheta))) length(testThetavalues)]);
%    orientThetaext=repmat(orientTheta,[ones(1,length(size(orientTheta))) length(testThetavalues)]);
%   
%    energy = -Bext.*(cos(orientThetaext-testTheta)) + .5*Banext.*sin(2*testTheta);%

%    [minenergy,minindex]=min(energy,[],length(size(energy)));
    
%    theta=testThetavalues(minindex);

%    clear testThetavalues testTheta Bext Banext K2toK1ext orientThetaext;
    
%    factor=.5*(Ban./(B+1e-9)).*sin(2*theta);
%    offset=(factor<=1).*asin(factor)+(factor>1).*orientTheta;
%    dtheta=orientTheta-offset-theta;
%    theta=real(theta+dtheta);
    
%    maxdtheta=max(reshape(dtheta,[prod(size(dtheta)) 1]));
    
%    while maxdtheta > 0.02
%        factor=.5*(Ban./(B+1e-9)).*sin(2*theta);
%        offset=(factor<=1).*asin(factor)+(factor>1).*orientTheta;
%        dtheta=orientTheta-offset-theta;
%        theta=real(theta+dtheta);

%        maxdtheta=max(reshape(dtheta,[prod(size(dtheta)) 1]));
%    end
   
    [B, Ban, K2toK1, orientTheta] = ndgrid(Bvalues, Banvalues, K2toK1values, orientThetavalues);

    orientTheta=orientTheta+(orientTheta==0)*1e-9;
    
    %gammalarge = 1 - (1./(2*orientTheta)) .* (Ban./(B+1e-9)) .* sin(2*orientTheta)
    %gammasmall = tan(orientTheta)./orientTheta
    %gamma=B./(B+Ban);
    gamma = orientTheta .* 1.14 .* (B./(B+Ban)).^1.12;
    gamma = (gamma>1) + (gamma<=1).*(gamma>0).*gamma;

    
    maxdgamma = 1;
    times=1;
    sumgamma=gamma;
    while maxdgamma > .02
        gamma = 1 + (1./orientTheta).*asin(-.5*Ban.*sin(2*gamma.*orientTheta)./B);
        gamma = (gamma>1) + (gamma<=1).*gamma;
        maxdgamma=max(reshape(abs(sumgamma/(times)-gamma),[prod(size(gamma)) 1]));
        times = times+1;
        sumgamma= sumgamma+gamma;
        gamma=sumgamma/times;
        if times > 1000, break; end;
    
    end
    times
    theta=real(gamma.*orientTheta);
end
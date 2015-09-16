function theta=equilibriumAngleUniaxial(Bvalues,Banvalues,K2toK1values,orientThetavalues)
    if length(size(Bvalues))==1
        [B, Ban, K2toK1, orientTheta] = ndgrid(Bvalues, Banvalues, K2toK1values, orientThetavalues);
    else
        B=Bvalues;
        Ban=Banvalues;
        K2toK1=K2toK1values;
        orientTheta=orientThetavalues;
    end
    
        
        
    h=B./(Ban+1e-9.*(Ban==0));
    matsize=size(B);
    for i=1:matsize(1)
        for j=1:matsize(2)
            for k=1:matsize(3)
                for l=1:matsize(4)
                    theta(i,j,k,l)=minimization(h(i,j,k,l),K2toK1(i,j,k,l),orientTheta(i,j,k,l));
                    %theta.appx(i,j,k,l)=largeBapproximation(h(i,j,k,l),K2toK1(i,j,k,l),orientTheta(i,j,k,l));
                    %theta.appx2(i,j,k,l)=largeBapproximation2(h(i,j,k,l),K2toK1(i,j,k,l),orientTheta(i,j,k,l));
                    %theta.tan(i,j,k,l)=tanapproximate(h(i,j,k,l),K2toK1(i,j,k,l),orientTheta(i,j,k,l));
                end
            end
        end
    end
    
    
end

function y=minimization(h,K2toK1,orientTheta)
    y=fminbnd(@(theta) -h*(cos(orientTheta-theta)) + .5*sin(theta)^2, 0, orientTheta);
end

function y=largeBapproximation(h,K2toK1,orientTheta)
    y=orientTheta-(1/(2*h))*(sin(2*orientTheta));
end


function y=largeBapproximation2(h,K2toK1,orientTheta)
    y=orientTheta-((1/(2*h))*(sin(2*orientTheta)))/(1+(1/(2*h)*cos(2*orientTheta)));
end


function y=tanapproximate(h,K2toK1,orientTheta)
    y=tanh(.6*h)*orientTheta;
end
    
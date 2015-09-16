function y=AnisotropyFieldUniaxial(K1values,Msvalues, qvalues)
    [K1,Ms,q]=ndgrid(K1values,Msvalues,qvalues);
    mu0 = 4*pi*10^-7;
    
    y=2*K1./Ms + dN(q) .* Ms * mu0;

end

% calculate demagnetizing factors following Diaz-Ricci & Kirschvink, 1992

function y = F(pvalues,qvalues)
    p=pvalues;
    q=qvalues;
    % eliminate divides by zero
    if length(find(q==0))>0
        if length(find(q~=0)>0)
            q=q+1e-9*min(q(find(q~=0)));
        else
            q=q+1e-9;
        end
        
    end

	y= (p.^2 - q.^2) .* asinh (1./sqrt(p.^2 + q.^2)) +  p .* (1 - q.^2) .* asinh(p./sqrt(1 + q.^2)) + p .* q.^2 .* asinh(p./q)  + q.^2 .* asinh (1./q) + 2 * p .* q .* atan((q./p) .* sqrt(1 + p.^2 + q.^2)) - (1/3) * (1 + p.^2 - 2*q.^2) .* sqrt(1 +p.^2 + q.^2) + (1/3)* (1 - 2*q.^2) .* sqrt (1 + q.^2) + (1/3) * (p.^2 - 2*q.^2) .* sqrt(p.^2 + q.^2) - pi * p .* q + (2/3) * q.^3;
end

function y = g(p,q)
	y = (F(p,zeros(size(q))) - F(p, q))./(p.*q);
end

function y = dN(q)
	%dNcgs = 2*pi - 6*g(1,q);
	%y = dNcgs/(4*pi);
    
    % For a prolate or oblate spheroid
    
    if q==1
        Nperp=1/3;
    elseif q>1
        Nperp = (.5*q./(q.^2-1)).*(q-.5*(log((q+sqrt(q.^2-1))./(q-sqrt(q.^2-1))))./sqrt(q.^2-1));
    elseif q<1
        m=1./q;
        Nperp = (.5./(m.^2-1)) .* ( (m.^2./sqrt(m.^2-1)) .* asin(sqrt(m.^2-1)/m) - 1);
    end
    
    y = 3 * Nperp - 1;
end
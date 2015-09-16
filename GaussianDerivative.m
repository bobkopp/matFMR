function y=GaussianDerivative(xvalues,meanvalues,stdvalues)

    xvalues=reshape(xvalues,length(xvalues),1);
    dX = diff(xvalues);
    xvalues2=[xvalues-[.5*dX ; .5*dX(end)] ; xvalues(end) + .5*dX(end)];

    [x,means,stds] = ndgrid (xvalues2, meanvalues,stdvalues);
    w = sqrt(1/(2*pi)) ./stds .* exp (-((x - means).^2)./(2*stds.^2));

    matsize=size(w);
    matsize=matsize(2:end);
    
    x = repmat(xvalues2,[1 matsize]);
 
    y=diff(w)./diff(x);
    
end
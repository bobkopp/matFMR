function y=GaussianDerivativeMulti(xvalues,varargin)
    params=varargin;
    if length(params)==1
        for i=1:length(params{1})
            params{i}=varargin{1}(i);
        end
    end
    N=floor((length(params))/3);
    y=params{1}*GaussianDerivative(xvalues,params{2},params{3});
    for i=2:N
        y=y+params{(i-1)*3+1}*GaussianDerivative(xvalues,params{(i-1)*3+2},params{(i-1)*3+3});
    end
end
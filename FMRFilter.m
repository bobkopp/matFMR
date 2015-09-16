function y=FMRFilter(a,criterion,varargin)

    y=a;
    y.fields=a.fields(criterion);
    y.datDeriv=a.datDeriv(criterion);
    y.datInt=a.datInt(criterion);
    y.datDerivSmooth=a.datDerivSmooth(criterion);
    y.datIntSmooth=a.datIntSmooth(criterion);
    y.gvalues=a.gvalues(criterion);
    y.ginverse=1./y.gvalues;
    
    if nargin>2
        for i=1:length(varargin)
            if strcmp(varargin{i},'smooth')
                y.datDeriv=y.datDerivSmooth;
                y.datInt=y.datIntSmooth;
            end
        end
    end
    

end


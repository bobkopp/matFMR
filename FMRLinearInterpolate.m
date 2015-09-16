function y=FMRLinearInterpolate(a,lowField,highField,varargin)

	doIntegration = 0;

	[m,low]=min(abs(a.fields-lowField));
	[m,high]=min(abs(a.fields-highField));
	
	if nargin>3
		for i=1:length(varargin)
			if strcmpi(varargin{i},'integrate')
				doIntegration = 1;
			end
		end
	end

    y=a;
    y.datDeriv(low+1:high-1)=interpolate(a.fields([low high]),a.datDeriv([low high]),a.fields(low+1:high-1));
    y.datDerivSmooth(low+1:high-1)=interpolate(a.fields([low high]),a.datDerivSmooth([low high]),a.fields(low+1:high-1));

    if doIntegration
    	y.datInt=cumtrapz(y.fields,y.datDeriv-sum(y.datDeriv)/length(y.datDeriv));
    	y.datIntSmooth=cumtrapz(y.fields,y.datDerivSmooth-sum(y.datDeriv)/length(y.datDeriv));
    else
	y.datInt(low+1:high-1)=interpolate(a.fields([low high]),a.datInt([low high]),a.fields(low+1:high-1));
    	y.datIntSmooth(low+1:high-1)=interpolate(a.fields([low high]),a.datIntSmooth([low high]),a.fields(low+1:high-1));
    end
    

end


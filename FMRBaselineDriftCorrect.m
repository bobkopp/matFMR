function y=FMRBaselineDriftCorrect(spectrum,varargin)

ends=0;
if nargin>1
    if strcmp(varargin{1},'highend'), ends=1; elseif strcmp(varargin{1},'lowend'),ends=-1; end
end

y=spectrum;
smoothSpan=double(int32(0.03/((max(y.fields)-min(y.fields))/length(y.fields)))*2+1);

%slopes=smooth(diff(y.datDerivSmooth)./diff(y.fields),double(int32(smoothSpan/4)*2+1),'rlowess');
%slopes=.5*([slopes(1) ; slopes] + [slopes ; slopes(end)]);

%significantpoints=find((abs(y.datIntSmooth)/max(y.datIntSmooth)>0.0005) .* (abs(slopes)>400));
%significantpoints=find((abs(y.datIntSmooth)/max(y.datIntSmooth)>0.005))
significantpoints=15:(length(y.datInt)-14);

if ends == 1, significantpoints=sort(union([1:min(significantpoints)-1],significantpoints));
elseif ends==-1, significantpoints=sort(union(significantpoints,[significantpoints+1:length(y.datIntSmooth)])); end

[linefit,linefitgoodness]=fit(y.fields([1:(min(significantpoints)-1) (max(significantpoints)+1):end]),y.datDeriv([1:(min(significantpoints)-1) (max(significantpoints)+1):end]),'poly1')

%meanBaselineSlope=mean(slopes([1:min(significantpoints) max(significantpoints):end]))

%slopeRemovedSmooth=y.datDerivSmooth-meanBaselineSlope*y.fields;
%slopeRemoved=y.datDeriv-meanBaselineSlope*y.fields;

%meanIntercept=mean(slopeRemoved([1:5]))

if linefitgoodness.rsquare > 0.9, y.datDeriv=y.datDeriv-linefit(y.fields); end
%y.datDeriv=slopeRemoved - meanIntercept;

%y.derivInterp=fit(y.fields,y.datDeriv,'linearinterp');
y.datInt=cumtrapz(y.fields,y.datDeriv);
y.datIntSmooth=smooth(y.datInt,smoothSpan,'rlowess');
%y.intInterp=fit(y.fields,y.datIntSmooth,'linearinterp');
y.datDerivSmooth=smooth(y.datDeriv,double(int32(smoothSpan/4)*2+1),'rlowess');
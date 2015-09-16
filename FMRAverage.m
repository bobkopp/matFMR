function y=FMRAverage(FMRdata)

y.sample=FMRdata(1).sample;

y.frequency=mean([FMRdata.frequency]);
y.pwr=mean([FMRdata.pwr]);
y.sweeps=mean([FMRdata.sweeps]);
y.fields=FMRdata(1).fields;
y.datDeriv=zeros(size(y.fields));
for i=1:length(FMRdata)
    y.rawSampleInterp=fit(FMRdata(i).fields,FMRdata(i).datDeriv,'linearinterp');
    y.datDeriv=y.datDeriv+(1/length(FMRdata))*y.rawSampleInterp(y.fields);
end

y.datInt=cumtrapz(y.fields,y.datDeriv);

smoothSpan=double(int32(0.03/((max(y.fields)-min(y.fields))/length(y.fields)))*2+1);
y.datIntSmooth=smooth(y.datInt,smoothSpan,'rlowess');
y.intInterp=fit(y.fields,y.datIntSmooth,'linearinterp');
y.datDerivSmooth=smooth(y.datDeriv,double(int32(smoothSpan/4)*2+1),'rlowess');
y.TotalSampleToBlank=NaN;

function y=FMRBlankSub(FMRdata, FMRblank)

fracDiffFrequency=(FMRdata.frequency-FMRblank.frequency)/FMRdata.frequency;
fracDiffPower=(FMRdata.pwr-FMRblank.pwr)/FMRdata.pwr;
diffSweeps=FMRdata.sweeps-FMRblank.sweeps;

if (abs(fracDiffFrequency)) > 0.01
    display 'WARNING: Greater than 1% difference in frequency between sample and blank.'
    FMRdata.frequency
    FMRblank.frequency
end

if (abs(fracDiffPower)) > 0.01
    display 'WARNING: Greater than 1% difference in power between sample and blank.'
    FMRdata.pwr
    FMRblank.pwr
end

if (abs(diffSweeps)>0)
    display 'WARNING: Mismatch in number of sweeps between sample and blank.'
    FMRdata.sweeps
    FMRblank.sweeps
end

y.sample=FMRdata.sample;

y.frequency=FMRdata.frequency;
y.pwr=FMRdata.pwr;
y.sweeps=FMRdata.sweeps;
y.fields=FMRdata.fields;

y.rawSampleInterp=fit(FMRdata.fields,FMRdata.datDeriv,'linearinterp');
y.rawBlankInterp=fit(FMRblank.fields,FMRblank.datDeriv,'linearinterp');

y.datDeriv=feval(y.rawSampleInterp,y.fields)-feval(y.rawBlankInterp,y.fields);
y.datInt=cumtrapz(y.fields,y.datDeriv);

smoothSpan=double(int32(0.03/((max(y.fields)-min(y.fields))/length(y.fields)))*2+1);
y.datIntSmooth=smooth(y.datInt,smoothSpan,'rlowess');
y.intInterp=fit(y.fields,y.datIntSmooth,'linearinterp');
y.datDerivSmooth=smooth(y.datDeriv,double(int32(smoothSpan/4)*2+1),'rlowess');

dataParameters=FMRQuantify(FMRdata);
blankParameters=FMRQuantify(FMRblank);

y.TotalSampleToBlank=abs(dataParameters.totalAbs/blankParameters.totalAbs);
y.MeanSampleToBlank=mean(abs(FMRdata.datDerivSmooth./FMRblank.datDerivSmooth));
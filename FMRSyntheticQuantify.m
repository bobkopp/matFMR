function y=FMRSyntheticQuantify(fields,datDeriv,frequency)

PlanckConst=6.62607e-34;
BohrMagneton=-9.28476e-24;

y.sample='synthetic';

y.fields=fields';
y.datDeriv=datDeriv;
y.frequency=frequency;

y.gvalues=-PlanckConst*y.frequency./(BohrMagneton*y.fields);
y.ginverse=1./y.gvalues;


smoothSpan=double(int32(0.03/((max(y.fields)-min(y.fields))/length(y.fields)))*2+1);
y.derivInterp=fit(y.fields,y.datDeriv,'linearinterp');
y.datInt=integrate(y.derivInterp,y.fields,max(y.fields));
y.datIntSmooth=smooth(y.datInt,smoothSpan,'rlowess');
y.intInterp=fit(y.fields,y.datIntSmooth,'linearinterp');
y.datDerivSmooth=smooth(y.datDeriv,double(int32(smoothSpan/4)*2+1),'rlowess');
y.TotalSampleToBlank=NaN;

[maxAbs,fieldIndex]=max(y.datIntSmooth);
y.maxAbsUnsmooth=max(y.datInt);
y.maxDeriv=max(y.datDerivSmooth);
y.maxDerivUnsmooth=max(y.datDeriv);
y.totalAbs=integrate(y.intInterp,max(y.fields),0);

y.Beff=y.fields(fieldIndex);
y.geff=-PlanckConst*y.frequency/(BohrMagneton*y.Beff);

topHalfFields=y.fields(fieldIndex+1:length(y.fields));
topHalfAbs=y.datIntSmooth(fieldIndex+1:length(y.fields));

lowHalfFields=y.fields(1:fieldIndex);
lowHalfAbs=y.datIntSmooth(1:fieldIndex);

distFromHalf=abs(y.datIntSmooth/maxAbs -.5);
topHalfDist=distFromHalf(fieldIndex+1:length(y.fields));
lowHalfDist=distFromHalf(1:fieldIndex);

[tmpA,tmpB]=min(topHalfDist);
[tmpC,tmpD]=min(lowHalfDist);

y.dBlow=y.Beff-lowHalfFields(tmpD);
y.dBhigh=topHalfFields(tmpB)-y.Beff;
y.A=y.dBhigh/y.dBlow;
y.dBFWHM=y.dBhigh+y.dBlow;

[y.derivMax,MaxIndex]=max(y.datDeriv);
[y.derivMin,MinIndex]=min(y.datDeriv);
y.derivMaxField=y.fields(MaxIndex);
y.derivMinField=y.fields(MinIndex);
y.dBpeaktopeak=y.derivMinField-y.derivMaxField;

smoothParams=SmoothedLinearBaselineRemove(y.fields,y.datDeriv,[0.2 0.3]);
y.derivPeakRatio200To300=smoothParams.AreaUnderToPeakMax;

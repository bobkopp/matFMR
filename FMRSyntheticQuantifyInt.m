function y=FMRSyntheticQuantifyInt(fields,datInt,frequency)

PlanckConst=6.62607e-34;
BohrMagneton=-9.28476e-24;

y.sample='synthetic';

y.fields=fields';
y.datInt=datInt;
y.datIntSmooth=y.datInt;
y.intInterp=fit(y.fields,y.datInt,'linearinterp')
y.datDeriv=differentiate(y.intInterp,y.fields);
y.datDerivSmooth=y.datDeriv;
y.frequency=frequency;

y.gvalues=-PlanckConst*y.frequency./(BohrMagneton*y.fields);
y.ginverse=1./y.gvalues;

y.TotalSampleToBlank=NaN;

[maxAbs,fieldIndex]=max(y.datInt);
y.maxAbsUnsmooth=max(y.datInt);
y.maxDeriv=max(y.datDeriv);
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


y.factorAlpha=.17*y.A + .98*y.dBFWHM;
y.factorBeta=.98*y.A - .17*y.dBFWHM;


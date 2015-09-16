function y=FMRQuantify(FMRdata)

if length(FMRdata) == 1
	
	PlanckConst=6.62607e-34;
	BohrMagneton=-9.28476e-24;
	
	[maxAbs,fieldIndex]=max(FMRdata.datIntSmooth);
	
	y.sample=FMRdata.sample;
	
	y.maxAbs=maxAbs;
	y.maxAbsUnsmooth=max(FMRdata.datInt);
	y.maxDeriv=max(FMRdata.datDerivSmooth);
	y.maxDerivUnsmooth=max(FMRdata.datDeriv);
	%y.totalAbs=integrate(FMRdata.intInterp,max(FMRdata.fields),0);
	y.totalAbs=trapz(FMRdata.fields,FMRdata.datInt);
	
	sub=find((FMRdata.fields>.2).*(FMRdata.fields<.5));
	y.abs200to500 = trapz(FMRdata.fields(sub),FMRdata.datInt(sub));
	y.fracabs200to500 = y.abs200to500/y.totalAbs;
	
	y.TotalSampleToBlank=FMRdata.TotalSampleToBlank;
	y.derivSumSquares=norm(FMRdata.datDeriv)^2;
	
	
	y.Beff=FMRdata.fields(fieldIndex);
	y.geff=-PlanckConst*FMRdata.frequency/(BohrMagneton*y.Beff);
	
	topHalfFields=FMRdata.fields(fieldIndex+1:length(FMRdata.fields));
	topHalfAbs=FMRdata.datIntSmooth(fieldIndex+1:length(FMRdata.fields));
	
	lowHalfFields=FMRdata.fields(1:fieldIndex);
	lowHalfAbs=FMRdata.datIntSmooth(1:fieldIndex);
	
	distFromHalf=abs(FMRdata.datIntSmooth/maxAbs -.5);
	topHalfDist=distFromHalf(fieldIndex+1:length(FMRdata.fields));
	lowHalfDist=distFromHalf(1:fieldIndex);
	
	[tmpA,tmpB]=min(topHalfDist);
	[tmpC,tmpD]=min(lowHalfDist);

	y.dBlow=y.Beff-lowHalfFields(tmpD);
	y.dBhigh=topHalfFields(tmpB)-y.Beff;
	y.A=y.dBhigh/y.dBlow;
	y.dBFWHM=y.dBhigh+y.dBlow;
	
	
	subset=find((FMRdata.fields>(y.Beff-y.dBlow)).*(FMRdata.fields<(y.Beff+y.dBhigh)));
	y.absHalfMaximumPeak=trapz(FMRdata.fields(subset),FMRdata.datInt(subset));
	
	[y.derivMax,MaxIndex]=max(FMRdata.datDeriv);
	[y.derivMin,MinIndex]=min(FMRdata.datDeriv);
	y.derivMaxField=FMRdata.fields(MaxIndex);
	y.derivMinField=FMRdata.fields(MinIndex);
	y.dBpeaktopeak=y.derivMinField-y.derivMaxField;
	
	[y.derivMaxSmooth,MaxSmoothIndex]=max(FMRdata.datDerivSmooth);
	[y.derivMinSmooth,MinSmoothIndex]=min(FMRdata.datDerivSmooth);
	y.derivMaxSmoothField=FMRdata.fields(MaxSmoothIndex);
	y.derivMinSmoothField=FMRdata.fields(MinSmoothIndex);
	
	%y.midslope=(y.derivMax-y.derivMin)/(y.derivMax*(y.derivMaxField-y.derivMinField));
	%y.sloperatio=-((y.maxAbsUnsmooth-FMRdata.datInt(MaxIndex))/(y.Beff-y.derivMaxField))/((y.maxAbs-FMRdata.datInt(MinIndex))/(y.Beff-y.derivMinField));
	
	%y.midslopeSmooth=(y.derivMaxSmooth-y.derivMinSmooth)/(y.derivMaxSmooth*(y.derivMaxSmoothField-y.derivMinSmoothField));
	%y.sloperatioSmooth=-((y.maxAbs-FMRdata.datIntSmooth(MaxSmoothIndex))/(y.Beff-y.derivMaxSmoothField))/((y.maxAbs-FMRdata.datIntSmooth(MinSmoothIndex))/(y.Beff-y.derivMinSmoothField));
	
	%y.factorC = 4 - 12.5 * y.dBFWHM - 3 * y.sloperatioSmooth;
	%y.factorE= 12.5 * y.dBFWHM - 3 * y. sloperatioSmooth;
	
	
	cumulativeDist=cumtrapz(FMRdata.fields,FMRdata.datInt)/trapz(FMRdata.fields,FMRdata.datInt);
	y.B50=interpolate(cumulativeDist,FMRdata.fields,.5);
	y.dB67=interpolate(cumulativeDist,FMRdata.fields,.5+.333)-interpolate(cumulativeDist,FMRdata.fields,.5-.333);
	y.dB95=interpolate(cumulativeDist,FMRdata.fields,.975)-interpolate(cumulativeDist,FMRdata.fields,.025);
	y.dB99=interpolate(cumulativeDist,FMRdata.fields,.995)-interpolate(cumulativeDist,FMRdata.fields,.005);
	
	y.meanfield=trapz(FMRdata.fields,FMRdata.datInt.*FMRdata.fields)/trapz(FMRdata.fields,FMRdata.datInt);
	y.dispersion=sqrt(trapz(FMRdata.fields,((FMRdata.fields-y.meanfield).^2).*FMRdata.datInt)/trapz(FMRdata.fields,FMRdata.datInt));
	y.skewness=trapz(FMRdata.fields,((FMRdata.fields-y.meanfield).^3).*FMRdata.datInt)/(y.dispersion^3);
	
	y.factorAlpha=.17*y.A + .98*y.dBFWHM;
	y.factorBeta=.98*y.A - .17*y.dBFWHM;
	
	y.factorAlphaPrime=0.0365*y.A + 0.2107 * y.dBFWHM - 0.01622;
	
	y.BaselineSlope=[FMRdata.datDerivSmooth(end)-FMRdata.datDerivSmooth(1)]./[FMRdata.fields(end)-FMRdata.fields(1)];
else
	for i=1:length(FMRdata)
		y(i) = FMRQuantify(FMRdata(i));
	end
end

end
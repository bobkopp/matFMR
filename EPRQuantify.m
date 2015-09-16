function y=EPRQuantify(FMRdata)

if length(FMRdata)==1
	
	dDeltaDerivSmooth=[0 ; diff(FMRdata.datDeriv-FMRdata.datDerivSmooth)];
	dDeltaSign=(dDeltaDerivSmooth>0) - (dDeltaDerivSmooth<0); 
	y.totalAbs=trapz(FMRdata.fields,cumtrapz(FMRdata.fields,FMRdata.datDeriv));
	
	B0=FMRdata.fieldtogfactor/2;
	Bcenter=B0+.025;
	dB=.005;
	
	s=find((FMRdata.fields>Bcenter-dB-.003).*(FMRdata.fields<Bcenter+dB));
	[t1,i1]=max(FMRdata.datDeriv(s));
	[t2,i2]=min(FMRdata.datDeriv(s));
	
	y.MnII=t1-t2;
	y.MnIIdB=FMRdata.fields(s(i2))-FMRdata.fields(s(i1));
	
	if y.MnIIdB>2*dB
		[t1,i1]=max(FMRdata.datDeriv(s)-FMRdata.datDerivSmooth(s));
		[t2,i2]=min(FMRdata.datDeriv(s)-FMRdata.datDerivSmooth(s));
		y.MnII=t1-t2;
		y.MnIIdB=FMRdata.fields(s(i2))-FMRdata.fields(s(i1));
	end
	
	if or(or(dDeltaSign(s(i1))<0,dDeltaSign(s(i1)+1)>0),or(dDeltaSign(s(i2))>0,dDeltaSign(s(i2)+1)<0))
		y.MnII=0; y.MnIIdB=0;
	end
	
	Bcenter=FMRdata.fieldtogfactor/4.3;
	dB=.03;
	
	s=find((FMRdata.fields>Bcenter-dB).*(FMRdata.fields<Bcenter+dB));
	b=FMRdata.datDeriv(s(1));
	m=(FMRdata.datDeriv(s(end))-FMRdata.datDeriv(s(1)))/((FMRdata.fields(s(end))-FMRdata.fields(s(1))));
	n=FMRdata.datDeriv(s)-b-m*(FMRdata.fields(s)-FMRdata.fields(s(1)));
	[t1,i1]=max(n);
	y.FeIII=t1;
	
	Bcenter=FMRdata.fieldtogfactor/2;
	dB=.0015;
	
	s=find((FMRdata.fields>Bcenter-dB).*(FMRdata.fields<Bcenter+dB));
	[t1,i1]=max(FMRdata.datDeriv(s));
	[t2,i2]=min(FMRdata.datDeriv(s));
	y.FRC=t1-t2;
	y.FRCdB=FMRdata.fields(s(i2))-FMRdata.fields(s(i1));
	
	if or(y.MnIIdB>1.9*dB,y.MnIIdB==0)
		[t1,i1]=max(FMRdata.datDeriv(s)-FMRdata.datDerivSmooth(s));
		[t2,i2]=min(FMRdata.datDeriv(s)-FMRdata.datDerivSmooth(s));
		y.FRC=t1-t2;
		y.FRCdB=FMRdata.fields(s(i2))-FMRdata.fields(s(i1));
	end
	
	if or(or(dDeltaSign(s(i1))<0,dDeltaSign(s(i1)+1)>0),or(dDeltaSign(s(i2))>0,dDeltaSign(s(i2)+1)<0))
		y.FRC=0; y.FRCdB=0;
	end
else
	for i=1:length(FMRdata)
		y(i) = EPRQuantify(FMRdata(i));
	end

end

end

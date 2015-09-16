function y=FMRImport2(filetoload)

%
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Oct 10 09:03:39 EDT 2012

PlanckConst=6.62607e-34;
BohrMagneton=-9.28476e-24;

periodsInName=regexpi(filetoload,'\.');
extensionPeriod=periodsInName(length(periodsInName));

y.sample=filetoload(1:extensionPeriod-1);
%parameterFile=[y.sample,'.par'];

[y.fields,y.datDeriv,info] = BrukerRead([y.sample,'.DTA']);
y.pwr = info.MWPW;
y.frequency = info.MWFQ;
y.sweeps = info.NbScansDone;

y.fieldtogfactor=(-PlanckConst*y.frequency/BohrMagneton);

if (y.frequency<9e9) || (y.frequency > 1e10)
    warningString = ['Warning! Sample ' y.sample ' has an unusual frequency of ' num2str(y.frequency) '!']
end


%[y.fields,y.datDeriv]=textread(filetoload,'%f %f','headerlines',3);
y.fields=y.fields*1e-4;
if y.fields(1)==0
    y.fields(1)=5e-5;
end
y.gvalues=-PlanckConst*y.frequency./(BohrMagneton*y.fields);
y.ginverse=1./y.gvalues;


smoothSpan=double(int32(0.03/((max(y.fields)-min(y.fields))/length(y.fields)))*2+1);
%y.derivInterp=fit(y.fields,y.datDeriv,'linearinterp');
y.datInt=cumtrapz(y.fields,y.datDeriv);
y.datIntSmooth=smooth(y.datInt,smoothSpan,'rlowess');
%y.intInterp=fit(y.fields,y.datIntSmooth,'linearinterp');
y.datDerivSmooth=smooth(y.datDeriv,double(int32(smoothSpan/4)*2+1),'rlowess');

y.TotalSampleToBlank=NaN;
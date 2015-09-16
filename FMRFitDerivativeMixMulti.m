function  y=FMRFitDerivativeMixMulti(a,varargin)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Dec 27 16:46:12 EST 2012

PlanckConst=6.62607e-34;
BohrMagneton=-9.28476e-24;

for i=1:length(a)
	y.aQuant=FMRQuantify(a(i));
	y.totalAbs(i)=y.aQuant.totalAbs;
	y.totalAbsNorm(i)=1;
	y.geff(i)=y.aQuant.geff;
	y.A(i)=y.aQuant.A;
	y.dBFWHM(i)=y.aQuant.dBFWHM;
end

numspectra=size(a,2);
fieldsMulti=[a.fields];
derivMulti=[a.datDeriv];
%normFactors=sqrt(sum(derivMulti.*derivMulti));
normFactors=repmat(1,1,numspectra);
normFactors=y.totalAbs;
y.normFactors=normFactors;
derivMultiNorm = derivMulti./repmat(normFactors,length(fieldsMulti(:,1)),1);

fields=fieldsMulti(1:2:end,1);
deriv=derivMultiNorm(1:2:end,:);

y.fittype='mixed-derivative-multi';

y.frequency=mean([a.frequency]);
y.fieldtogfactor=(-PlanckConst*y.frequency/BohrMagneton);



y.components=0;
y.fitequation='0';
y.modelmatrix='';
y.coeffs={};
y.startingPoint=[];
y.Lower=[];
y.Upper=[];
y.secondOrder=[];

y.a=[];
y.g=[];
y.Ban=[];
y.K2toK1=[];
y.lw=[];

y.aparams=[];
y.gparams=[];
y.Banparams=[];
y.K2toK1params=[];
y.lwparams=[];

nextparam=1;

for i=1:2:length(varargin)
	if strcmp(varargin(i),'uniaxial')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'uniaxial2nd')
		y.components=y.components+1;
	elseif strcmp(varargin(i),'fixeduniaxial')
		y.components=y.components+1;        
    elseif strcmp(varargin(i),'uniaxial-fixedBan')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'uniaxial-fixedg')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'cubic')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'cubic2nd')
		y.components=y.components+1;
	elseif strcmp(varargin(i),'fixedcubicl')
		y.components=y.components+1;        
    elseif strcmp(varargin(i),'cubic-fixedBan')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'cubic-fixedg')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'gaussian')
		y.components=y.components+1;
	elseif strcmp(varargin(i),'fixedgaussian')
		y.components=y.components+1;
    elseif strcmp(varargin(i),'startingAmplitudes')
        startingAmplitudes=varargin{i+1};
        startingAmplitudes = bsxfun(@rdivide,startingAmplitudes,y.normFactors);
        y.startingPoint=startingAmplitudes(:)';
    end
end

y.aparams=1:(y.components*numspectra);
y.amatrix=[];
for i=1:y.components
	for j=1:numspectra
		y.amatrix=[y.amatrix ' x(' num2str((j-1)*y.components+i) ') '];
	end
	y.amatrix = [y.amatrix ' ; '];
end


nextparam=max(y.aparams)+1;

if length(y.startingPoint)==0
    y.startingPoint=reshape(repmat([y.totalAbsNorm],y.components,1)*.2/y.components,1,y.components*numspectra);
end
    
y.Lower=repmat(0,1,y.components*numspectra);
y.Upper=repmat(Inf,1,y.components*numspectra);

j=1;
for i=1:2:length(varargin)
	if strcmp(varargin(i),'uniaxial')
        
        y.fittypes{j}='uniaxial';
        y.fixed(j)=0;
        y.secondOrder(j)=0;
        
        
        paramspec=varargin{i+1};
        startpoint=paramspec(1,:);

        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeUniaxial(fields,x(' num2str(nextparam) '),' num2str(y.frequency) ',x(' num2str(nextparam+1) ,'),0,x(' num2str(nextparam+2) ') )'];
        
        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams nextparam+1];
        y.lwparams=[y.lwparams nextparam+2];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+3;
        
        y.startingPoint=[y.startingPoint startpoint];
        y.Lower=[y.Lower 1 -1 0.001];
        y.Upper=[y.Upper 10 1 1];
        
        y.g=[y.g y.startingPoint(end-2)];
        y.Ban=[y.Ban y.startingPoint(end-1)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];
        
    elseif strcmp(varargin(i),'uniaxial2nd')
        
        y.fittypes{j}='uniaxial';
        y.fixed(j)=0;
        y.secondOrder(j)=1;
        
        
        paramspec=varargin{i+1};
        startpoint=paramspec(1,:);

        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeUniaxial(fields,x(' num2str(nextparam) '),' num2str(y.frequency) ',x(' num2str(nextparam+1) ,'),x(',num2str(nextparam+3) ,'),x(' num2str(nextparam+2) ') )'];
        
        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams nextparam+1];
        y.lwparams=[y.lwparams nextparam+2];
        y.K2toK1params=[y.K2toK1params nextparam+3];
        nextparam=nextparam+4;
        
        
        y.startingPoint=[y.startingPoint startpoint];
        y.Lower=[y.Lower 1 -1 0.001 -10];
        y.Upper=[y.Upper 10 1 1 10];
        
        y.g=[y.g y.startingPoint(end-3)];
        y.Ban=[y.Ban y.startingPoint(end-2)];
        y.K2toK1=[y.K2toK1 y.startingPoint(end)];
        y.lw=[y.lw y.startingPoint(end-1)];

    elseif strcmp(varargin(i),'fixeduniaxial')
        
        fixedvalues=varargin{i+1};
        
        y.fittypes{j}='uniaxial';
        y.secondOrder(j)=1;
        y.fixed(j)=1;
        
        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeUniaxial(fields,'  num2str(fixedvalues(1)) ',' num2str(y.frequency) ',' num2str(fixedvalues(2)) ',' num2str(fixedvalues(3)) ',' num2str(fixedvalues(4)) ') '];

        y.gparams=[y.gparams NaN];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams NaN];
        y.K2toK1params=[y.K2toK1params NaN];

        y.g=[y.g fixedvalues(1)];
        y.Ban=[y.Ban fixedvalues(2)];
        y.K2toK1=[y.K2toK1 fixedvalues(3)];
        y.lw=[y.lw fixedvalues(4)];

    elseif strcmp(varargin(i),'uniaxial-fixedBan')
        
        startpoint=varargin{i+1};
        
        y.fittypes{j}='uniaxial';
        y.secondOrder(j)=0;
        y.fixed(j)=0;
        
        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeUniaxial(fields,x(' num2str(nextparam) '),' num2str(y.frequency) ',' num2str(startpoint(2)) ,',0,x(' num2str(nextparam+1) ') )'];

        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams nextparam+1];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+2;

        y.startingPoint=[y.startingPoint startpoint([1 3])];
        y.Lower=[y.Lower 1 0.001];
        y.Upper=[y.Upper 10 1];
        
        y.g=[y.g y.startingPoint(end-1)];
        y.Ban=[y.Ban startpoint(2)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];
        
        
    elseif strcmp(varargin(i),'uniaxial-fixedg')
        
        y.fittypes{j}='uniaxial';
        y.fixed(j)=0;
        y.secondOrder(j)=0;
        
        
        paramspec=varargin{i+1};
        startpoint=paramspec(1,:);

        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeUniaxial(fields,' num2str(startpoint(1)) ',' num2str(y.frequency) ',x(' num2str(nextparam) ,'),0,x(' num2str(nextparam+1) ') )'];
        
        y.gparams=[y.gparams NaN];
        y.Banparams=[y.Banparams nextparam];
        y.lwparams=[y.lwparams nextparam+1];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+2;
        
        
        y.startingPoint=[y.startingPoint startpoint(2:3)];
        y.Lower=[y.Lower -1 0.001];
        y.Upper=[y.Upper 1 1];
        
        y.g=[y.g startpoint(1)];
        y.Ban=[y.Ban y.startingPoint(end-1)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];

        
    elseif strcmp(varargin(i),'cubic')
        
        y.fittypes{j}='cubic';
        y.fixed(j)=0;
        y.secondOrder(j)=0;
        
        
        paramspec=varargin{i+1};
        startpoint=paramspec(1,:);

        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeCubic(fields,x(' num2str(nextparam) '),' num2str(y.frequency) ',x(' num2str(nextparam+1) ,'),0,x(' num2str(nextparam+2) ') )'];
        
        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams nextparam+1];
        y.lwparams=[y.lwparams nextparam+2];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+3;
        
        
        y.startingPoint=[y.startingPoint startpoint];
        y.Lower=[y.Lower 1 -1 0.001];
        y.Upper=[y.Upper 10 1 1];
        
        y.g=[y.g y.startingPoint(end-2)];
        y.Ban=[y.Ban y.startingPoint(end-1)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];
        
    elseif strcmp(varargin(i),'fixedcubic')
        
        fixedvalues=varargin{i+1};
        
        y.fittypes{j}='cubic';
        y.secondOrder(j)=1;
        y.fixed(j)=1;
        
        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeCubic(fields,'  num2str(fixedvalues(1)) ',' num2str(y.frequency) ',' num2str(fixedvalues(2)) ',' num2str(fixedvalues(3)) ',' num2str(fixedvalues(4)) ') '];

        y.gparams=[y.gparams NaN];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams NaN];
        y.K2toK1params=[y.K2toK1params NaN];

        y.g=[y.g fixedvalues(1)];
        y.Ban=[y.Ban fixedvalues(2)];
        y.K2toK1=[y.K2toK1 fixedvalues(3)];
        y.lw=[y.lw fixedvalues(4)];

    elseif strcmp(varargin(i),'cubic-fixedBan')
        
        startpoint=varargin{i+1};
        
        y.fittypes{j}='cubic';
        y.secondOrder(j)=0;
        y.fixed(j)=0;
        
        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeCubic(fields,x(' num2str(nextparam) '),' num2str(y.frequency) ',' num2str(startpoint(2)) ,',0,x(' num2str(nextparam+1) ') )'];

        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams nextparam+1];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+2;

        y.startingPoint=[y.startingPoint startpoint([1 3])];
        y.Lower=[y.Lower 1 0.001];
        y.Upper=[y.Upper 10 1];
        
        y.g=[y.g y.startingPoint(end-1)];
        y.Ban=[y.Ban startpoint(2)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];
        
        
    elseif strcmp(varargin(i),'cubic-fixedg')
        
        y.fittypes{j}='cubic';
        y.fixed(j)=0;
        y.secondOrder(j)=0;
        
        
        paramspec=varargin{i+1};
        startpoint=paramspec(1,:);

        y.modelmatrix=[y.modelmatrix ' FMRSpectrumDerivativeCubic(fields,' num2str(startpoint(1)) ',' num2str(y.frequency) ',x(' num2str(nextparam) ,'),0,x(' num2str(nextparam+1) ') )'];
        
        y.gparams=[y.gparams NaN];
        y.Banparams=[y.Banparams nextparam];
        y.lwparams=[y.lwparams nextparam+1];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+2;
        
        
        y.startingPoint=[y.startingPoint startpoint(2:3)];
        y.Lower=[y.Lower -1 0.001];
        y.Upper=[y.Upper 1 1];
        
        y.g=[y.g startpoint(1)];
        y.Ban=[y.Ban y.startingPoint(end-1)];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];

        
     elseif strcmp(varargin(i),'gaussian')
        y.fittypes{j}='gaussian';
        y.fixed(j)=0;
        y.secondOrder(j)=0;
        y.modelmatrix=[y.modelmatrix ' (1./(x(' num2str(nextparam+1) ').*sqrt(2*pi))).*exp(-((fields-'  num2str(y.fieldtogfactor) './x(' num2str(nextparam) '))./x(' num2str(nextparam+1) ')).^2).*-2.*((fields-'  num2str(y.fieldtogfactor) './x(' num2str(nextparam) '))./x(' num2str(nextparam+1) '))'];

        
        y.gparams=[y.gparams nextparam];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams nextparam+1];
        y.K2toK1params=[y.K2toK1params NaN];
        nextparam=nextparam+2;
        
        startpoint=varargin{i+1};
        
        y.startingPoint=[y.startingPoint startpoint];
        y.Lower=[y.Lower 1 0.001];
        y.Upper=[y.Upper 10 1];
        
        y.g=[y.g y.startingPoint(end-1)];
        y.Ban=[y.Ban 0];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw y.startingPoint(end)];
       
    elseif strcmp(varargin(i),'fixedgaussian')
                fixedvalues=varargin{i+1};

        y.fittypes{j}='gaussian';
        y.fixed(j)=1;
        y.secondOrder(j)=0;
        y.modelmatrix=[y.modelmatrix ' (1/(' num2str((fixedvalues(2))) '*sqrt(2*pi))).*exp(-((fields-'  num2str(y.fieldtogfactor) '/' num2str(fixedvalues(1)) ')/(' num2str(fixedvalues(2)) ')).^2)*-2*((fields-'  num2str(y.fieldtogfactor) '/(' num2str(fixedvalues(1)) '))/(' num2str(fixedvalues(2)) '))'];

 
        y.gparams=[y.gparams NaN];
        y.Banparams=[y.Banparams NaN];
        y.lwparams=[y.lwparams NaN];
        y.K2toK1params=[y.K2toK1params NaN];

        y.g=[y.g fixedvalues(1)];
        y.Ban=[y.Ban 0];
        y.K2toK1=[y.K2toK1 0];
        y.lw=[y.lw fixedvalues(2)];
    end
    j=j+1;
end


y.equation = ['[' y.modelmatrix '] * [ ' y.amatrix ' ] '];
y.equation

fitoptions=optimset('Display','iter','TolX',1e-4,'MaxFunEvals',4000);
[y.coeffs,y.resnorm,y.residual,y.exitflag,y.output,y.lambda,y.jacobian]=lsqnonlin(@(x)reshape(eval(y.equation)-deriv,[],1),y.startingPoint,y.Lower,y.Upper,fitoptions);

R=qr(y.jacobian,0);
rinv = R \ eye(length(R));
dfe = prod(size(deriv)) - (nextparam - 1);
sse = y.resnorm;
v=sum(rinv.^2,2) * (sse / dfe);
alpha = .05/2;
t = -tinv(alpha,dfe);
activebounds=[];
y.db = t * sqrt(v');

y.a=[];

for i=1:y.components
	for j=1:numspectra
		y.a(i,j) = y.coeffs((j-1)*y.components+i) * y.normFactors(j);
		y.aerror(i,j) = y.normFactors(j) * y.db((j-1)*y.components+i);
	end
end

for i=1:y.components
	if isfinite(y.gparams(i))
		y.g(i)=y.coeffs(y.gparams(i));
		y.gerror(i)=y.db(y.gparams(i));
	else
		%y.g(i)=NaN;
		y.gerror(i)=0;
	end
	
	if isfinite(y.Banparams(i))
		y.Ban(i)=y.coeffs(y.Banparams(i));
		y.Banerror(i)=y.db(y.Banparams(i));
	else
		%y.Ban(i)=NaN;
		y.Banerror(i)=0;
	end

	if isfinite(y.K2toK1params(i))
		y.K2toK1(i)=y.coeffs(y.K2toK1params(i));
		y.K2toK1error(i)=y.db(y.K2toK1params(i));
	else
		%y.K2toK1(i)=0;
		y.K2toK1error(i)=0;
	end
	
	if isfinite(y.lwparams(i))
		y.lw(i)=y.coeffs(y.lwparams(i));
		y.lwerror(i)=y.db(y.lwparams(i));
	else
		%y.lw(i)=NaN;
		y.lwerror(i)=0;
	end
end


fields=fieldsMulti(:,1);
deriv = derivMulti;

for i=1:y.components
	for j=1:numspectra
        if strcmp(y.fittypes{i},'uniaxial')
            y.comp(:,j,i)=y.a(i,j)*FMRSpectrumDerivativeUniaxial(fields,y.g(i),y.frequency,y.Ban(i),y.K2toK1(i),y.lw(i));        
            y.compInt(:,j,i)=cumtrapz(fields,y.comp(:,j,i));
            y.compArea(i,j)=trapz(fields,y.compInt(:,j,i));
            y.compSumSquares(i,j)=sum(y.comp(:,j,i).*y.comp(:,j,i));
            y.compIntSumSquares(i,j)=sum(y.compInt(:,j,i).*y.compInt(:,j,i))^2;
        elseif strcmp(y.fittypes{i},'cubic')
            y.comp(:,j,i)=y.a(i,j)*FMRSpectrumDerivativeCubic(fields,y.g(i),y.frequency,y.Ban(i),y.K2toK1(i),y.lw(i));        
            y.compInt(:,j,i)=cumtrapz(fields,y.comp(:,j,i));
            y.compArea(i,j)=trapz(fields,y.compInt(:,j,i));
            y.compSumSquares(i,j)=sum(y.comp(:,j,i).*y.comp(:,j,i));
            y.compIntSumSquares(i,j)=sum(y.compInt(:,j,i).*y.compInt(:,j,i))^2;
        elseif strcmp(y.fittypes{i},'gaussian')
            y.comp(:,j,i)=y.a(i,j).*(1./(y.lw(i).*sqrt(2*pi))).*exp(-((fields-y.fieldtogfactor./y.g(i))./y.lw(i)).^2)*-2.*((fields-y.fieldtogfactor./y.g(i))./y.lw(i));        
            y.compInt(:,j,i)=cumtrapz(fields,y.comp(:,j,i));
            y.compArea(i,j)=trapz(fields,y.compInt(:,j,i));
            y.compSumSquares(i,j)=sum(y.comp(:,j,i).*y.comp(:,j,i));
            y.compIntSumSquares(i,j)=sum(y.compInt(:,j,i).*y.compInt(:,j,i))^2;

        end
    end
end

y.fitValues=sum(y.comp,3);
y.fitIntValues=cumtrapz(fields,y.fitValues);

y.sse=sum((y.fitValues-deriv).*(y.fitValues-deriv),1);
y.sseInt=sum((y.fitIntValues-deriv).*(y.fitIntValues-deriv),1);

y.fitArea=sum(y.compArea,1);
y.f=y.compArea./repmat(y.fitArea,y.components,1);

y.DerivTotalSumSquares=sum(deriv.*deriv,1);
y.DerivFitSumSquares=sum(y.fitValues.*y.fitValues,1);

y.Nparameters = nextparam - 1;
y.dfe = dfe;
y.n = prod(size(deriv));
y.AIC = y.n .* log(y.resnorm/y.n) + 2 * y.Nparameters;

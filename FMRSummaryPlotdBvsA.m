function y=FMRSummaryPlotdBvsA(datasets,varargin);

colors = 'brgmcky';
symbols = 'osdv^<>h';

labels='';

if nargin > 1
    labels = varargin{1};
end

if nargin > 2
    colors = varargin{2};
end

if nargin > 3
    symbols = varargin{3};
end

axes1=axes('FontSize',12);

hold on;

axis tight;

miny=50;
maxy=210;

minx=1;
maxx=1.3;

for i=1:length(datasets)
    clear quant;
    length(datasets{i})
    for j=1:length(datasets{i})
        quant(j)=FMRQuantify(datasets{i}(j));
    end
    maxx=max([maxx quant.A]);
    maxy=1000*max([maxy/1000 quant.dBFWHM]);
    minx=min([minx quant.A]);
    miny=1000*min([miny/1000 quant.dBFWHM]);
    
    plot([quant.A],[quant.dBFWHM]*1000,[colors(mod(i-1,length(colors))+1) symbols(mod(i-1,length(symbols))+1)],'MarkerFaceColor',colors(mod(i-1,length(colors))+1));
    
end

if length(labels)>0, legend(labels); end;


A=[0 2];
plot(A,(.25-.17*A)/.98e-3,'k--');
plot(A,(.30-.17*A)/.98e-3,'k--');
plot(A,(.35-.17*A)/.98e-3,'k--');
plot(A,(.40-.17*A)/.98e-3,'k--');

text(1.1,70,'\alpha = 0.25');
text(1.172,110,'\alpha = 0.30');
text(1.213,150,'\alpha = 0.35');
text(1.277,195,'\alpha = 0.40');

xlabel(axes1,'A');
ylabel(axes1,'\DeltaB_{FWHM} (mT)');

xlim([.95*minx 1.05*maxx]);
ylim([.95*miny 1.05*maxy]);

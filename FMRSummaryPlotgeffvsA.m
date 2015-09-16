function y=FMRSummaryPlotgeffvsA(datasets,varargin);

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

for i=1:length(datasets)
    clear quant;
    for j=1:length(datasets{i})
        quant(j)=FMRQuantify(datasets{i}(j));
    end
    plot([quant.A],[quant.geff],[colors(mod(i-1,length(colors))+1) symbols(mod(i-1,length(symbols))+1)],'MarkerFaceColor',colors(mod(i-1,length(colors))+1));
    
end

if length(labels)>0, legend(labels); end;

axis manual;

plot([1 1],[0 2.12],'k--');
plot([0 1],[2.12 2.12],'k--');

xlabel(axes1,'A');
ylabel(axes1,'g_{eff}');

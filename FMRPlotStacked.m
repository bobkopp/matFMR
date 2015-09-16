function y=FMRPlotStacked(FMRdatasets,varargin)

scaleCurves = 1;
xaxisgfactors = 0;
withIntegratedFits=0;
smoothCurves = 0;
withFits=0;
withAnnotations = 1;
style='';
withCustomTitles = 0;
if nargin > 1
    for i=1:nargin-1
        if strcmp(varargin{i},'unscaled'), scaleCurves = 0; end
        if strcmp(varargin{i},'xaxisgfactors'),xaxisgfactors=1; end
        if strcmp(varargin{i},'integratedfits'),fits=varargin{i+1}; i=i+1; withIntegratedFits=1; end
        
        if strcmp(varargin{i},'fits'),fits=varargin{i+1}; i=i+1; withFits=1; end
        if strcmp(varargin{i},'smooth'), smoothCurves=1; end
        if strcmp(varargin{i},'style'),style=varargin{i+1}; i=i+1; end
        if strcmp(varargin{i},'noAnnotations'), withAnnotations=0; end
        if strcmp(varargin{i},'CustomTitles'), withCustomTitles=1; CustomTitles=varargin{i+1}; i=i+1; end
    end
end

offset  = 0;
top = 0;
bottom = 0;
maxfield = 0;

for i=1:length(FMRdatasets)
	FMRdata=FMRdatasets(i);
    if smoothCurves, FMRdata.datDeriv=FMRdata.datDerivSmooth; end
    if length(FMRdata.datDeriv)>0
        if scaleCurves==1, scaleFactor(i)=1/max(abs(FMRdata.datDeriv)); else, scaleFactor(i)=1; end
        offset = -bottom + .5 * max(abs(FMRdata.datDeriv)) * scaleFactor(i);
        offsets(i)=offset;
        normedDeriv=FMRdata.datDeriv *scaleFactor(i) - offset;
        top=max(top,max(normedDeriv))+.1;
        bottom=min(bottom,min(normedDeriv))-.3;
        tops(i)=max(normedDeriv)+.2;
        bottoms(i)=min(normedDeriv)-.3;
        if xaxisgfactors ==0, xvalues=FMRdata.fields*1000; else xvalues=FMRdata.ginverse; end
        
        plot (xvalues,normedDeriv,style)
        hold on;
        xlim([0 max(xvalues)]);
        ylim([bottom top]);
        
        if withIntegratedFits
            if fits(i).frequency > 0
                fields=FMRdata.fields;
                ginverse=FMRdata.ginverse;
                if xaxisgfactors ==0, xvalues2=1000*(.5*fields(1:end-1)+.5*fields(2:end)); else xvalues2=(.5*ginverse(1:end-1)+.5*ginverse(2:end)); end

                plot(xvalues2,(diff(fits(i).fit(fields))./diff(fields))*scaleFactor(i)-offsets(i),[style '--']);

                if fits(i).components>1
                    for j=1:fits(i).components
                        if strcmp(fits(i).fittypes{j},'cubic')
                            plot(1000*(.5*fields(1:end-1)+.5*fields(2:end)),(diff(fits(i).a(j)*FMRSpectrumCubic(fields,fits(i).g(j),fits(i).frequency,fits(i).Ban(j),fits(i).K2toK1(j),fits(i).lw(j)))./diff(fields))*scaleFactor(i)-offsets(i),[style ':']);
                        elseif strcmp(fits(i).fittypes{j}, 'uniaxial')
                            plot(1000*(.5*fields(1:end-1)+.5*fields(2:end)),(diff(fits(i).a(j)*FMRSpectrumUniaxial(fields,fits(i).g(j),fits(i).frequency,fits(i).Ban(j),fits(i).K2toK1(j),fits(i).lw(j)))./diff(fields))*scaleFactor(i)-offsets(i),[style ':']);
                        end
                    end
                end
            end
                
        end
        
        if withFits
            if fits(i).frequency > 0
                fields=FMRdata.fields;
                ginverse=FMRdata.ginverse;
                if xaxisgfactors ==0, xvalues2=1000*fields; else xvalues2=.5*ginverse; end

                plot(xvalues2,fits(i).fit(fields)*scaleFactor(i)-offsets(i),[style '--']);

                if fits(i).components>1
                    for j=1:fits(i).components
                        if strcmp(fits(i).fittypes{j},'cubic')
                            plot(xvalues2,fits(i).a(j)*FMRSpectrumDerivativeCubic(fields,fits(i).g(j),fits(i).frequency,fits(i).Ban(j),fits(i).K2toK1(j),fits(i).lw(j))*scaleFactor(i)-offsets(i),[style ':']);
                        elseif strcmp(fits(i).fittypes{j}, 'uniaxial')
                            plot(xvalues2,fits(i).a(j)*FMRSpectrumDerivativeUniaxial(fields,fits(i).g(j),fits(i).frequency,fits(i).Ban(j),fits(i).K2toK1(j),fits(i).lw(j))*scaleFactor(i)-offsets(i),[style ':']);
                        end
                    end
                end
            end
                
        end
        
        maxfield=max(maxfield,max(xvalues)); 
        if xaxisgfactors ==0,  xlabel('B (mT)'); else, xlabel('1/g'); end
        
        
    end
end

for i=1:length(FMRdatasets)
    midheight(i)=.5*(tops(i)+bottoms(i));
    totalheight(i)=tops(i)-bottoms(i);

    if withCustomTitles
        text(.02*maxfield,tops(i),CustomTitles{i},'FontWeight','bold','FontSize',12,'VerticalAlignment','top');
    else
        text(.02*maxfield,tops(i),FMRdatasets(i).sample,'FontWeight','bold','FontSize',12,'VerticalAlignment','top');
    end
       
    FMRQuant=FMRQuantify(FMRdatasets(i));

    if withAnnotations
        dataString(1)={['g_{eff} = ', sprintf('%0.2f', FMRQuant.geff) ' / ' 'A = ',sprintf('%0.2f',FMRQuant.A)]};
        dataString(2)={['\DeltaB_{FWHM} = ',sprintf('%0.0f',FMRQuant.dBFWHM*1000), ' mT' ' / ' '\alpha = ' sprintf('%0.2f',FMRQuant.factorAlpha)]};

    text(maxfield,tops(i),dataString,'HorizontalAlignment','right','VerticalAlignment','top');
    end
        
end
set(gca,  'YTick',offsets,  'YTickLabel',{''},'FontSize',12)
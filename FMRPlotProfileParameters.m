function y=FMRPlotProfileParameters(FMRdata,depth)

    for i=1:length(FMRdata)
        FMRq(i)=FMRQuantify(FMRdata(i));
    end
   
    subplot(1,5,1)
    plot([FMRq.totalAbs],depth,'x')
    ylabel('Depth')
    xlabel('Total Abs.')
    
    subplot(1,5,2)
    plot([FMRq.geff],depth,'x')
    xlabel('g_{eff}')
    
    subplot(1,5,3)
    plot([FMRq.dBFWHM],depth,'x')
    xlabel('\DeltaB_{FWHM}')
    
    subplot(1,5,4)
    plot([FMRq.A],depth,'x')
    xlabel('A')
    
    subplot(1,5,5)
    plot([FMRq.factorAlpha],depth,'x')
    xlabel('\alpha')
    
end
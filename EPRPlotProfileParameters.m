function y=EPRPlotProfileParameters(FMRdata,depth)

    for i=1:length(FMRdata)
        FMRq(i)=EPRQuantify(FMRdata(i));
    end
   
    subplot(1,6,1)
    plot([FMRq.totalAbs],depth,'x')
    ylabel('Depth')
    xlabel('Total Abs.')
    
    subplot(1,6,2)
    plot([FMRq.MnII],depth,'x')
    xlabel('Mn(II)')
    
    subplot(1,6,3)
    plot([FMRq.FeIII],depth,'x')
    xlabel('Fe(III)')
    
    subplot(1,6,4)
    plot([FMRq.FRC],depth,'x')
    xlabel('FRC (g=2)')

    subplot(1,6,5)
    plot([FMRq.MnII]./[FMRq.totalAbs],depth,'x')
    xlabel('Mn(II)/Total')

    subplot(1,6,6)
    plot([FMRq.FRC]./[FMRq.MnII],depth,'x')
    xlabel('FRC/Mn(II)')
end
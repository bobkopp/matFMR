function y=FMRPlotSpectrumAndCoercivity(FMRData,RmgData)

    for i=1:length(RmgData)
        if length(RmgData(i).samplename) > 0 
            subplot(length(FMRData),2,2*i-1)
            RmgSIRMDerivativePlot(RmgData(i),'AcqOnly')
        end
    end

    for i=1:length(FMRData)
        if length(FMRData(i).sample) > 0 
            subplot(length(FMRData),2,2*i)
            plot(FMRData(i).fields*1000,FMRData(i).datDerivSmooth,'b');
            title(FMRData(i).sample);
            xlabel('B (mT)');

        end
    end
    
    for i=1:length(FMRData)
        if length(FMRData(i).sample) > 0 
            subplot(length(FMRData),2,2*i)
            hold on; axis manual
            plot(FMRData(i).fields*1000,FMRData(i).datDeriv,'g');
        end
    end
end
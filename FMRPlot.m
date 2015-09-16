function y=FMRPlot(FMRdata)

FMRQuant=FMRQuantify(FMRdata)
plot (FMRdata.fields*1000,FMRdata.datInt/max(FMRdata.datInt),'g');
hold on;
plot (FMRdata.fields*1000,FMRdata.datDeriv/max(FMRdata.datDeriv),'b');
plot (FMRdata.fields*1000,FMRdata.datIntSmooth/max(FMRdata.datInt),'r');
plot (FMRdata.fields*1000,FMRdata.datDerivSmooth/max(FMRdata.datDeriv),'y');
title(FMRdata.sample);
xlim([0 max(FMRdata.fields)*1000]);
ylim([1.1*min(FMRdata.datDeriv)/max(FMRdata.datDeriv) 1.1]);
xlabel('B (mT)');

dataString(1)={['g_{eff} = ', num2str(FMRQuant.geff)]};
dataString(2)={['A = ',num2str(FMRQuant.A)]};
dataString(3)={['\DeltaB_{FWHM} = ',num2str(FMRQuant.dBFWHM*1000), ' mT']};
dataString(4)={['total = ', num2str(FMRQuant.totalAbs)]};
dataString(5)={['\alpha = ', num2str(FMRQuant.factorAlpha)]};
dataString(6)={['\beta = ', num2str(FMRQuant.factorBeta)]};
dataString(7)={['\alpha` = ', num2str(FMRQuant.factorAlphaPrime*1000) ' mT']};


statBox=annotation('textbox',[.2,.45,0,0]);
set(statBox,'LineStyle','none');
set(statBox,'FitHeightToText','on');
set(statBox,'string',dataString);
hold off;

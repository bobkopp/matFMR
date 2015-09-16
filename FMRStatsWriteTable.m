function y=FMRStatsTableWrite(FMRData,filename,varargin)

headers = ['Sample' '\tTotal Abs' '\tgeff' '\tA' '\tdBFWHM (mT)' '\talpha' '\tbeta'];


fid=fopen([filename '.asc'],'w');
fprintf(fid,[headers '\n']);

for i=1:length(FMRData)
    clear dataString;
    quant = FMRQuantify(FMRData(i));
    dataString=[FMRData(i).sample '\t'];
    dataString=[dataString num2str(quant.totalAbs) '\t'];
    dataString=[dataString num2str(quant.geff) '\t'];
    dataString=[dataString num2str(quant.A) '\t'];
    dataString=[dataString num2str(quant.dBFWHM*1000) '\t'];
    dataString=[dataString num2str(quant.factorAlpha) '\t'];
    dataString=[dataString num2str(quant.factorBeta) '\t'];
    dataString=[dataString '\n'];
    fprintf(fid,dataString);
end

fclose(fid);

end
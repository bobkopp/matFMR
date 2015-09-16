function y=EPRStatsTableWrite(FMRData,filename,varargin)

headers = ['Sample' '\tTotal Abs' '\tMnII' '\tMnII dB' '\tFeIII' '\tFRC (g=2)' '\tMnII/Total' '\tFRC/MnII'];


fid=fopen([filename '.asc'],'w');
fprintf(fid,[headers '\n']);

for i=1:length(FMRData)
    clear dataString;
    quant = EPRQuantify(FMRData(i));
    dataString=[FMRData(i).sample '\t'];
    dataString=[dataString num2str(quant.totalAbs) '\t'];
    dataString=[dataString num2str(quant.MnII) '\t'];
    dataString=[dataString num2str(quant.MnIIdB) '\t'];
    dataString=[dataString num2str(quant.FeIII) '\t'];
    dataString=[dataString num2str(quant.FRC) '\t'];
    dataString=[dataString num2str(quant.MnII/quant.totalAbs) '\t'];
    dataString=[dataString num2str(quant.FRC/quant.MnII) '\t'];
    dataString=[dataString '\n'];
    fprintf(fid,dataString);
end

fclose(fid);

end
function y=FMRWriteSpectrum(spectrum,filename)

fid=fopen([filename '.asc'],'w');
fprintf(fid,'Filename:\t');
fprintf(fid,[filename '.par']);
fprintf(fid,'\n\n');
fprintf(fid,'X [G]\tIntensity\n');
fclose(fid);

M=[spectrum.fields*10000 spectrum.datDeriv];
dlmwrite([filename '.asc'],M,'delimiter','\t','precision','%.6f','-append');

fid=fopen([filename '.par'],'w');
fprintf(fid,['MF  ',num2str(spectrum.frequency/1e9)]);
fprintf(fid,'\n');
fprintf(fid,['MP  ',num2str(spectrum.pwr)]);
fprintf(fid,'\n');
fprintf(fid,['JSD ',num2str(spectrum.sweeps)]);

fprintf(fid,'\n');
fclose(fid);
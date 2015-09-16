function y=FMRFitsMultiSingleWriteTable(datasets,fitsets,filename)

headers = ['Sample' '\tFit Type' '\tComponent' '\tComponent Type' '\tComponent Fraction' '\tComponent Area' '\ta' '\ta error' '\tg' '\tg error' ...
    '\tBan' '\tBan error' '\tK2/K1' '\tK2/K1 error' '\tlw' '\tlw error' '\tResidual''];


fid=fopen([filename '.asc'],'w');
fprintf(fid,[headers '\n']);

for i=1:length(datasets)
    dataset=datasets{i};
    fits=fitsets{i};
    for j=1:length(dataset)
        for k=1:fits(j).components
            fprintf(fid,[dataset(j).sample '\t']);
            fprintf(fid,[fits(j).fittype '\t']);
            fprintf(fid,[num2str(k) '\t']);
            fprintf(fid,[fits(j).fittypes{k} '\t']);
            fprintf(fid,[num2str(fits(j).f(k)) '\t']);
            fprintf(fid,[num2str(fits(j).compArea(k)) '\t']);
            fprintf(fid,[num2str(fits(j).a(k)) '\t']);
            fprintf(fid,[num2str(fits(j).aerror(k)) '\t']);
            fprintf(fid,[num2str(fits(j).g(k)) '\t']);
            fprintf(fid,[num2str(fits(j).gerror(k)) '\t']);
            fprintf(fid,[num2str(fits(j).Ban(k)) '\t']);
            fprintf(fid,[num2str(fits(j).Banerror(k)) '\t']);
            fprintf(fid,[num2str(fits(j).K2toK1(k)) '\t']);
            fprintf(fid,[num2str(fits(j).K2toK1error(k)) '\t']);
            
            fprintf(fid,[num2str(fits(j).lw(k)) '\t']);
            
            fprintf(fid,[num2str(fits(j).lwerror(k)) '\t']);
            
            fprintf(fid,[num2str(fits(j).resnorm) '\t']);
            
            fprintf(fid,'\n');
        end
    end
    
end

fclose(fid);

end
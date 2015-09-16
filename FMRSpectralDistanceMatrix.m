function y=FMRSpectralDistanceMatrix(spectrumArray)

for i=1:length(spectrumArray)-1
    for j=i+1:length(spectrumArray)
        distij=FMRSpectralDistance(spectrumArray(i),spectrumArray(j));
        y(i,j)=distij.NormalizedToAbs;
%        y(i,j)=distij.NormalizedToDer;
        y(j,i)=y(i,j);
    end
end

function y=FMRSpectralDistance(a,b)

aQuantify=FMRQuantify(a);
bQuantify=FMRQuantify(b);

aNormAbs=aQuantify.totalAbs;
bNormAbs=bQuantify.totalAbs;

aNormDer=norm(a.datDeriv);
bNormDer=norm(b.datDeriv);

aDerivNormAbs=a.datDeriv/aNormAbs;
aDerivNormDer=a.datDeriv/aNormDer;
bDerivNormAbs=b.datDeriv/bNormAbs;
bDerivNormDer=b.datDeriv/bNormDer;

y.NormalizedToAbs=norm(aDerivNormAbs-bDerivNormAbs);
y.NormalizedToDer=norm(aDerivNormDer-bDerivNormDer);
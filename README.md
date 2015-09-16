# matFMR 1.1

by Robert E. Kopp
March 28, 2007
Licensed under the terms of the GNU General Public License

Released as electronic supplementary material for:

	R. E. Kopp (2007), The Identification and Interpretation of Microbial Biogeomagnetism.
	Ph.D. Thesis, Division of Geological and Planetary Sciences,
California Institute of Technology.

-----

## Contents:
	1. Requirements
	2. Basic Function Definitions
	3. Advanced Function Definitions
	4. Example usage

*****

## REQUIREMENTS

matFMR was produced in MATLAB 7.1. It extensively uses the Curve Fitting toolbox.

## BASIC FUNCTION DEFINITIONS

FMRBlankSub(FMRdata, FMRblank)

	FMRBlankSub produces a FMR spectrum structure by subtracting FMRblank from
	FMRdata.
	
FMRImport(filetoload)

	FMRImport reads in a Bruker EMX-format .asc and .par file to produce a FMR
	spectrum structure.
	
	filetoload: .asc file to load
	
FMRPlot(FMRdata)

	Plots spectrum FMRdata in integrated and derivative forms, both smoothed
	and unsmoothed, and displays empirical parameters derived from spectrum.
	
FMRPlotStacked(FMRdatasets,[option1, [option2, ...]])

	Produces a stacked plot of the FMR derivative spectra in cell array FMR
	datasets, and displays empirical parameters.
	
	Possible options:
		'CustomTitles' - uses the spectrum titles specified by the cell array
			in the next parameter
		'noAnnotations' - does not display empirical parameters
		'smooth' - smooths the spectra
		'style' - uses the MATLAB line style specified by the next parameter
		'unscaled' - do not normalize FMR spectra
		'xaxisgfactors' - use g rather than B as the x axis

FMRQuantify(FMRdata)

	Generates a structure containing empirical parameters derived from the
	spectrum FMRdata.
	
FMRStatsTableWrite(FMRData,filename)

	Writes an ASCII table with file filename containing FMR parameters.

FMRSummaryPlotdBvsA(datasets,[labels])
FMRSummaryPlotgeffvsA(datasets,[labels])
FMRSummaryPlotTotalAbsvsA(datasets,[labels])

	Produces summary cross-plots of empirical parameters dB or geff or total
	absorption (integral of integrated spectrum) against A.
	
	datasets: cell array of arrays of FMR spectrum structures
	labels: optional array of strings containing a label for each array of FMR
		spectrum structures in datasets

FMRWriteSpectrum(spectrum,filename)

	Write spectrum to Bruker EMX-format .asc and .par files. Useful for saving
	blank-corrected spectra.

EPRQuantify(FMRdata)

	Generates a structure containing empirical EPR parameters derived from the
	spectrum FMRdata. Parameters include estimates of paramagnetic Fe(III) from the height of
	the g = 4.3 peak above a linearly interpolated baseline, Mn(II) from the
	difference between the maximum and minimum of the peak 25 mT above g = 2
	(which should be the sixth peak of the Mn(II) sextet), and "free radicals" estimated
	from the difference between the maximum and the minimum of the g = 2
	peak.
	
EPRStatsTableWrite(FMRData,filename)

	Writes an ASCII table with file filename containing EPR parameters.

3. ADVANCED FUNCTION DEFINITIONS

FMRFilter(FMRdata,criterion,[opts])

	Returns a FMR spectrum structure identical to FMRdata but containing
	only the points identified (by position number within array, not by field)
	in criterion. If opts is set to 'smooth', the data will also be smoothed.

FMRFitDerivativeMixMulti(FMRdata,fit1type,fit1param,[fit2type],[fit2param],...)

	FMRFitDerivativeMixMulti produces a FMR fit structure. It returns the best
	fit to the derivate spectrum or spectra in FMRdata. Whem multiple spectra
	are in FMRdata, it seeks to optimize the fit components so as to provide
	the best fit to all normalized spectra.
	
	As general advice, uniaxial models are much more efficient computationally.
	The results you get can be highly dependent upon the seed values picked.
	It is a good idea to compare the quality of fit and the confidence intervals
	on the fit parameters arrived at from different seeds.
	
	FMRdata: spectrum or spectra to be fit
	
	fit1type, [fit2type, ...]: possible options are
		'cubic' - spectrum generated using cubic anisotropy condition, K2/K1=0
		'cubic2nd' - spectrum generated using cubic anisotropy condition
		'cubic-fixedBan' - cubic anisotropy, K2/K1 = 0, Ban fixed
		'cubic-fixedg' - cubic anisotropy, K2/K1 = 0, g fixed
		'uniaxial' - spectrum generated using uniaxial anisotropy, K2/K1 = 0
		'uniaxial2nd' - spectrum generated using uniaxial anisotropy
		'uniaxial-fixedBan' - uniaxial anisotropy, K2/K1 = 0, Ban fixed
		'uniaxial-fixedg' - uniaxial anisotropy, K2/K1 = 0, g fixed
		'gaussian' - derivative Gaussian spectrum with peak location specified by g
					 and standard deviation specified by lw
		'fixedcubic' - cubic spectrum with known parameters
		'fixeduniaxial' - uniaxial spectrum with known parameters
		'fixedgaussian' - derivative Gaussian spectrum with known parameters
		
	fit1param, [fit2param, ...]: array specifying parameters for preceeding
								 fittype
								 
		for 'cubic' and 'uniaxial', takes form of  [g Ban lw]
			where g, Ban, and lw are seed values for g, Ban, and
			linewidth; e.g. [2.12 0 0.03]
			
		for 'cubic2nd' and 'uniaxial2nd', takes form of [g Ban K2toK1 lw]
			e.g. [1 2.12 0 0 0.03]
			
		for 'cubic-fixedBan', 'cubic-fixedg', 'uniaxial-fixedBan', and 'uniaxial-fixedg'
			parameters are as for 'cubic' and 'uniaxial'
			
		for 'fixeduniaxial' and 'fixedcubic', parameters are as for 'uniaxial2nd'
			and 'cubic2nd'
			
		for 'gaussian' and 'fixedgaussian', takes form of [g lw] where
			g and lw are seed or fixed values for peak and linewidth; e.g. [2.12 .03]
	
FMRLinearInterpolate(FMRdata,lowfield,highfield)

	Returns a FMR data structure identical to FMR data but with both the derivative
	and integrated spectra linearly interpolated between the field values specified
	by lowfield and highfield.
	

FMRRmgStatDepthProfiles(RmgData,Rmgdepths,FMRdata,FMRdepths,plots,[massfactor,[linespec,[opts]]])

	Requires that matRockmag be installed in the MATLAB path if rock magnetic
	parameters are called.
	
	Constructs profile subplots of parameters vs. depth.
	
	RmgData: Array of rock magnetic data (produced by matRockmag through RmgImport). Can be
			 left as empty array [].
	Rmgdepths: Depth of each member of RmgData.
	FMRdata: Array of FMR data. Can be left as empty array [].
	FMRdepths: Depth of each member of FMRdata.
	massfactor: Mass factor by which to multiply rock magnetic parameter to produce kg
				normalized parameters. Defaults to 10^6 (i.e., assumes 'vol' field of
				rock magnetic data is mass in mg).
	linespec: MATLAB-style line spec
	opts: 'smooth' will cause the curve to be smoothed.
	plots: Cell array of strings identifying each parameter to be profiled. Any 
		   string not recognized will produced a blank subplot. Options are:
		   		'sirm/vol' - ratio of SIRM to volume/mass
		   		'logsirm/vol' - ratio of SIRM to volume/mass, plotted on log scale
		   		'armsuscep' - ARM susceptibility
		   		'arm100/irm' - ratio of ARM at 100 uT to IRM
		   		'arm500/irm' - ratio of ARM at 500 uT to IRM
		   		'dfirm/db' - value of the derivative of IRM at the max field value
		   		'armsuscep/irm' - ratio of ARM susceptibility to IRM
		   		'cisowskir' - Cisowski R parameter
		   		'hcr' - Hcr, determined from IRM acquisition and AF
		   		'irm30/irm100' - ratio of IRM at 30 mT to IRM at 100 mT
		   		'irm100/irm300' - ratio of IRM at 100 mT to IRM at 300 mT
		   		'maf-irm' - median acquisition field of IRM
		   		'mdf-irm' - median destructive field of IRM
		   		'mdf-arm' - median destructive field of ARM
		   		'mdf-irmtoarm' - ratio of MDF of IRM to MDF of ARM
		   		'mdf-arm-mdf-irm' - difference between MDF of ARM and MDF of IRM
		   		'dp-irmacq' - dispersion parameter of IRM acquisition spectrum
		   		
		   		'fmrabs' - total FMR absorption 
				'geff' - effective g value
				'dbfwhm' - full-width half maximum
				'db95' - width around median integrated FMR absorption containing 95% of absorption
				'a' - asymmetry ratio A
				'alpha' - factor alpha
				'mnii' - estimated Mn(II)
				'feiii' - estimated paramagnetic Fe(III)
				'frc' - estimated free radicals
				'mnii/total' - ratio of estimated Mn(II) to total absorption
				'feiii/total' - ratio of estimated paramagnetic Fe(III) to total absorption
				'frc/total' - ratio of estimated free radicals to total absorption
				'frc/mnii' - ratio of estimated free radicals to estimated Mn(II)
				'feiii/mnii' - ratio of estimated Fe(III) to Mn(II)
				'sirm/fmrabs' - ratio of SIRM to total absorption
				'logsirm/fmrabs' - ratio of SIRM to total absorption, log scale
				'armsuscep/fmrabs' - ratio of ARM susceptibility to total absorption

FMRSpectralDistance(a,b)

	Determines the norm of the difference between spectra a and b, normalized
	either to total absorption or to the norm of the derivative spectra.
	
FMRSpectralDistanceMatrix(spectrumArray)

	From the one-dimensional array of FMR data structures in spectrumArray,
	constructs a distance matrix between absorption-normalized spectra.

FMRSpectrumCubic(Bvalues,gtrue,frequency,Ban,K2toK1,linewidth)
FMRSpectrumUniaxial(Bvalues,gtrue,frequency,Ban,K2toK1,linewidth)
FMRSpectrumDerivativeCubic(Bvalues,gtrue,frequency,Ban,K2toK1,linewidth)
FMRSpectrumDerivativeUniaxial(Bvalues,gtrue,frequency,Ban,K2toK1,linewidth)

	Generates an integrated or derivative spectrum, evaluated at the points
	in array Bvalues (in T), using the specified parameters.
	
	Bvalues: field values (T)
	gtrue: true g value
	frequency: measuring frequency (Hz)
	Ban: anisotropy field (T)
	K2toK1: K2/K1
	linewidth: Gaussian linewidth (T)

FMRSyntheticQuantifyInt(fields,datInt,frequency)
FMRSyntheticQuantify(fields,datDeriv,frequency)

	Generates a merged FMR spectrum structrure and FMRQuantify structure for
	a synthetic spectrum.
	
	fields: fields at which synthetic spectrum was generated (T)
	datInt: value of integrated spectrum at points in field
	datDeriv: value of integrated spectrum at points in field
	frequency: frequency used in generating spectrum (Hz)


## EXAMPLE USAGE

Load FMR data from files, fit uniaxial spectra to the data, and plot up
the data.

	filesV = dir('V*.asc')
	for i=1:length(filesV)
		FMRdataV(i)=FMRImport(filesV(i).name);
		FMRfitV(i)=FMRFitDerivativeMixMulti(FMRdataV,'uniaxial2nd',[2 0 0 .03])
	end
	
	filesA = dir(A*.asc')
	for i=1:length(filesA)
		FMRdataA(i)=FMRImport(filesA(i).name);
		FMRfitA(i)=FMRFitDerivativeMixMulti(FMRdataA,'uniaxial2nd',[2 0 0 .03])
	end
	
	h=figure;
	for i=1:length(FMRdataV)
		subplot(length(FMRdataV),1,i);
		plot(FMRdataV(i).fields,FMRdataV(i).datDeriv,'k-');
		hold on;
		plot(FMRdataV(i).fields,FMRfitV(i).fitValues,'r--');
	end
	
	h=figure;
	for i=1:length(FMRdataA)
		subplot(length(FMRdataA),1,i);
		plot(FMRdataA(i).fields,FMRdataA(i).datDeriv,'k-');
		hold on;
		plot(FMRdataA(i).fields,FMRfitA(i).fitValues,'r--');
	end
	
	h=figure;
	FMRSummaryPlotdBvsA({FMRdataV,FMRdataA},{'V','A'});
	
	h=figure;
	FMRSummaryPlotgeffvsA({FMRdataV,FMRdataA},{'V','A'});
	
	h=figure;
	FMRRmgStatDepthProfiles([],[],FMRdataA,[.1 .2 .3],{'fmrabs','geff','dbfwhm','a','alpha'});

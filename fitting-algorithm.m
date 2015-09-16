for i=1:length(FMR)

	FMRi(i)=FMRLinearInterpolate(FMR(i),0.15,0.185)
	
	FMRFit1(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[2.076 .0906 .0373]) 
	FMRFit1A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[2.076 -.0906 .0373]) 
	for j=1:3
		FMRFit1(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit1(i).g FMRFit1(i).Ban FMRFit1(i).lw],'startingAmplitudes',FMRFit1(i).a)
		FMRFit1A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit1A(i).g FMRFit1A(i).Ban FMRFit1A(i).lw],'startingAmplitudes',FMRFit1A(i).a)
	end
	
	FMRFit2(i)=FMRFitDerivativeMixMulti(FMRi(i),'fixeduniaxial',[FMRFit1(i).g FMRFit1(i).Ban 0 FMRFit1(i).lw],'uniaxial',[2.12 .1 .1],'startingAmplitudes',[FMRFit1(i).a(1) .3*FMRFit1(i).a(1)])
	
	FMRFit2A(i)=FMRFitDerivativeMixMulti(FMRi(i),'fixeduniaxial',[FMRFit1A(i).g FMRFit1A(i).Ban 0 FMRFit1A(i).lw],'uniaxial',[2.12 -.1 .1],'startingAmplitudes',[FMRFit1A(i).a(1) .3*FMRFit1A(i).a(1)])
	
	
	
	for j=1:3
		FMRFit2(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit2(i).g(1) FMRFit2(i).Ban(1) FMRFit2(i).lw(1)],'uniaxial',[FMRFit2(i).g(2) FMRFit2(i).Ban(2) FMRFit2(i).lw(2)],'startingAmplitudes',FMRFit2(i).a)
		
		FMRFit2A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit2A(i).g(1) FMRFit2A(i).Ban(1) FMRFit2A(i).lw(1)],'uniaxial',[FMRFit2A(i).g(2) FMRFit2A(i).Ban(2) FMRFit2A(i).lw(2)],'startingAmplitudes',FMRFit2A(i).a)
		
	end
	
	FMRFit3(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial2nd',[FMRFit2(i).g(1) FMRFit2(i).Ban(1) FMRFit2(i).lw(1) 0],'uniaxial2nd',[FMRFit2(i).g(2) FMRFit2(i).Ban(2) FMRFit2(i).lw(2) 0],'startingAmplitudes',FMRFit2(i).a)
	FMRFit3A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial2nd',[FMRFit2A(i).g(1) FMRFit2A(i).Ban(1) FMRFit2A(i).lw(1) 0],'uniaxial2nd',[FMRFit2A(i).g(2) FMRFit2A(i).Ban(2) FMRFit2A(i).lw(2) 0],'startingAmplitudes',FMRFit2A(i).a)
	
	for j=1:3
		FMRFit3(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial2nd',[FMRFit3(i).g(1) FMRFit3(i).Ban(1) FMRFit3(i).lw(1) FMRFit3(i).K2toK1(1)],'uniaxial2nd',[FMRFit3(i).g(2) FMRFit3(i).Ban(2) FMRFit3(i).lw(2) FMRFit3(i).K2toK1(2)],'startingAmplitudes',FMRFit3(i).a)
		FMRFit3A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial2nd',[FMRFit3A(i).g(1) FMRFit3A(i).Ban(1) FMRFit3A(i).lw(1) FMRFit3A(i).K2toK1(1)],'uniaxial2nd',[FMRFit3A(i).g(2) FMRFit3A(i).Ban(2) FMRFit3A(i).lw(2) FMRFit3A(i).K2toK1(2)],'startingAmplitudes',FMRFit3A(i).a)
	end
	
	FMRFit4(i)=FMRFitDerivativeMixMulti(FMRi(i),'fixeduniaxial',[FMRFit2(i).g(1) FMRFit2(i).Ban(1) 0 FMRFit2(i).lw(1)],'fixeduniaxial',[FMRFit2(i).g(2) FMRFit2(i).Ban(2) 0 FMRFit2(i).lw(2)],'uniaxial',[2.12 .2 .2],'startingAmplitudes',[FMRFit2(i).a ; FMRFit2(i).a(1)/3])
	FMRFit4A(i)=FMRFitDerivativeMixMulti(FMRi(i),'fixeduniaxial',[FMRFit2(i).g(1) FMRFit2(i).Ban(1) 0 FMRFit2(i).lw(1)],'fixeduniaxial',[FMRFit2(i).g(2) FMRFit2(i).Ban(2) 0 FMRFit2(i).lw(2)],'uniaxial',[2.12 -.2 .2],'startingAmplitudes',[FMRFit2(i).a ; FMRFit2(i).a(1)/3])
	
	
	for j=1:3
		if max(FMRFit4(i).g)<3.5
			FMRFit4(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit4(i).g(1) FMRFit4(i).Ban(1) FMRFit4(i).lw(1)],'uniaxial',[FMRFit4(i).g(2) FMRFit4(i).Ban(2) FMRFit4(i).lw(2)],'uniaxial',[FMRFit4(i).g(3) FMRFit4(i).Ban(3) FMRFit4(i).lw(3)],'startingAmplitudes',[FMRFit4(i).a])
		end
	
		if max(FMRFit4A(i).g)<3.5
			FMRFit4A(i)=FMRFitDerivativeMixMulti(FMRi(i),'uniaxial',[FMRFit4A(i).g(1) FMRFit4A(i).Ban(1) FMRFit4A(i).lw(1)],'uniaxial',[FMRFit4A(i).g(2) FMRFit4A(i).Ban(2) FMRFit4A(i).lw(2)],'uniaxial',[FMRFit4A(i).g(3) FMRFit4A(i).Ban(3) FMRFit4A(i).lw(3)],'startingAmplitudes',[FMRFit4A(i).a])
		end
	
	end
	FMRFitsMultiWriteTable({repmat(FMR(1:i),8,1)},{[FMRFit1(1:i) FMRFit1A(1:i) FMRFit2(1:i) FMRFit2A(1:i) FMRFit3(1:i) FMRFit3A(1:i) FMRFit4(1:i) FMRFit4A(1:i)]},'table-fits')
	save fits
	
end
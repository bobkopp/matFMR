function y=FMRRmgStatDepthProfiles(RmgData,Rmgdepths,FMRdata,FMRdepths,plots,varargin)

	linespec={'x'};

	[Rmgdepths,i]=sort(Rmgdepths);
	RmgData=RmgData(i);

	[FMRdepths,i]=sort(FMRdepths);
	FMRdata=FMRdata(i);

	mindepth=min([Rmgdepths FMRdepths]);
	maxdepth=max([Rmgdepths FMRdepths]);
	
	smoothdata=0;

	if length(varargin)>0
		massfactor=varargin{1};
	else
		massfactor = 10^6;
		% assume massin mg;
	end
	
	if length(varargin)>1
		linespec=varargin{2};
	end

	if length(varargin)>2
		for i=3:length(varargin)
			if strcmp(varargin{i},'smooth')
				smoothdata=1;
			end
		end
	end

    for i=1:length(RmgData)
        RmgStat(i)=RmgStats(RmgData(i));
    end

    for i=1:length(FMRdata)
        FMRq(i)=FMRQuantify(FMRdata(i));
        EPRq(i)=EPRQuantify(FMRdata(i));
    end
    
    [overlapDepths,overlapDepthsRmgI,overlapDepthsFMRI]=intersect(Rmgdepths,FMRdepths);
    
    n=length(plots);
    
    for i=1:n
    
    	subplot(1,n,i);
    	hold on;
    
    	currentplot=lower(plots{i});
    
    	test = strcmp(currentplot,'sirm/vol');
    	if test
    		data=[RmgStat.sIRMperkg]*massfactor;
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel('A m^2 / kg');
			title(['IRM @ ' num2str(RmgStat(1).dfIRMdBatField) ' ' RmgStat(1).units.Hcr]);
			ylabel('depth');
    	end
 
     	test = strcmp(currentplot,'logsirm/vol');
    	if test
    		data=[RmgStat.sIRMperkg]*massfactor;
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel('A m^2 / kg');
			title(['IRM @ ' num2str(RmgStat(1).dfIRMdBatField) ' ' RmgStat(1).units.Hcr]);
			ylabel('depth');
			ylim([mindepth maxdepth]);
			set(gca,'XScale','log');
			
			xlim([min(data) max(data)]);
    	end
 
    	test = strcmp(currentplot,'armsucep');
    	if test
    		data=[RmgStat.ARMsusceptibility];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.ARMsusceptibility);
			title('\chi_{ARM}');		
    	end
 
     	test = strcmp(currentplot,'arm100/irm');
    	if test
    		data=[RmgStat.ARMtoIRMat100uT];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			title('ARM / IRM at 100 uT');		
		end

     	test = strcmp(currentplot,'arm500/irm');
    	if test
    		data=[RmgStat.ARMtoIRMat500uT];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			title('ARM / IRM at 500 uT');		
		end

      	test = strcmp(currentplot,'dfirm/db');
    	if test
    		data=[RmgStat.dfIRMdB];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			title(['df_{IRM}/dB at ' num2str(RmgStat(1).dfIRMdBatField) ' mT']);		
		end
 
 		test = strcmp(currentplot,'armsuscep/irm');
    	if test
    		data=[RmgStat.ARMsusceptibilityToIRM];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.ARMsusceptibilityToIRM);
			title('\chi_{ARM} / IRM');		
    	end
    	
    	test = strcmp(currentplot,'cisowskir');
    	if test
    		data=[RmgStat.CisowskiR];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel('R');
			title('Cisowski R');
			
    	end

    	test = strcmp(currentplot,'hcr');
    	if test
    		data=[RmgStat.Hcr];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.Hcr);
			title('Bcr');			
		end
    	
    	test = strcmp(currentplot,'irm30/irm100');
    	if test
    		data=[RmgStat.IRM30toIRM100];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			title('IRM_{30}/IRM_{100}');			
		end
		
		test = strcmp(currentplot,'irm100/irm300');
    	if test
    		data=[RmgStat.IRM100toIRM300];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			title('IRM_{100}/IRM_{300}');			
		end
    	

    	test = strcmp(currentplot,'maf-irm');
    	if test
    		data=[RmgStat.MAFofIRM];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.MDF);
			title('MAF of IRM');			
    	end
    	
    	test = strcmp(currentplot,'mdf-irm');
    	if test
    		data=[RmgStat.MDFofIRM];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.MDF);
			title('MDF of IRM');			
		end

    	test = strcmp(currentplot,'mdf-arm');
    	if test
    		data=[RmgStat.MDFofARM];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.MDF);
			title('MDF of IRM');			
    	end

    	test = strcmp(currentplot,'mdf-irmatarm');
    	if test
    		data=[RmgStat.MDFofIRMatARM];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.MDF);
			title('MDF of IRM at ARM Field');			
		end

    	test = strcmp(currentplot,'mdf-arm-mdf-irm');
    	if test
    		data=[RmgStat.MDFofARM]-[RmgStat.MDFofIRMatARM];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel(RmgStat(1).units.MDF);
			title('MDF_{ARM} - MDF_{IRM @ ARM field}');			
		end


    	test = strcmp(currentplot,'dp-irmacq');
    	if test
    		data=[RmgStat.DPofIRMAcq];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,Rmgdepths,linespec{:});
			xlabel('DP');
			title('DP of IRM Acq.');			
		end


    	
    	test = strcmp(currentplot,'fmrabs');
    	if test
    		data=[FMRq.totalAbs];
     		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,FMRdepths,linespec{:})
    		title('Total ESR Abs.')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

    	
    	test = strcmp(currentplot,'geff');
    	if test
    		data=[FMRq.geff];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,FMRdepths,linespec{:})
    		title('g_{eff}')
    		xlabel('g_{eff}')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

		test = strcmp(currentplot,'dbfwhm');
    	if test
    		data=[FMRq.dBFWHM]*1000;
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
			title('\DeltaB_{FWHM}')
			xlabel('\DeltaB_{FWHM} (mT)')
			
		
 
    	end

		test = strcmp(currentplot,'db95');
    	if test
    		data=[FMRq.dB95]*1000;
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
			title('\DeltaB_{95}')
			xlabel('\DeltaB_{95} (mT)')
			
		
 
    	end

		test = strcmp(currentplot,'a');
    	if test
    		data=[FMRq.A];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,FMRdepths,linespec{:})
			title('A')
			xlabel('A')
    	end

		test = strcmp(currentplot,'alpha');
    	if test
    		data=[FMRq.factorAlpha];
    		if (smoothdata), data=smooth(data,5); end;
    	    plot(data,FMRdepths,linespec{:})
			title('\alpha')
    	end

		test = strcmp(currentplot,'mnii');
    	if test
    		data=[EPRq.MnII];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
    		title('Mn(II)')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

		test = strcmp(currentplot,'feiii');
    	if test
    		data=[EPRq.FeIII];
     		if (smoothdata), data=smooth(data,5); end;
     		plot(data,FMRdepths,linespec{:})
    		title('Fe(III)')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

		test = strcmp(currentplot,'frc');
    	if test
    		data=[EPRq.FRC];
     		if (smoothdata), data=smooth(data,5); end;
     		plot(data,FMRdepths,linespec{:})
    		title('FRC (g = 2.00)')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
     	end

		test = strcmp(currentplot,'mnii/total');
    	if test
    		data=[EPRq.MnII]./[EPRq.totalAbs];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
    		title('Mn(II)/Total')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

		test = strcmp(currentplot,'feiii/total');
    	if test
    		data=[EPRq.FeIII]./[EPRq.totalAbs];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
    		title('Fe(III)/Total')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

	    test = strcmp(currentplot,'frc/total');
    	if test
    		data=[EPRq.FRC]./[EPRq.totalAbs];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
    		title('FRC/Total')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

	    test = strcmp(currentplot,'frc/mnii');
    	if test
    		data=[EPRq.FRC]./[EPRq.MnII];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths,linespec{:})
    		title('FRC/Mn(II)')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
        end
    	
        test = strcmp(currentplot,'feiii/mnii');
    	if test
    		data=[EPRq.FeIII]./[EPRq.MnII];
     		if (smoothdata), data=smooth(data,5); end;
     		plot(data,FMRdepths,linespec{:})
    		title('Fe(III)/Mn(II)')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end
        
    	test = strcmp(currentplot,'sirm/fmrabs');
    	if test
    		data=[RmgStat(overlapDepthsRmgI).sIRMperkg]./[FMRq(overlapDepthsFMRI).totalAbs];
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths(overlapDepthsFMRI),linespec{:})
    		title('SIRM / FMR Abs.')
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end

    	test = strcmp(currentplot,'logsirm/fmrabs');
    	if test
    		data=[RmgStat(overlapDepthsRmgI).sIRMperkg]./[FMRq(overlapDepthsFMRI).totalAbs];
     		if (smoothdata), data=smooth(data,5); end;
     		plot(data,FMRdepths(overlapDepthsFMRI),linespec{:})
    		title('SIRM / FMR Abs.')
    		ylim([mindepth maxdepth]);
    		
			set(gca,'XScale','log');
			xlim([min(data) max(data)]);
			
    	end


    	test = strcmp(currentplot,'armsuscep/fmrabs');
    	if test
    		data=[RmgStat(overlapDepthsRmgI).ARMsusceptibility]./[FMRq(overlapDepthsFMRI).totalAbs]./[RmgData(overlapDepthsRmgI).volume]*massfactor;
    		if (smoothdata), data=smooth(data,5); end;
    		plot(data,FMRdepths(overlapDepthsFMRI),linespec{:})
    		title('\chi_{ARM} / FMR Abs.')
    
    		ylim([mindepth maxdepth]);
			xlim([min(data) max(data)]);
    	end
    
    end
    

	for i=1:n
		subplot(1,n,i);
		children=get(gca,'Children');
		xdata=get(children,'XData');
		ydata=get(children,'YData');
		if length(children)>1
			for j=1:length(children)
				minxs(j)=min(xdata{j});
				maxxs(j)=max(xdata{j});
				minys(j)=min(ydata{j});
				maxys(j)=max(ydata{j});
			end
		else
				minxs=min(xdata);
				maxxs=max(xdata);
				minys=min(ydata);
				maxys=max(ydata);
		end			
		mindepth=min([mindepth minys ]);
		maxdepth=max([maxdepth maxys]);
		subplot(1,n,i); ylim([mindepth maxdepth]);
		subplot(1,n,i); xlim([min(minxs) max(maxxs)]);
	end
    
end


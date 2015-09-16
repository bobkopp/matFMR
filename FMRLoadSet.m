function [sh,depth]=FMRLoadSet(code,varargin)
	writefiles = 1;
	sep='_';
	codesep=' ';
	if length(varargin)>0
		for i=1:length(varargin)
			if strcmpi(varargin{i},'nowrite')
				writefiles = 0;
			elseif strcmp(varargin{i},'separator')
				sep=varargin{i+1};
				i=i+1;
			elseif strcmp(varargin{i},'codesep')
				codesep=varargin{i+1};
				i=i+1;
			end
		end
	end
	code=[code codesep];
	files=dir(['*' code '*.asc'])
	for i=1:length(files)
		sh(i)=FMRImport(files(i).name);
		ind = strfind(sh(1).sample,code);
		d = regexp(sh(i).sample, [code '\s*(\d+)' sep '(\d+)'],'tokens');
		if length(d)>0
			depth(i) = str2num(d{1}{1}) + str2num(d{1}{2})*10^-length(d{1}{2})
		else
			d = regexp(sh(i).sample, [code '\s*(\d+)'],'tokens');
			depth(i) = str2num(d{1}{1});
		end
		if writefiles
			clf; FMRPlot(sh(i));
			saveas(gcf,[sh(i).sample '.eps'],'epsc2');
		end
	end
	[depth,sorti] = sort(depth);
	sh=sh(sorti);
	
	if writefiles
		FMRStatsWriteTable(sh,['FMRStats-' code(1:end-1)]);
	end
	FMRRmgStatDepthProfiles([],[],sh,-depth,{'fmrabs','alpha','geff','dbfwhm','a','mnii','feiii'})
end
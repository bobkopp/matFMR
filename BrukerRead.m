function varargout = BrukerRead(varargin)

% BRUKERREAD Open Bruker BE3ST files (.DTA / .DSC , .spc / .par) into the
% MATLAB workspace.
%
% BRUKERREAD ()
% BRUKERREAD ('/path/to/file')
% BRUKERREAD ('plot')
% BRUKERREAD ('/path/to/file','plot')
% [x, y] = BRUKERREAD (...)
% [x, y, info] = BRUKERREAD (...)
%
% BRUKERREAD when run without any inputs, opens a GUI so that the user can
% open the file themselves. BRUKERREAD can also accept a path to a file as
% an input if the path is put in 'quotes' and the extension (.DTA) is left
% off.
%
% BRUKERREAD now also works with Bruker .spc/.par from EMX machines -
% thanks to Muege Aksoyoglu (Uni. Freiburg, DE) who kindly donated some
% .spc/.par files for testing
%
% BRUKERREAD can be run with the optional 'plot' input, to plot the file
% being loaded
%
% BRUKERREAD outputs a x matrix (magnetic field / time), a y matrix
% (intensity) and an optional info field.
%
% If no outputs are selected then the x and y values are plotted
% With the graph option (graph is assigned value 1) the data is also
% plotted
%
% BRUKERREAD is fully functional with cw experiments and 3 dimensional cw
% experiments where there is an additional .YGF file. BRUKERREAD provides
% basic functionality for pulsed experiments, it has been tested with field
% swept echoes and fourier induced decays. More complex experiments such as
% ENDOR/ESEEM are untested and may have varying results but the
% information array should be complete.
%
% Inputs:
%    input1     - a string input to the path of a file
%    input2     - 'plot' draws a plot of the imported file
%
% Outputs:
%    output0    - plot
%                   A new figure showing the data
%    output1    - x axis
%                   Magnetic field
%    output2    - y axis
%                   Intensity
%    output3    - information
%                   Array of information about the loaded file
%
% Example: 
%    [x,y] = BrukerRead
%               GUI load a file
%
%    [x,y,info] = BrukerRead('/path/to/file.DTA','plot')
%               load x,y and info of file.DTA with to the workspace and 
%               plot x,y as a new figure
%
% Other m-files required:   none
%
% Subfunctions:             none
%
% MAT-files required:       none
%
%
% See also: EPRTOOLBOX CWPLOTTER


%                                        _                             _   
%                                       | |                           | |  
%  _ __ ___   ___  _ __ __ _  __ _ _ __ | |__  _   _  ___   _ __   ___| |_ 
% | '_ ` _ \ / _ \| '__/ _` |/ _` | '_ \| '_ \| | | |/ _ \ | '_ \ / _ \ __|
% | | | | | | (_) | | | (_| | (_| | | | | |_) | |_| |  __/_| | | |  __/ |_ 
% |_| |_| |_|\___/|_|  \__, |\__,_|_| |_|_.__/ \__, |\___(_)_| |_|\___|\__|
%                       __/ |                   __/ |                      
%                      |___/                   |___/                       
%
%
% M. Bye v12.7
%
% Author:       Morgan Bye
% Work address: Henry Wellcome Unit for Biological EPR
%               University of East Anglia
%               NORWICH, UK
% Email:        morgan.bye@uea.ac.uk
% Website:      http://www.morganbye.net/eprtoolbox/brukerread
% Jul 2012;     Last revision: 9-July-2012
%
% Approximate coding time of file:
%               10 hours
%
%
% Version history:
% Jun 12        > Muege Aksoyoglu (Uni. Freiburg, DE) kindly donated some
%                   .spc/.par files from an EMX machine allowing for
%                   support for these files
%               > GUI load window updated to "Bruker files" instead of
%                   seperate options for DTA and spc
%               > if statements for DTA/spc changed to switch statements
%                   for cleaner code
%
% Apr 12        > Better handling of pulse experiments
%               > Inputs changed so that 'plot' can be used by itself
%               > Figure plotting reconfigured and improved
%                       - imaginary channel for pulse experiments
%                       - legend is displayed for real/imaginary channels
%                       - x-axis corrected for time or mag field depending
%                       on FSE or FID experiment with correct units/label
%
% Feb 12        Minor edit - fix for occasional file load by command line
%               error
%
% Dec 11        Plotted figures have file names included in header
%
% Oct 11        > "Return" catch statement if user cancels a file load
%                   using the GUI
%               > Fixed bug with 'plot' input parameter
%
% Sept 11       Changed fid/fopen name reference system so that files are
%               still loaded even if not in the same location as MATLAB's
%               "Current Folder". Only a problem occasionally when called
%               from another script.

%                         Input arguments
% ========================================================================

switch nargin
    case 0
        % GUI get file
        [file , directory] = uigetfile({'*.DTA;*.spc','Bruker File (*.DTA,*.spc)'; ...
            '*.*',  'All Files (*.*)'},...
            'Load Bruker file');
        
        % if user cancels command nothing happens
        if isequal(file,0) %|| isequal(directory,0)
            return
        end
                
        % File name/path manipulation
        address = [directory,file];
        [~, name,extension] = fileparts(address);
        
        graph = '';
        
        
    case 1
        % If input is 'plot'
        if strcmp(varargin{1},'plot')
            graph = 'plot';
            
            % GUI get file
            [file , directory] = uigetfile({'*.DTA','Bruker BES3T File (*.DTA)'; ...
            '*.spc','Bruker Spc File (*.spc)'; ...
            '*.*',  'All Files (*.*)'},'Load Bruker file');
        
            % if user cancels command nothing happens
            if isequal(file,0) %|| isequal(directory,0)
                return
            end
                
            % File name/path manipulation
            address = [directory,file];
            [~, name,extension] = fileparts(address);
            
        % Or if input is a path
        else
            address = varargin{1};
            [directory,name,extension] = fileparts(address);
            graph = '';
        end
        
    case 2
        address = varargin{1};
        [directory,name,extension] = fileparts(address);
        graph = varargin{2};
end

%                         Parameter files
% ========================================================================

% Load .dsc/.par file
% ===================
switch extension
    case '.spc'
    filepar = [ directory '/' name '.par'];
    fid = fopen( filepar , 'r');
    
    delimiter = 10;
    
    case '.DTA'
    filedsc = [directory '/' name '.DSC'];
    fid = fopen( filedsc , 'r');
    
    delimiter = 10;
end

% Get characters from file
[string , string_characters] = fscanf( fid , '%c' );

% Close file, free up memory
fclose(fid);

% Convert character array into useful string array (insert line breaks)
parameters = [];

while isempty(string) == 0;
    [token, string] = strtok(string,char(delimiter));
    parameters = str2mat(parameters,token);
end

% Do basic reading of info file to decide on how to proceed with number
% crunching.

switch extension
    case '.spc'
        info.EXPT = 'CW';
        
    case '.DTA'
        
        info.EXPT  = strtrim(regexprep(parameters((strmatch('EXPT',parameters)),:),'EXPT',''));
        info.IRNAM = strtrim(regexprep(parameters((strmatch('IRNAM',parameters)),:),'IRNAM',''));
        info.XNAM  = strtrim(regexprep(parameters((strmatch('XNAM',parameters)),:),'XNAM',''));
        info.XUNI  = strtrim(regexprep(parameters((strmatch('XUNI',parameters)),:),'XUNI',''));
        info.TITL  = strtrim(regexprep(parameters((strmatch('TITL',parameters)),:),'TITL',''));
end


%                            Data files
% ========================================================================


% Load .dta/.spc file
switch extension
    case '.spc'
        fid   = fopen( [directory '/' name '.spc'], 'r');
        [y,n] = fread(fid,inf,'float');
        
    case '.DTA'
        fid   = fopen( [directory '/' name '.DTA'], 'r','ieee-be.l64');
        [y,n] = fread(fid,inf,'float64');
end

fclose(fid);


% Pulsed experiment splitting into real and imaginary channels
if strcmp(info.EXPT,'PLS')
    z = reshape(y,2,[]);
    y.real = z(1,:)';
    y.imag = z(2,:)';
end
    



%                            YGF files
% ========================================================================


% search directory for .YGF file
if exist([directory name '.YGF'],'file')
    
    % if exist, load .YGF , convert to usable matrix
    fid = fopen( [directory name '.YGF'], 'r','ieee-be.l64');
    [info.z_axis,info.z_axis_points] = fread(fid,inf,'float64');
    
    % find number of data points
    MagField.points = str2num(regexprep(parameters((strmatch('XPTS',parameters)),:),'XPTS',''));
    % reshape the y-axis into columns using number of data points
    y = reshape(y, MagField.points , []);
    
end




%                         Magnetic field
% ========================================================================


switch extension
    case '.spc'
        
        % Mag field centre point
        MagField.centre = str2num(regexprep(parameters((strmatch('HCF',parameters)),:),'HCF',''));
        
        % Mag field sweep width
        MagField.width  = str2num(regexprep(parameters((strmatch('HSW',parameters)),:),'HSW',''));
        
        % Mag field data points
        MagField.points = str2num(regexprep(parameters((strmatch('ANZ',parameters)),:),'ANZ',''));
        
        MagField.min      = MagField.centre - (MagField.width / 2);
        MagField.sampling = MagField.width / MagField.points;
        
        for k = 1:MagField.points
            x(k) = MagField.min + (MagField.sampling * ((k)-1));
        end
        
        x = x';
        
    case '.DTA'    

        % =================================
        % a = strmatch('XMIN',basic_info)  % find line of XMIN
        % mag_start = basic_info(a,:);     % read that line
        % regexprep(mag_start,'XMIN','')   % delete XMIN from that line
        % str2num                          % convert string to number
        %
        % in a shell this could easily be done with "awk"
        % =================================
        
        % Mag field start
        MagField.min = str2num(regexprep(parameters((strmatch('XMIN',parameters)),:),'XMIN',''));
        
        % Mag field width
        MagField.width = str2num(regexprep(parameters((strmatch('XWID',parameters)),:),'XWID',''));
        
        % Mag number of data points
        MagField.points = str2num(regexprep(parameters((strmatch('XPTS',parameters)),:),'XPTS',''));
        
        
        MagField.max = MagField.min + MagField.width;
        MagField.sampling = MagField.width / MagField.points;
        
        for k = 1:MagField.points
            x(k) = MagField.min + (MagField.sampling * ((k)-1));
        end
        
        x = x';
        
end



%                         Output arguments
% ========================================================================

% Output results according to requests
switch nargout
    case 0
        graph = 'plot';
    case 1
        varargout{1} = y;
    case 2
        varargout{1} = x;
        varargout{2} = y;
    case 3
        switch extension
            case '.spc'
                info.MWFQ = str2num(regexprep(parameters((strmatch('MF' ,parameters)),:),'MF' ,''));
                info.MWPW = str2num(regexprep(parameters((strmatch('MP' ,parameters)),:),'MP' ,''));
                
                info.CenterField    = str2num(regexprep(parameters((strmatch('HCF',parameters)),:),'HCF',''));
                info.SweepWidth     = str2num(regexprep(parameters((strmatch('HSW',parameters)),:),'HSW',''));
                
            case '.DTA'
                
                % create info parameters, only done if requested to minimise
                % computation time of function
                info.DSRC = strtrim(regexprep(parameters((strmatch('DSRC',parameters)),:),'DSRC',''));
                info.BSEQ = strtrim(regexprep(parameters((strmatch('BSEQ',parameters)),:),'BSEQ',''));
                info.IKKF = strtrim(regexprep(parameters((strmatch('IKKF',parameters)),:),'IKKF',''));
                info.XTYP = strtrim(regexprep(parameters((strmatch('XTYP',parameters)),:),'XTYP',''));
                info.YTYP = strtrim(regexprep(parameters((strmatch('YTYP',parameters)),:),'YTYP',''));
                info.ZTYP = strtrim(regexprep(parameters((strmatch('ZTYP',parameters)),:),'ZTYP',''));
                
                info.IRFMT = strtrim(regexprep(parameters((strmatch('IRFMT',parameters)),:),'IRFMT',''));
                
                info.XMIN = str2num(regexprep(parameters((strmatch('XMIN',parameters)),:),'XMIN',''));
                info.XWID = str2num(regexprep(parameters((strmatch('XWID',parameters)),:),'XWID',''));
                info.XPTS = str2num(regexprep(parameters((strmatch('XPTS',parameters)),:),'XPTS',''));
                
                % info.TITL = strtrim(regexprep(parameters((strmatch('TITL',parameters)),:),'TITL',''));
                info.IRUNI = strtrim(regexprep(parameters((strmatch('IRUNI',parameters)),:),'IRUNI',''));
                
                info.OPER = strtrim(regexprep(parameters((strmatch('OPER',parameters)),:),'OPER',''));
                info.DATE = strtrim(regexprep(parameters((strmatch('DATE',parameters)),:),'DATE',''));
                info.TIME = strtrim(regexprep(parameters((strmatch('TIME',parameters)),:),'TIME',''));
                info.CMNT = strtrim(regexprep(parameters((strmatch('CMNT',parameters)),:),'CMNT',''));
                info.SAMP = strtrim(regexprep(parameters((strmatch('SAMP',parameters)),:),'SAMP',''));
                info.SFOR = strtrim(regexprep(parameters((strmatch('SFOR',parameters)),:),'SFOR',''));
                info.STAG = strtrim(regexprep(parameters((strmatch('STAG',parameters)),:),'STAG',''));
                info.OXS1 = strtrim(regexprep(parameters((strmatch('OXS1',parameters)),:),'OXS1',''));
                info.AXS1 = strtrim(regexprep(parameters((strmatch('AXS1',parameters)),:),'AXS1',''));
                info.AXS2 = strtrim(regexprep(parameters((strmatch('AXS2',parameters)),:),'AXS2',''));
                info.AXS3 = strtrim(regexprep(parameters((strmatch('AXS3',parameters)),:),'AXS3',''));
                info.A1CT = str2num(regexprep(parameters((strmatch('A1CT',parameters)),:),'A1CT',''));
                info.A1SW = str2num(regexprep(parameters((strmatch('A1SW',parameters)),:),'A1SW',''));
                info.MWFQ = str2num(regexprep(parameters((strmatch('MWFQ',parameters)),:),'MWFQ',''));
                info.MWPW = str2num(regexprep(parameters((strmatch('MWPW',parameters)),:),'MWPW',''));
                info.AVGS = str2num(regexprep(parameters((strmatch('AVGS',parameters)),:),'AVGS',''));
                
                info.CenterField = strtrim(regexprep(parameters((strmatch('CenterField',parameters)),:),'CenterField',''));
                info.FieldFlyback = strtrim(regexprep(parameters((strmatch('FieldFlyback',parameters)),:),'FieldFlyback',''));
                info.FieldWait = strtrim(regexprep(parameters((strmatch('FieldWait',parameters)),:),'FieldWait',''));
                info.GFactor = str2num(regexprep(parameters((strmatch('GFactor',parameters)),:),'GFactor',''));
                info.SetToSampleG = strtrim(regexprep(parameters((strmatch('SetToSampleG',parameters)),:),'SetToSampleG',''));
                info.SweepDirection = strtrim(regexprep(parameters((strmatch('SweepDirection',parameters)),:),'SweepDirection',''));
                
                info.FrequencyMon = strtrim(regexprep(parameters((strmatch('FrequencyMon',parameters)),:),'FrequencyMon',''));
                info.QMonitBridge = strtrim(regexprep(parameters((strmatch('QMonitBridge',parameters)),:),'QMonitBridge',''));
                % Power is different because strmatch for 'Power' finds Power and
                % PowerAtten
                info.Power = parameters((strmatch('Power',parameters)),:);
                info.Power = strtrim(regexprep((info.Power(1,:)),'Power',''));
                info.PowerAtten = strtrim(regexprep(parameters((strmatch('PowerAtten',parameters)),:),'PowerAtten',''));
                info.BaselineCorr = strtrim(regexprep(parameters((strmatch('BaselineCorr',parameters)),:),'BaselineCorr',''));
                info.NbScansAcc = str2num(regexprep(parameters((strmatch('NbScansAcc',parameters)),:),'NbScansAcc',''));
                info.NbScansDone = str2num(regexprep(parameters((strmatch('NbScansDone',parameters)),:),'NbScansDone',''));
                info.NbScansToDo = str2num(regexprep(parameters((strmatch('NbScansToDo',parameters)),:),'NbScansToDo',''));
                
                
                % Experiment dependent parameters
                switch info.EXPT
                    case 'PLS'
                        % documentation text
                        info.IIFMT = strtrim(regexprep(parameters((strmatch('IIFMT',parameters)),:),'IIFMT',''));
                        info.IIUNI = strtrim(regexprep(parameters((strmatch('IIUNI',parameters)),:),'IIUNI',''));
                        
                        info.Delay = parameters((strmatch('Delay',parameters)),:);
                        info.Delay = strtrim(regexprep((info.Delay(1,:)),'Delay',''));
                        % cw bridge
                        info.AcqFineTuning = strtrim(regexprep(parameters((strmatch('AcqFineTuning',parameters)),:),'AcqFineTuning',''));
                        info.AcqScanFTuning = strtrim(regexprep(parameters((strmatch('AcqScanFTuning',parameters)),:),'AcqScanFTuning',''));
                        info.AcqSliceFTuning = strtrim(regexprep(parameters((strmatch('AcqSliceFTuning',parameters)),:),'AcqSliceFTuning',''));
                        % endor
                        info.EIEENDORFreq = strtrim(regexprep(parameters((strmatch('EIEENDORFreq',parameters)),:),'EIEENDORFreq',''));
                        info.EIEIsotope = strtrim(regexprep(parameters((strmatch('EIEIsotope',parameters)),:),'EIEIsotope',''));
                        info.EIERFSweepDir = strtrim(regexprep(parameters((strmatch('EIERFSweepDir',parameters)),:),'EIERFSweepDir',''));
                        info.EIEStaticField = strtrim(regexprep(parameters((strmatch('EIEStaticField',parameters)),:),'EIEStaticField',''));
                        info.EIEStaticRF = strtrim(regexprep(parameters((strmatch('EIEStaticRF',parameters)),:),'EIEStaticRF',''));
                        info.ENDORType = strtrim(regexprep(parameters((strmatch('ENDORType',parameters)),:),'ENDORType',''));
                        info.RF1Atten = strtrim(regexprep(parameters((strmatch('RF1Atten',parameters)),:),'RF1Atten',''));
                        info.RF1FreqPos = strtrim(regexprep(parameters((strmatch('RF1FreqPos',parameters)),:),'RF1FreqPos',''));
                        info.RF1StartFreq = strtrim(regexprep(parameters((strmatch('RF1StartFreq',parameters)),:),'RF1StartFreq',''));
                        info.RF1SweepWidth = strtrim(regexprep(parameters((strmatch('RF1SweepWidth',parameters)),:),'RF1SweepWidth',''));
                        info.RF2Atten = strtrim(regexprep(parameters((strmatch('RF2Atten',parameters)),:),'RF2Atten',''));
                        info.RF2FreqPos = strtrim(regexprep(parameters((strmatch('RF2FreqPos',parameters)),:),'RF2FreqPos',''));
                        info.RF2StartFreq = strtrim(regexprep(parameters((strmatch('RF2StartFreq',parameters)),:),'RF2StartFreq',''));
                        info.RF2SweepWidth = strtrim(regexprep(parameters((strmatch('RF2SweepWidth',parameters)),:),'RF2SweepWidth',''));
                        info.RFSrcMixing = strtrim(regexprep(parameters((strmatch('RFSrcMixing',parameters)),:),'RFSrcMixing',''));
                        info.SumAtten = parameters((strmatch('SumAtten',parameters)),:);
                        info.SumAtten = strtrim(regexprep((info.Power(1,:)),'SumAtten',''));
                        info.SumAttenStart = strtrim(regexprep(parameters((strmatch('SumAttenStart',parameters)),:),'SumAttenStart',''));
                        info.SumAttenWidth = strtrim(regexprep(parameters((strmatch('SumAttenWidth',parameters)),:),'SumAttenWidth',''));
                        % fieldCtrl
                        info.FieldResol = str2num(regexprep(parameters((strmatch('FieldResol',parameters)),:),'FieldResol',''));
                        % ftBridge
                        info.Attenuation = strtrim(regexprep(parameters((strmatch('Attenuation',parameters)),:),'Attenuation',''));
                        info.ELDORAtt = strtrim(regexprep(parameters((strmatch('ELDORAtt',parameters)),:),'ELDORAtt',''));
                        info.FrequencyA = strtrim(regexprep(parameters((strmatch('FrequencyA',parameters)),:),'FrequencyA',''));
                        info.VideoBW = strtrim(regexprep(parameters((strmatch('VideoBW',parameters)),:),'VideoBW',''));
                        info.VideoGain = strtrim(regexprep(parameters((strmatch('VideoGain',parameters)),:),'VideoGain',''));
                        % ftEpr
                        info.AutoTimeOut = strtrim(regexprep(parameters((strmatch('AutoTimeOut',parameters)),:),'AutoTimeOut',''));
                        info.AveragesPerScan = str2num(regexprep(parameters((strmatch('AveragesPerScan',parameters)),:),'AveragesPerScan',''));
                        info.ELDORFreqStart = strtrim(regexprep(parameters((strmatch('ELDORFreqStart',parameters)),:),'ELDORFreqStart',''));
                        info.ELDORFreqWidth = strtrim(regexprep(parameters((strmatch('ELDORFreqWidth',parameters)),:),'ELDORFreqWidth',''));
                        info.FTAcqModeSlct = strtrim(regexprep(parameters((strmatch('FTAcqModeSlct',parameters)),:),'FTAcqModeSlct',''));
                        info.PPExtTrg = parameters((strmatch('PPExtTrg',parameters)),:);
                        info.PPExtTrg = strtrim(regexprep((info.PPExtTrg(1,:)),'PPExtTrg',''));
                        info.PPExtTrgSlope = strtrim(regexprep(parameters((strmatch('PPExtTrgSlope',parameters)),:),'PPExtTrgSlope',''));
                        info.PlsSPELEXPSlct = strtrim(regexprep(parameters((strmatch('PlsSPELEXPSlct',parameters)),:),'PlsSPELEXPSlct',''));
                        info.ReplaceMode = parameters((strmatch('ReplaceMode',parameters)),:);
                        info.ReplaceMode_ftEpr = strtrim(regexprep((info.ReplaceMode(1,:)),'ReplaceMode',''));
                        info.ReplaceMode_recorder = strtrim(regexprep((info.ReplaceMode(2,:)),'ReplaceMode',''));
                        clear info.ReplaceMode
                        info.ShotRepTime = strtrim(regexprep(parameters((strmatch('ShotRepTime',parameters)),:),'ShotRepTime',''));
                        info.ShotsPLoop = strtrim(regexprep(parameters((strmatch('ShotsPLoop',parameters)),:),'ShotsPLoop',''));
                        info.SptProgress = strtrim(regexprep(parameters((strmatch('SptProgress',parameters)),:),'SptProgress',''));
                        info.SweepsPExp = strtrim(regexprep(parameters((strmatch('SweepsPExp',parameters)),:),'SweepsPExp',''));
                        info.TriggerTimeOut = strtrim(regexprep(parameters((strmatch('TriggerTimeOut',parameters)),:),'TriggerTimeOut',''));
                        info.XAxisQuant = strtrim(regexprep(parameters((strmatch('XAxisQuant',parameters)),:),'XAxisQuant',''));
                        info.XSpecRes = str2num(regexprep(parameters((strmatch('XSpecRes',parameters)),:),'XSpecRes',''));
                        info.YAxisQuant = strtrim(regexprep(parameters((strmatch('YAxisQuant',parameters)),:),'YAxisQuant',''));
                        info.YSpecRes = str2num(regexprep(parameters((strmatch('YSpecRes',parameters)),:),'YSpecRes',''));
                        
                        
                    case 'CW'
                        % document text
                        info.YFMT = strtrim(regexprep(parameters((strmatch('YFMT',parameters)),:),'YFMT',''));
                        info.YNAM = strtrim(regexprep(parameters((strmatch('YNAM',parameters)),:),'YNAM',''));
                        info.YUNI = strtrim(regexprep(parameters((strmatch('YUNI',parameters)),:),'YUNI',''));
                        % standard parameter
                        info.RESO = strtrim(regexprep(parameters((strmatch('RESO',parameters)),:),'RESO',''));
                        info.SPTP = str2num(regexprep(parameters((strmatch('SPTP',parameters)),:),'SPTP',''));
                        info.RCAG = str2num(regexprep(parameters((strmatch('RCAG',parameters)),:),'RCAG',''));
                        info.RCHM = str2num(regexprep(parameters((strmatch('RCHM',parameters)),:),'RCHM',''));
                        info.B0MA = str2num(regexprep(parameters((strmatch('B0MA',parameters)),:),'B0MA',''));
                        info.B0MF = str2num(regexprep(parameters((strmatch('B0MF',parameters)),:),'B0MF',''));
                        info.RCPH = str2num(regexprep(parameters((strmatch('RCPH',parameters)),:),'RCPH',''));
                        info.RCOF = str2num(regexprep(parameters((strmatch('RCOF',parameters)),:),'RCOF',''));
                        info.A1RS = str2num(regexprep(parameters((strmatch('A1RS',parameters)),:),'A1RS',''));
                        info.RCTC = str2num(regexprep(parameters((strmatch('RCTC',parameters)),:),'RCTC',''));
                        % mwBridge
                        info.AcqFineTuning = strtrim(regexprep(parameters((strmatch('AcqFineTuning',parameters)),:),'AcqFineTuning',''));
                        % fieldCtrl
                        %            info.SweepWidth = parameters((strmatch('SweepWidth',parameters)),:);
                        %            info.SweepWidth_fieldCtrl = strtrim(regexprep((info.SweepWidth(1,:)),'SweepWidth',''));
                        %            info.SweepWidth_ramp2 = strtrim(regexprep((info.SweepWidth(2,:)),'SweepWidth',''));
                        clear info.SweepWidth
                        info.Delay = strtrim(regexprep(parameters((strmatch('Delay',parameters,'exact')),:),'Delay',''));
                        info.ReplaceMode = strtrim(regexprep(parameters((strmatch('ReplaceMode',parameters)),:),'ReplaceMode',''));
                        % signalChannel
                        info.AFCTrap = strtrim(regexprep(parameters((strmatch('AFCTrap',parameters)),:),'AFCTrap',''));
                        info.Calibrated = strtrim(regexprep(parameters((strmatch('Calibrated',parameters)),:),'Calibrated',''));
                        info.ConvTime = strtrim(regexprep(parameters((strmatch('ConvTime',parameters)),:),'ConvTime',''));
                        info.DModAFCTrap = strtrim(regexprep(parameters((strmatch('DModAFCTrap',parameters)),:),'DModAFCTrap',''));
                        info.DModAmp = strtrim(regexprep(parameters((strmatch('DModAmp',parameters)),:),'DModAmp',''));
                        info.DModCalibrated = strtrim(regexprep(parameters((strmatch('DModCalibrated',parameters)),:),'DModCalibrated',''));
                        info.DModDetectSCT = strtrim(regexprep(parameters((strmatch('DModDetectSCT',parameters)),:),'DModDetectSCT',''));
                        info.DModEliDelay = strtrim(regexprep(parameters((strmatch('DModEliDelay',parameters)),:),'DModEliDelay',''));
                        info.DModExtLockIn = strtrim(regexprep(parameters((strmatch('DModExtLockIn',parameters)),:),'DModExtLockIn',''));
                        info.DModExtTrigger = strtrim(regexprep(parameters((strmatch('DModExtTrigger',parameters)),:),'DModExtTrigger',''));
                        info.DModFieldMod = strtrim(regexprep(parameters((strmatch('DModFieldMod',parameters)),:),'DModFieldMod',''));
                        info.DModGain = strtrim(regexprep(parameters((strmatch('DModGain',parameters)),:),'DModGain',''));
                        info.DModHighPass = strtrim(regexprep(parameters((strmatch('DModHighPass',parameters)),:),'DModHighPass',''));
                        info.DModIntegrator = strtrim(regexprep(parameters((strmatch('DModIntegrator',parameters)),:),'DModIntegrator',''));
                        info.DModModOutput = strtrim(regexprep(parameters((strmatch('DModModOutput',parameters)),:),'DModModOutput',''));
                        info.DModSignalInput = strtrim(regexprep(parameters((strmatch('DModSignalInput',parameters)),:),'DModSignalInput',''));
                        info.DModTimeConst = strtrim(regexprep(parameters((strmatch('DModTimeConst',parameters)),:),'DModTimeConst',''));
                        info.DoubleModFreq = strtrim(regexprep(parameters((strmatch('DoubleModFreq',parameters)),:),'DoubleModFreq',''));
                        info.DoubleModPhase = strtrim(regexprep(parameters((strmatch('DoubleModPhase',parameters)),:),'DoubleModPhase',''));
                        info.DoubleMode = strtrim(regexprep(parameters((strmatch('DoubleMode',parameters)),:),'DoubleMode',''));
                        info.EliDelay = strtrim(regexprep(parameters((strmatch('EliDelay',parameters)),:),'EliDelay',''));
                        info.ExtLockIn = strtrim(regexprep(parameters((strmatch('ExtLockIn',parameters)),:),'ExtLockIn',''));
                        info.ExtTrigger = strtrim(regexprep(parameters((strmatch('ExtTrigger',parameters)),:),'ExtTrigger',''));
                        info.Gain = strtrim(regexprep(parameters((strmatch('Gain',parameters)),:),'Gain',''));
                        info.Harmonic = str2num(regexprep(parameters((strmatch('Harmonic',parameters)),:),'Harmonic',''));
                        info.HighPass = strtrim(regexprep(parameters((strmatch('HighPass',parameters)),:),'HighPass',''));
                        info.Integrator = strtrim(regexprep(parameters((strmatch('Integrator',parameters)),:),'Integrator',''));
                        info.ModAmp = strtrim(regexprep(parameters((strmatch('ModAmp',parameters)),:),'ModAmp',''));
                        info.ModFreq = strtrim(regexprep(parameters((strmatch('ModFreq',parameters)),:),'ModFreq',''));
                        info.ModInput = strtrim(regexprep(parameters((strmatch('ModInput',parameters)),:),'ModInput',''));
                        info.ModOutput = strtrim(regexprep(parameters((strmatch('ModOutput',parameters)),:),'ModOutput',''));
                        info.ModPhase = str2num(regexprep(parameters((strmatch('ModPhase',parameters)),:),'ModPhase',''));
                        info.Offset = strtrim(regexprep(parameters((strmatch('Offset',parameters)),:),'Offset',''));
                        info.QuadMode = strtrim(regexprep(parameters((strmatch('QuadMode',parameters)),:),'QuadMode',''));
                        info.QuadPhase = str2num(regexprep(parameters((strmatch('QuadPhase',parameters)),:),'QuadPhase',''));
                        info.Resolution = str2num(regexprep(parameters((strmatch('Resolution',parameters)),:),'Resolution',''));
                        info.Resonator = str2num(regexprep(parameters((strmatch('Resonator',parameters)),:),'Resonator',''));
                        info.SctNorm = strtrim(regexprep(parameters((strmatch('SctNorm',parameters)),:),'SctNorm',''));
                        info.SctRevision = strtrim(regexprep(parameters((strmatch('SctRevision',parameters)),:),'SctRevision',''));
                        info.SignalInput = strtrim(regexprep(parameters((strmatch('SignalInput',parameters)),:),'SignalInput',''));
                        info.SweepTime = strtrim(regexprep(parameters((strmatch('SweepTime',parameters)),:),'SweepTime',''));
                        info.TimeConst = strtrim(regexprep(parameters((strmatch('TimeConst',parameters)),:),'TimeConst',''));
                        info.TuneCaps = str2num(regexprep(parameters((strmatch('TuneCaps',parameters)),:),'TuneCaps',''));
                end
        end
        
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = info;
    otherwise
        varargout{1} = x;
        varargout{2} = y;
end

%                         Figure plotting
% ========================================================================

% Check if figure plotting is selected
if isequal(graph,'plot')
    figure('name' , ['BrukerRead: ' name], 'NumberTitle','off');
    
    % PULSE EXPERIMENTS
    switch info.EXPT
        case 'PLS'
            switch info.XNAM
                    % For FSEs
                case '''Field'''
                    plot(x/10,y.real,'k-' , x/10,y.imag,'r-');
                    xlabel('Magnetic Field / mT');
                    
                    % For FIDs
                case '''Time'''
                    plot(x,y.real,'k-' , x,y.imag,'r-');
                    xlabel('Time / ns');
                    
                % case '''Frequency'''
                    % PLACEHOLDER
                    
                    % Try plotting for other experiments (ie ENDOR)
                otherwise
                    try
                        plot(x,y.real,'k-' , x,y.imag,'r-');
                    end
            end
            
            legend('Real','Imaginary')
            
    % CW EXPERIMENTS
        case 'CW'
            plot(x/10,y);
            xlabel('Magnetic Field / mT');
    end
    
    % Figure formatting
    ylabel('Intensity');
    set(gcf,'color', 'white');
    axis tight;
    
end

end
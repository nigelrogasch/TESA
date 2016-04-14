% tesa_findpulsepeak() - finds TMS pulses by detecting the large TMS artifacts
%                   peaks in the data. This script works by extracting a
%                   single channel and finding the time points in which a peak
%                   above a certain threshold is detected. Different
%                   methods for setting the threshold can be determined and
%                   either the positive or negative peak used (or an
%                   interactive gui which allows the user to select which
%                   peaks are included).
%                   Paired pulses and repetitive TMS trains can
%                   also be deteceted.
%                   This script is an alternative to tesa_findpulse
% 
%                   Note that this script requires a toolbox developed by
%                   Daniel Wagenaar (WagenaarMBL) which is available for
%                   download at: http://www.its.caltech.edu/~daw/teach.html
%                   Once the toolbox is unzipped, please make sure this is
%                   added to your matlab path.
%
% Usage:
%   >>  EEG = tesa_findpulsepeak( EEG, elec, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string with electrode to use for finding artifact
% 
% Optional input pairs:
%   'dtrnd','str'  'poly'|'gradient'|'median'|'linear'|'none'. Defines the type of detrend used
%                   to centre the data, to generate an analytic signal.
%                   default = poly
%   'thrshtype','str'/int - 'dynamic'|'median'|value. Defines the type of
%                   threshold used to determine peaks. Dynamic sets threshold 
%                   to the range points above/below 99.9 percent of data
%                   trace. Median sets threshold as median of points above/below 
%                   99.9 percent of data trace. Value is a user defined integer 
%                   for setting the threshold (in uV). (e.g. 1000)
%                   default = 'dymanic'
%   'wpeaks','str' - 'pos'|'neg'|'gui'. Defines whether to use the
%                   positive or negative peak to define the artifact, or to
%                   use an interactive GUI. For the GUI, a box is generated
%                   by clicking and dragging the mouse cursor. Peaks within
%                   the box will be included as events.
%                   default = 'pos'
%   'plots','str' - 'on'|'off'. Brings up a plot showing the detected
%                   peaks. Black = detected, pink = selected for
%                   definition.
%                   default = 'on'
%   'tmsLabel','str'- 'str' is a string for the single TMS label.  
%                   default = 'TMS'
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - required. 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
%   'ISI', [int]    - required. [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [1,2,...]. 
%                   default = []
%   'pairLabel',{'str'} - required if more than 1 ISI. {'str'} is a cell array
%                   containing string labels for different ISI conditions.  
%                   Multiple labels can be defined as {'SICI','LICI',...}.
%                   The number of labels defined must equal the number of
%                   ISI conditions defined.
%                   default = {'TMSpair'}
% 
%  Input pairs for detecting repetitive TMS trains
%  'repetitive','str' - required. 'str' - type 'yes' to turn on repetitive detection
%                   default = 'no'
%   'ITI', int      - required. int defines the inter-train interval in ms.
%                   For example, if a 10 Hz rTMS condition is used with 4s
%                   of stimulation (40 pulses) and 26s of rest, ITI = 2600;
%                   default = []
%   'pulseNum', int - required. int defines the number of pulses in a
%                   train. Using the above example, this would be 40. 
%                   deafult = []
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_findpulsepeak( EEG, 'Cz' ); %default use
%   EEG = tesa_findpulsepeak( EEG, 'Fz', 'dtrnd', 'linear', 'thrshtype', 'median', 'wpeaks', 'gui', 'plots', 'off', 'tmsLabel', 'single' ); %user defined
%   EEG = tesa_findpulsepeak( EEG, 'Cz', 'paired', 'yes', 'ISI', [100],'pairLabel', {'LICI'}); %paired pulse use
%   EEG = tesa_findpulsepeak( EEG, 'Cz', 'repetitive', 'yes', 'ITI', 26, 'pulseNum', 40 ); %rTMS use 
%
% See also:
%   tesa_findpulse, tesa_fixevent 

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% Authors:
% Caley Sullivan, Monash University, calley.sullivan@monash.edu
% 
% Based on functions developed by Daniel Wagenaar 
%                                  http://www.its.caltech.edu/~daw/teach.html
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function EEG = tesa_findpulsepeak( EEG, elec, varargin )

if nargin < 2
	error('Not enough input arguments.');
end

%define defaults
options = struct('dtrnd','poly','thrshtype','dynamic','wpeaks','pos','plots','on','tmsLabel','TMS','paired','no','ISI',[],'pairLabel',{'TMSpair'},'repetitive','no','ITI',[],'pulseNum',[]);
options.pairLabel = cellstr(options.pairLabel);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs key/value pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = pair{1}; % make case insensitive

   if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

%check that data is continuous, not epoched
if size(size(EEG.data),2) > 2
    error('tesa_findpulsepeak only works on continuous data. Use tesa_fixtrigger for epoched data.')
end

%check that paired and repetitive have been correctly called
if ~(strcmp(options.paired,'no') || strcmp(options.paired,'yes'))
    error('paired must be either ''yes'' or ''no''.');
end
if ~(strcmp(options.repetitive,'no') || strcmp(options.repetitive,'yes'))
    error('repetitive must be either ''yes'' or ''no''.');
end

%checks that the WagenaarMBL toolbox is added to the matlab path
if exist('detectspike','file') ~= 2
    error('Please add the Wagenaar MBL toolbox to your matlab path. See help or the TESA manual for details.');
end

%finds channel for thresholding
for z = 1:EEG.nbchan;
    chan{1,z} = EEG.chanlocs(1,z).labels;
end;
num = find(strcmpi(elec,chan));%defines row number of channel used for thresholding

%check that channel for thresholding exists
if isempty(num)
    error('Electrode not found. Please enter a channel that is present.');
end

%Extracts channel
signal = EEG.data(num,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINDS PEAKS OF TMS ARTIFACTS
% Detrend the data
if strcmp(options.dtrnd,'poly')
    %baseline detrend  useing polyfit
    [p3,~,mu3] = polyfit((1:numel(signal'))',signal',6);
    f_y3 = polyval(p3,(1:numel(signal'))',[],mu3);
    sig  = signal'-f_y3;       
elseif strcmp(options.dtrnd,'gradient')
    % use the 1st derivative of the signal
    sig =gradient(signal); 
elseif strcmp(options.dtrnd,'median')
    sig = signal - median(signal) ; 
elseif strcmp(options.dtrnd,'linear')
    sig = detrend(signal,'linear',10) ; 
    % Y = detrend(X,'linear',BP) removes a continuous, piecewise linear trend.
    %     Breakpoint indices for the linear trend are contained in the vector BP.
    %     The default is no breakpoints, such that one single straight line is
    %     removed from each column of X.
else 
   sig=signal; 
end
clear signal; 

% Defines theshold for artifact peaks
if  strcmp(options.thrshtype,'dynamic'); 
% Set threshold to the range points above/below 99.9% of data 
    prcnt= prctile(sig,[0.1 99.9]);
    dynamic_thrsh=min(abs(  [range((sig(sig>[prcnt(2)])))/10,range((sig(sig<[prcnt(1)])))/10])); 
    thrsh = dynamic_thrsh; 
elseif strcmp(options.thrshtype,'median'); 
% Set threshold as median of points above/below 99.9% of data 
    prcnt= prctile(sig,[0.1 99.9]);
    med_thrsh = min(abs(   [median(sig(sig>[prcnt(2)])) ,  median(sig(sig<[prcnt(1)]))]  ) );
    thrsh= med_thrsh; 
elseif  isnumeric(options.thrshtype); 
% set threshold to ??????uV
    static_thrsh = options.thrshtype;
    thrsh = static_thrsh;  
end ;
 
dat = (double(sig));
tms = (1:numel(sig))';
spk = detectspike((dat),tms,thrsh,1*EEG.srate,10*EEG.srate);
artlat= squeeze(spk.tms);


if strcmp(options.wpeaks,'gui'); 
    options.plots='off';
end 
%Sanity plot
if strcmp(options.plots,'on'); 
    figure; 
    plot(tms,dat,'b');
    hold on;
    plot(spk.tms, spk.amp,'k.');
end; 

% Defines peaks to use
if strcmp(options.wpeaks,'neg') ;
%   select negative peaks
    spk.tms(spk.amp>-thrsh)=[];
    spk.amp(spk.amp>-thrsh)=[];
elseif strcmp(options.wpeaks,'gui'); 
    
global redo
    redo=true;
    %gdspk = selectspike(spk);
while redo==true
    f= figure; 
    h1= plot((tms),(dat),'b');
    hold on;
    h2= plot(spk.tms, spk.amp,'k.');  
    title('Inspect the detected peaks...','BackgroundColor',[0.729411780834198 0.831372559070587 0.95686274766922])
        h = uicontrol('Position', [5 5 150 20], 'String', 'continue with selection', ...
                      'Callback', 'uiresume(gcbf)');
        %disp('This will print immediately');
        uiwait(gcf);
        %disp('This will print after you click Continue');
    title('   Click, Hold & Drag the mousepointer to select which peaks to keep...','BackgroundColor',[0.729411780834198 0.831372559070587 0.95686274766922])
    select = tesa_selectdata('Return','selected','selectionmode','lasso','action','list','ignore',h1 ,'Verify','off' );
    gdspk.tms=spk.tms(select);
    gdspk.amp=spk.amp(select);
    hold on; plot(gdspk.tms, gdspk.amp,'ro');  
ButtonName = questdlg('Are you happy with the peaks selected?', ...
                         'Verify Selection', ...
                         'No-Redo', 'Yes-Continue', 'Cancel', 'No-Redo');
   switch ButtonName,
     case 'No-Redo',
      redo=true;
     case 'Yes-Continue',
      redo=false;
      case 'Cancel',
          redo=false; close all ; return ;
   end % switch
close
end
%       f= figure; 
%     h1= plot((tms(1:end)),(dat(1:end,1)),'b');
%     hold on;
%     h2= plot(spk.tms, spk.amp,'k.');
%     h3= plot(spk.tms(select), spk.amp(select),'g.');   
    spk= gdspk; 

elseif strcmp(options.wpeaks,'pos');
    %select positive peaks
    spk.tms(spk.amp<thrsh)=[];
    spk.amp(spk.amp<thrsh)=[];
end 

if strcmp(options.plots,'on'); 
    hold on;
    h2= plot(spk.tms, spk.amp,'m.');
end; 
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimAll = spk.tms';

%Trims any white space from the single label
options.tmsLabel = strtrim(options.tmsLabel);

%Create label master 
for a = 1:size(stimAll,2)
    stimLabel{1,a} = options.tmsLabel;
end

%Marks conditioning pulses in paired pulse paradigms
if strcmp(options.paired,'yes')
    
    %Check that ISI has been provided
    if isempty(options.ISI)
        error('Please provide the interstimulus interval (ISI) for detecting paired pulse.');
    end
    
    %If more than one ISI provided, check that an equal number of labels
    %have been provided
    if size(options.ISI,2) > 1
        if size(options.ISI,2) ~= size(options.pairLabel,2)
            error('Number of paired labels must equal the number of interstimulus intervals provided. Use ''pairLabel'',{}, in function to provide names')
        end
    end
    
    %Check that refractory period is less that the ISI
%     if options.refract > options.ISI
%         error('The refractory period is shorter than the interstimulus interval. This will result in inaccurate detection of the test pulse. Script termninated.');
%     end
    
    %Check that single and paired labels are unique
    for a = 1:size(options.pairLabel,2)
        if strcmp(options.tmsLabel,options.pairLabel{1,a})
            error('Paired label is the same as single label. Please ensure that labels are unique.');
        end
    end
    
    %Check that paired labels are unique
    if ~(size(unique(options.pairLabel),2) == size(options.pairLabel,2))
        error('Paired labels are not unique. Please ensure each label is different.')
    end
    
    %Trims any white space from the label names
    options.pairLabel = strtrim(options.pairLabel);
    
    for a = 1:size(options.ISI,2)
        sISI = ceil(EEG.srate./1000.*options.ISI(1,a)); %converts ISI to samples
        prec = ceil(EEG.srate./1000.*0.5); %precision for searching for second pulse = +/-0.5 ms
        diffStim = diff(stimAll);
        for b = 1:size(diffStim,2)
            if diffStim(1,b) > sISI-prec && diffStim(1,b) < sISI+prec
                stimLabel{1,b} = 'con';
                stimLabel{1,b+1} = options.pairLabel{1,a};
            end
        end
    end
end

%Marks repetitive pulses for rTMS
if strcmp(options.repetitive,'yes')
    
    %Check that ITI has been provided
    if isempty(options.ITI)
        error('Please provide the inter-train interval (ITI - the time between trains of pulses). The ITI is in ms.');
    end
    
    %Check that pulseNum has been provided
    if isempty(options.pulseNum)
        error('Please provide the number of pulses in a train (pulseNum).');
    end
    
    %check if the number of pulses is divisible by the train number given
    testDiv = size(stimAll,2)/options.pulseNum;
    if ~floor(testDiv) == testDiv
        warning('The number of pulses in a train is not divisible by the total number of pulses. There may be some incorrectly identified pulses or the number of pulses in a train may be incorrect.');
    end
    
    sITI = ceil(EEG.srate./1000.*options.ITI); %converts ITI to samples
    prec = ceil(EEG.srate./1000.*5); %precision for searching for second pulse = +/-5 ms
    diffStim = diff([stimAll(1,1)-options.ITI stimAll]);
    for a = 1:size(diffStim,2)
        if diffStim(1,a) > options.ITI-prec
            if stimAll(1,a+(options.pulseNum-1))-stimAll(1,a) > diffStim(1,a+1).*(options.pulseNum-1)+10
                error('Number of pulses in detected train are different from that specified. Script terminated');
            else
                for b = 1:options.pulseNum
                    stimLabel{1,a+(b-1)} = ['TMS',num2str(b)];
                end
            end
        end
    end
end

%Create/insert new events in to event and urevent fields
if isempty(EEG.event)
    for a=1:size(stimAll,2);
        EEG.event(1,a).type=stimLabel{1,a};
        EEG.event(1,a).latency=stimAll(1,a);
        EEG.event(1,a).urevent=a;
        EEG.urevent(1,a).type=stimLabel{1,a};
        EEG.urevent(1,a).latency=stimAll(1,a);
    end;
else
    n=size(EEG.event,2);
    for a=1:size(stimAll,2);
        EEG.event(1,a+n).type=stimLabel{1,a};
        EEG.event(1,a+n).latency=stimAll(1,a);
        EEG.event(1,a+n).urevent=a;
        EEG.urevent(1,a+n).type=stimLabel{1,a};
        EEG.urevent(1,a+n).latency=stimAll(1,a);
    end;
end

%Count the number of each stimuli and display
uniCount = sum(strcmp(options.tmsLabel,stimLabel));
fprintf('%d single pulses detected and labelled as ''%s''.\n',uniCount,options.tmsLabel);

if strcmp(options.paired,'yes')
    for a = 1:size(options.pairLabel,2)
        uniCount(1,a) = sum(strcmp(options.pairLabel{1,a},stimLabel));
        if uniCount(1,a) == 0
            warning('No pulses for condition ''%s'' were detected.',options.pairLabel{1,a})
        else
            fprintf('%d paired test pulses detected and labelled as ''%s''. Conditioning pulses labelled as ''con''.\n',uniCount(1,a),options.pairLabel{1,a});
        end
    end
end

if strcmp(options.repetitive,'yes')
    uniLabel = unique(stimLabel);
    uniCount = sum(strcmp(uniLabel{1,1},stimLabel));   
    fprintf('%d repetitive trains detected with %d pulses in each train. The first pulse of each train labelled as ''%s''.\n',uniCount,size(uniLabel,2),uniLabel{1,1});
end

end
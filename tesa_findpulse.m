% tesa_findpulse() - finds TMS pulses by detecting the large TMS artifacts
%                   present in the data. This script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses and repetitive TMS trains can
%                   also be deteceted. 
%
% Usage:
%   >>  EEG = tesa_findpulse( EEG, elec );
%   >>  EEG = tesa_findpulse( EEG, elec, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string defining electrode to use for finding artifact
%                   example: 'Cz'
% 
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artifact to recover below the rate of change). int
%                   is an integer in ms. 
%                   default = 3
% 
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e4
% 
%   'tmsLabel','str'- 'str' is a string for the single TMS label.  
%                   default = 'TMS'
%
%   'plots','str' - 'on'|'off'. Brings up a plot showing the detected
%                   peaks. Red = detected 
%                   default = 'on'
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - required. 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
% 
%   'ISI', [int]    - required. [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [3,100,...]. 
%                   default = []
% 
%   'pairLabel',{'str'} - required if more than 1 ISI. {'str'} is a cell array
%                   containing string labels for different ISI conditions.  
%                   Multiple labels can be defined as {'SICI','LICI',...}.
%                   The number of labels defined must equal the number of
%                   ISI conditions defined.
%                   default = {'TMSpair'}
% 
% Input pairs for detecting repetitive TMS trains
%  'repetitive','str' - required. 'str' - type 'yes' to turn on repetitive detection
%                   default = 'no'
% 
%   'ITI', int      - required. int defines the inter-train interval in ms.
%                   For example, if a 10 Hz rTMS condition is used with 4s
%                   of stimulation (40 pulses) and 26s of rest, ITI = 2600;
%                   default = []
% 
%   'pulseNum', int - required. int defines the number of pulses in a
%                   train. Using the above example, this would be 40. 
%                   deafult = []
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = tesa_findpulse( EEG, 'Cz' ); %default use
%   EEG = tesa_findpulse( EEG, 'Fz', 'refract', 4, 'rate', 2e5, 'tmsLabel', 'single','plots','off' ); %user defined input
%   EEG = tesa_findpulse( EEG, 'Cz', 'paired', 'yes', 'ISI', [100],'pairLabel', {'LICI'}); %paired pulse use
%   EEG = tesa_findpulse( EEG, 'Cz', 'repetitive', 'yes', 'ITI', 26, 'pulseNum', 40 ); %rTMS use 
%
% See also:
%   tesa_findpulsepeak, tesa_fixevent

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
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

function EEG = tesa_findpulse( EEG, elec, varargin )

if nargin < 2
	error('Not enough input arguments.');
end

%define defaults
options = struct('refract',3,'rate',1e4,'tmsLabel','TMS','plots','on','paired','no','ISI',[],'pairLabel',{'TMSpair'},'repetitive','no','ITI',[],'pulseNum',[]);
options.pairLabel = cellstr(options.pairLabel);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('Script needs key/value pairs')
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
    error('tesa_findpulse only works on continuous data. Use tesa_fixtrigger for epoched data.')
end

%check that paired and repetitive have been correctly called
if ~(strcmp(options.paired,'no') || strcmp(options.paired,'yes'))
    error('paired input must be either ''yes'' or ''no''.');
end
if ~(strcmp(options.repetitive,'no') || strcmp(options.repetitive,'yes'))
    error('repetitive input must be either ''yes'' or ''no''.');
end

%if events are present, check that the EEG.event.type is a string and that
%the event name is different from that specified
if ~isempty(EEG.event)
    for a = 1:size(EEG.event,2)
        if isnumeric(EEG.event(a).type)
            EEG.event(a).type = num2str(EEG.event(a).type);
        end
        if strcmp(EEG.event(a).type,options.tmsLabel)
            error('The label ''%s'' already exists in the data. Please choose a different trigger label.',options.tmsLabel)
        end
    end
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
data = EEG.data(num,:);

%Calculate derivatives
h = (1/EEG.srate)*1000;     %step size in ms
der1 = diff(data)/h;        %calculates first derivative

%finds artifact (defined as first derivative > rate)
rateS = options.rate.*h; %Convert rate in to change in uV per sample
logstim = abs(der1)>rateS;
samp =(1:size(data,2));
stim = samp(logstim);
stim = stim-1; %Makes start of artifact the defining point

%Remove triggers within refractory period
sRef = ceil(EEG.srate./1000.*options.refract); %converts refractory period to samples
refPer = stim(1,1)+sRef; %defines refractory period following stimulus
stimAll(1,1) = stim(1,1);
stimAmp(1,1) = data(1,stim(1,1));
for a = 2:size(stim,2)
    if stim(1,a) > refPer
        stimAll(1,size(stimAll,2)+1) = stim(1,a);
        stimAmp(1,size(stimAmp,2)+1) = data(1,stim(1,a));
        refPer = stim(1,a)+sRef;
    end 
end

%Sanity plot
if strcmp(options.plots,'on'); 
    figure; 
    plot(1:1:EEG.pnts,data(1,1:end),'b');
    hold on;
    plot(stimAll(1,1:end), stimAmp(1,1:end),'r.');
    title(['TMS pulse artifacts identified in electrode ', elec]);
end; 

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
    if options.refract > options.ISI
        error('The refractory period is shorter than the interstimulus interval. This will result in inaccurate detection of the test pulse. The refractory period can be altered using the ''refract'' input.');
    end
    
    %Check that single and paired labels are unique
    for a = 1:size(options.pairLabel,2)
        if strcmp(options.tmsLabel,options.pairLabel{1,a})
            error('Paired label is the same as single label. Please ensure that labels are unique.');
        end
    end
    
    %Check that paired labels are unique
    if size(unique(options.pairLabel),2) ~= size(options.pairLabel,2)
        error('Paired labels are not unique. Please ensure each label is different following ''pairLabel''.')
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
                error('Number of pulses in detected train are different from that specified. Check the ITI input is correct. Note the ITI is in ms.');
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

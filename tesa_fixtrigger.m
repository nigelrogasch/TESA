% tesa_fixtrigger() - finds TMS pulses by detecting the large TMS artifacts
%                   present in already epoched data. This script is
%                   designed for instances when the recoreded triggers do
%                   not correspond with when the TMS pulse was given.
%                   The script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses and repetitive TMS trains can
%                   also be deteceted.
% 
%                   IMPORTANT: If you need to use this function, make sure
%                   that the initial epochs you use are larger than the
%                   final epoch size you desire. If the initial epoch size is
%                   too small, the new epoch window will be out of range of
%                   with the new triggers. E.g. inital epoch is -100 to 100
%                   and the trigger is shifted 10 ms so the new 0 now sits
%                   at 10 ms. Re-epoching the data to -100 to 100 won't
%                   work as the new range is effectively -90 to 110. In this case, 
%                   run the initial epoch at -120 to 120 and the epoch can now
%                   be taken.
%                   
%
% Usage:
%   >>  EEG = tesa_fixtrigger( EEG, elec, newEpoch );
%   >>  EEG = tesa_fixtrigger( EEG, elec, newEpoch, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string with electrode to use for finding artifact
%   newEpoch        - vector with start and end time of new epoch in
%                   seconds (following pop_epoch convention). 
%                   Example: [-1,1] %For -1 s to 1s epoch
% 
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artefact to recover below the rate of change). int
%                   is in ms. 
%                   default = 3
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e5
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
% See also:
%   SAMPLE, EEGLAB 

% Copyright (C) 2015  Nigel Rogasch, Monash University,
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

function EEG = tesa_fixtrigger( EEG, elec, newEpoch, varargin )

if nargin < 2
	error('Not enough input arguments.');
end

%define defaults
options = struct('refract',3,'rate',1e5,'tmsLabel','TMS','paired','no','ISI',[],'pairLabel',{'TMSpair'},'repetitive','no','ITI',[],'pulseNum',[]);
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

%check that paired and repetitive have been correctly called
if ~(strcmp(options.paired,'no') || strcmp(options.paired,'yes'))
    error('paired must be either ''yes'' or ''no''.');
end
if ~(strcmp(options.repetitive,'no') || strcmp(options.repetitive,'yes'))
    error('repetitive must be either ''yes'' or ''no''.');
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
data = EEG.data(num,:,:);
epTime = reshape(1:1:(size(EEG.data,2)*size(EEG.data,3)),[],size(EEG.data,3))';

%Trims any white space from the single label
options.tmsLabel = strtrim(options.tmsLabel);

%For single pulse data
if strcmp(options.paired,'no') && strcmp(options.repetitive,'no')
    for a = 1:size(data,3)

        %Calculate derivatives
        h = (1/EEG.srate)*1000;     %step size in ms
        der1 = diff(data(1,:,a))/h;        %calculates first derivative

        %finds artifact (defined as first derivative > rate)
        rateS = options.rate.*h; %Convert rate in to change in uV per sample
        logstim = abs(der1)>rateS;
        samp =(1:size(data,2));
        stim = samp(logstim);

        %Creates new event info
        n = size(EEG.event,2);
        p = size(EEG.urevent,2);
        EEG.event(1,n+1).type = options.tmsLabel;
        EEG.event(1,n+1).latency = epTime(a,stim(1,1));
        EEG.event(1,n+1).urevent = p+1;
        EEG.urevent(1,p+1).type = options.tmsLabel;
        EEG.urevent(1,p+1).latency = epTime(a,stim(1,1));

    end

    %Re-epochs data
    EEG = pop_epoch( EEG, {options.tmsLabel}, newEpoch, 'epochinfo', 'yes');
end
    
%For paired pulse data
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
            error('The refractory period is shorter than the interstimulus interval. This will result in inaccurate detection of the test pulse. Script termninated.');
        end

        %Check that single and paired labels are unique
        for a = 1:size(options.pairLabel,2)
            if strcmp(options.tmsLabel,options.pairLabel{1,a})
                error('Paired label is the same as single label. Please ensure that labels are unique.');
            end
        end

        %Check that paired labels are unique
        if ~size(unique(options.pairLabel),2) == size(options.pairLabel,2)
            error('Paired labels are not unique. Please ensure each label is different.')
        end

        %Trims any white space from the label names
        options.pairLabel = strtrim(options.pairLabel);

        for x = 1:size(data,3)

            %Calculate derivatives
            h = (1/EEG.srate)*1000;     %step size in ms
            der1 = diff(data(1,:,x))/h;        %calculates first derivative

            %finds artifact (defined as first derivative > rate)
            rateS = options.rate.*h; %Convert rate in to change in uV per sample
            logstim = abs(der1)>rateS;
            samp =(1:size(data,2));
            stim = samp(logstim);

            %Remove triggers within refractory period
            sRef = ceil(EEG.srate./1000.*options.refract); %converts refractory period to samples
            refPer = stim(1,1)+sRef; %defines refractory period following stimulus
            stimAll(1,1) = stim(1,1);
            for a = 2:size(stim,2)
                if stim(1,a) > refPer
                    stimAll(1,size(stimAll,2)+1) = stim(1,a);
                    refPer = stim(1,a)+sRef;
                end 
            end

            %Determine whether trial has single or paired data
            stimLabel = [];
            for a = 1:size(options.ISI,2)
                sISI = ceil(EEG.srate./1000.*options.ISI(1,a)); %converts ISI to samples
                prec = ceil(EEG.srate./1000.*0.5); %precision for searching for second pulse = +/-0.5 ms
                diffStim = diff(stimAll);
                for b = 1:size(diffStim,2)
                    if diffStim(1,b) > sISI-prec && diffStim(1,b) < sISI+prec
                        stimLabel{1,1} = 'con';
                        stimLabel{1,2} = options.pairLabel{1,a};
                    end
                end
            end

            %Creates new event info
            if isempty(stimLabel)
                n = size(EEG.event,2);
                p = size(EEG.urevent,2);
                EEG.event(1,n+1).type = options.tmsLabel;
                EEG.event(1,n+1).latency = epTime(x,stimAll(1,1));
                EEG.event(1,n+1).urevent = p+1;
                EEG.urevent(1,p+1).type = options.tmsLabel;
                EEG.urevent(1,p+1).latency = epTime(x,stimAll(1,1));
                stimKeep{1,x} = options.tmsLabel;
            else
                n = size(EEG.event,2);
                p = size(EEG.urevent,2);
                %Conditioning pulse
                EEG.event(1,n+1).type = stimLabel{1,1};
                EEG.event(1,n+1).latency = epTime(x,stimAll(1,1));
                EEG.event(1,n+1).urevent = p+1;
                EEG.urevent(1,p+1).type = stimLabel{1,1};
                EEG.urevent(1,p+1).latency = epTime(x,stimAll(1,1));
                %Test pulse
                EEG.event(1,n+2).type = stimLabel{1,2};
                EEG.event(1,n+2).latency = epTime(x,stimAll(1,2));
                EEG.event(1,n+2).urevent = p+2;
                EEG.urevent(1,p+2).type = stimLabel{1,2};
                EEG.urevent(1,p+2).latency = epTime(x,stimAll(1,2));
                stimKeep{1,x} = stimLabel{1,2};
            end
        end
        
        %Re-epochs data
        stimNames = unique(stimKeep);
        EEG = pop_epoch( EEG, {stimNames{:}}, newEpoch, 'epochinfo', 'yes');
end

%For repetitive data
if strcmp(options.repetitive,'yes') 
    
    %Check that ITI has been provided
    if isempty(options.ITI)
        error('Please provide the inter-train interval (ITI - the time between trains of pulses). The ITI is in ms.');
    end
    
    %Check that pulseNum has been provided
    if isempty(options.pulseNum)
        error('Please provide the number of pulses in a train (pulseNum).');
    end
    
    for x = 1:size(data,3)
        
        %Calculate derivatives
        h = (1/EEG.srate)*1000;     %step size in ms
        der1 = diff(data(1,:,x))/h;        %calculates first derivative

        %finds artifact (defined as first derivative > rate)
        rateS = options.rate.*h; %Convert rate in to change in uV per sample
        logstim = abs(der1)>rateS;
        samp =(1:size(data,2));
        stim = samp(logstim);

        %Remove triggers within refractory period
        sRef = ceil(EEG.srate./1000.*options.refract); %converts refractory period to samples
        refPer = stim(1,1)+sRef; %defines refractory period following stimulus
        stimAll(1,1) = stim(1,1);
        for a = 2:size(stim,2)
            if stim(1,a) > refPer
                stimAll(1,size(stimAll,2)+1) = stim(1,a);
                refPer = stim(1,a)+sRef;
            end 
        end
        
        sITI = ceil(EEG.srate./1000.*options.ITI); %converts ITI to samples
        prec = ceil(EEG.srate./1000.*5); %precision for searching for second pulse = +/-5 ms
        diffStim = diff([stimAll(1,1)-options.ITI stimAll]);
        stimLabel = [];
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
        
        n=size(EEG.event,2);
        p = size(EEG.urevent,2);
        for a=1:size(stimAll,2);
            EEG.event(1,a+n).type=stimLabel{1,a};
            EEG.event(1,a+n).latency=epTime(x,stimAll(1,a));
            EEG.event(1,a+p).urevent=p+a;
            EEG.urevent(1,a+n).type=stimLabel{1,a};
            EEG.urevent(1,a+p).latency=epTime(x,stimAll(1,a));
        end;
        
    end
    
    %Re-epochs data
    EEG = pop_epoch( EEG, {stimLabel{1,1}}, newEpoch, 'epochinfo', 'yes');
end

%Display
fprintf('New triggers added and data re-epoched.\n');

end

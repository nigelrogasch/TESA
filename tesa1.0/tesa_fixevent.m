% tesa_fixevent() - finds TMS pulses by detecting the large TMS artifacts
%                   present in already epoched data. This script is
%                   designed for instances when the recorded events do
%                   not correspond with when the TMS pulse was given.
%                   The script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses can also be detected.
% 
%                   IMPORTANT: If you need to use this function, make sure
%                   that the initial epochs you use are larger than the
%                   final epoch size you desire. If the initial epoch size is
%                   too small, the new epoch window will be out of range 
%                   with the new events. E.g. inital epoch is -100 to 100
%                   and the event is shifted 10 ms so the new 0 now sits
%                   at 10 ms. Re-epoching the data to -100 to 100 won't
%                   work as the new range is effectively -90 to 110. In this case, 
%                   run the initial epoch at -120 to 120 and the epoch can now
%                   be taken.
%                   
%
% Usage:
%   >>  EEG = tesa_fixevent( EEG, elec, newEpoch, tmsLabel );
%   >>  EEG = tesa_fixevent( EEG, elec, newEpoch, tmsLabel, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - [required] string with electrode to use for finding artifact
%   newEpoch        - [required] vector with start and end time of new epoch in
%                   seconds (following pop_epoch convention). 
%                   Example: [-1,1] %For -1 s to 1s epoch
%   tmsLabel        - [required] string indicating the event that requires 
%                   correcting (e.g. 'TMS')  
%                   
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artefact to recover below the rate of change). int
%                   is in ms. 
%                   default = 3
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e4
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - [required] 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
%   'ISI', [int]    - [required] [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [1,2,...]. 
%                   default = []
%     
% Outputs:
%   EEG             - EEGLAB EEG structure
%
%   % Examples
%   EEG = tesa_fixevent( EEG, 'Cz', [-0.8,0.8], 'TMS' ); %default use
%   EEG = tesa_fixevent( EEG, 'Fz', [-0.7,0.7], 'TMS', 'refract', 4, 'rate', 2e5 ); %user defined input
%   EEG = tesa_fixevent( EEG, 'Cz', [-0.8,0.8], 'LICI', 'paired', 'yes', 'ISI', 100 ); %paired pulse use
%
% See also:
%   tesa_findpulse, tesa_findpulsepeak

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

function EEG = tesa_fixevent( EEG, elec, newEpoch, tmsLabel, varargin )

if nargin < 4
	error('Not enough input arguments.');
end

%define defaults
options = struct('refract',3,'rate',1e4,'paired','no','ISI',[]);

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

%check new epoch inputs
if size(newEpoch,2) ~= 2
    error('Input for ''newEpoch'' needs to be in the format [start, end]. e.g. [-1,1]. Note time is in seconds.');
elseif newEpoch(1,1) < EEG.xmin || newEpoch(1,2) > EEG.xmax
    error('The new epoch size is outside of the range of the existing epoch. Note time is in seconds.');
end

%check that tmsLabel exists
for a = 1:size(EEG.event,2)
    eventName{a,1} = EEG.event(a).type;
end
eventNameU = unique(eventName);

if sum(strcmp(tmsLabel,eventNameU)) == 0
   eventNameU = eventNameU';
    for a = 1:size(eventNameU,2)
        eventNameU{2,a} = ' ';
    end 
    tempOut = reshape(eventNameU,1,[]);
    tempStr = [tempOut{1,:}];
    error('''tmsLabel'' input ''%s'' does not exist. Available labels are: %s.\n',tmsLabel,tempStr)
end
    
%check that paired and repetitive have been correctly called
if ~(strcmp(options.paired,'no') || strcmp(options.paired,'yes'))
    error('paired must be either ''yes'' or ''no''.');
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

%For single pulse data
if strcmp(options.paired,'no')
    for a = 1:size(data,3)

        %Calculate derivatives
        h = (1/EEG.srate)*1000;     %step size in ms
        der1 = diff(data(1,:,a))/h; %calculates first derivative

        %finds artifact (defined as first derivative > rate)
        rateS = options.rate.*h; %Convert rate in to change in uV per sample
        logstim = abs(der1)>rateS;
        samp =(1:size(data,2));
        stim = samp(logstim);
        
        if isempty(stim)
            error('No stimulus artifacts were found in trial %d. Please adjust ''rate'' trying a lower number. Note that ''rate'' default is 1e4.',a);
        end
            
        %Creates new event info
        for b = 1:size(EEG.event,2)
            if strcmp(EEG.event(b).type,tmsLabel) && EEG.event(b).epoch == a
                EEG.event(b).latency = epTime(a,stim(1,1));
            end
        end
    end

    %Re-epochs data
    EEG = pop_epoch( EEG, {tmsLabel}, newEpoch, 'epochinfo', 'yes');
end
    
%For paired pulse data
if strcmp(options.paired,'yes')
            
        %Check that ISI has been provided
        if isempty(options.ISI)
            error('Please provide the interstimulus interval (ISI) for detecting paired pulse.');
        end

        %Check that refractory period is less that the ISI
        if options.refract > options.ISI
            error('The refractory period ''%d'' is shorter than the interstimulus interval ''%d''. This will result in inaccurate detection of the test pulse. Please shorten refractory period.\n',options.refract,options.ISI);
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

            %Remove events within refractory period
            sRef = ceil(EEG.srate./1000.*options.refract); %converts refractory period to samples
            refPer = stim(1,1)+sRef; %defines refractory period following stimulus
            stimAll = [];
            stimAll(1,1) = stim(1,1);
            for a = 2:size(stim,2)
                if stim(1,a) > refPer
                    stimAll(1,size(stimAll,2)+1) = stim(1,a);
                    refPer = stim(1,a)+sRef;
                end 
            end

            % Check if artifact was found
            if isempty(stimAll)
                error('No stimulus artifacts were found in trial %d. Please adjust ''rate'' trying a lower number. Note that ''rate'' default is 1e4.',x);
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
                        stimLabel{1,2} = tmsLabel;
                    else
                        stimLabel{1,1} = tmsLabel;
                    end
                end
            end
            
            %Creates new event info
            for b = 1:size(EEG.event,2)
                if strcmp(EEG.event(b).type,tmsLabel) && EEG.event(b).epoch == x
                    if size(stimLabel,2) == 1
                        error('Only one artefact was found in trial %d. Please check the ''ISI'' or adjust ''rate'' trying a lower number. Note that ''rate'' default is 1e4.',x);
                    end
                    EEG.event(b).latency = epTime(x,stimAll(1,2));
                end
            end
            
        end
        
        %Re-epochs data
        EEG = pop_epoch( EEG, {tmsLabel}, newEpoch, 'epochinfo', 'yes');
end

%Display
fprintf('Latency of events ''%s'' adjusted and data re-epoched.\n',tmsLabel);

end

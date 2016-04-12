% tesa_edm()    - find artefactual components automatically by using the 
%                 enhanced deflation method (EDM) algorithm as advocated in 
%                 the following paper:
%
%                 Korhonen, Hernandez-Pavon et al (2011) Removal of 
%                 large muscle artifacts from transcranial magnetic
%                 stimulation-evoked EEG by independent component
%                 analysis. Med Biol Eng Compt, 49:397-407.
% 
%                 The function uses the same algorithmn as tesa_compselect
%                 to automatically detect components representing the
%                 TMS-evoked muscle artifact. Settings for this detection
%                 can be manually altered by the user.
% Usage:
%   >>  EEG = tesa_edm( EEG );
%   >>  EEG = tesa_edm( EEG, chanlocs, Nic, sf );
%   >>  EEG = tesa_edm( EEG, chanlocs, Nic, sf,  'key1', value1...);
%
% Inputs:
%   EEG                - EEGLAB EEG structure
%   chanlocs           - Channel locations 
%   Nic                - [Integer] Number of independent components to look for
%                        Default = rank(EEG)-5 (to make sure the
%                        algorithm converges).
%   sf                 - Sampling Frequency in Hz
% 
% Optional input pairs (varargin):
%   'comps',int         - int is an integer describing the number of components
%                       to perform selection on (e.g. first 10 components).
%                       Leave empty for all components. 
%                       Default: []
%
% Optional input pairs for detecting artifact components (varargin):
%   
%   TMS-evoked muscle activity
%   This type of artifact is detected by comparing the mean absolute amplitude
%   of the component time course within a target window ('tmsMuscleWin') 
%   and the mean absolute amplitude across the entire component time course. A threshold is
%   set by the user for detection (e.g. 8 means the mean absolute amplitude in the target window
%   is 8 times larger than the mean absolute amplitude across the entire time course). 
%   Components are stored under EEG.edmBadComp.
%   'tmsMuscle','str'   - 'on' | 'off'. A string which turns on TMS-evoked muscle 
%                       activity detection. 
%                       Default: 'on'
%   'tmsMuscleThresh',int - An integer determining the threshold for
%                       detecting components representing TMS-evoked muscle activity.
%                       Default: 8
%   'tmsMuscleWin',[start,end] - a vector describing the target window for
%                       TMS-evoked muscle activity (in ms).
%                       Default: [11,50]
%   'tmsMuscleFeedback','str' - 'on' | 'off'. String turning on feedback
%                       of TMS-evoked muscle threshold value for each component
%                       in the command window.
%                       (Useful for determining a suitable threshold).
%                       Default: 'off'
% 
% Outputs:
%   EEG                - EEGLAB EEG structure, data after removing
%                        artefactual ICs
%
% Examples:
%  EEG = tesa_edm( EEG ); %default use
%  EEG = tesa_edm( EEG, [], 30, 1000); % only look for 30 components, sampling rate is 1000 Hz
%  EEG = tesa_edm( EEG, [], [], [], 'comps', 10, 'tmsMuscleThresh',10,'tmsMuscleWin',[11,50],'tmsMuscleFeedback','on'); % only plot the top 10 components, change the threshold for artefact detection to 10, change the window for comparions to 11-50 ms and return threshold values for each component in the command window.
% 
% See also:
%   tesa_pcacompress, tesa_fastica, tesa_edmpreprocess, tesa_edmcompselect, icaedm

% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com; Nigel Rogasch, Monash University,
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


function EEG = tesa_edm(EEG, chanlocs, Nic, sf, varargin)

%define defaults
options = struct('comps',[], 'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11,50],'tmsMuscleFeedback','off');

% read the acceptable names
optionNames = fieldnames(options);

% if nargin < 4
%     error('Not enough input arguments.');
% end

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

if strcmpi(options.tmsMuscle,'on')

    if options.tmsMuscleThresh < 0
        error('Input for ''tmsMuscleThresh'' must be greater than 0.');
    elseif size(options.tmsMuscleWin,2)~=2
        error('Input for ''tmsMusclesWin'' must be in the following format: [start,end]. e.g. [11,50].');
    elseif options.tmsMuscleWin(1,1) > options.tmsMuscleWin(1,2)
        error('Input for ''tmsMusclesWin'' must be in the following format: [start,end]. e.g. [11,50].');
    elseif options.tmsMuscleWin(1,1) < EEG.times(1,1) || options.tmsMuscleWin(1,2) > EEG.times(1,end)
        error('Input for ''tmsMuscleWin'' must be within the limits of the data [%d to %d].',EEG.times(1,1),EEG.times(1,end));
    end
    
end

%Check for channel locations
if isempty(EEG.chanlocs(1).theta)
    error('Channel locations are required for this function. Please load channel locations: Edit -> Channel locations.')
end

% Preprocessing steps before running ICA 
[EEGwhite,inMat,Meanb] = tesa_edmpreprocess(EEG); 
T=size(inMat,2);

if nargin <2;
	chanlocs = EEG.chanlocs;
end
if isempty(chanlocs)
    chanlocs = EEG.chanlocs;
end

if nargin <3;
	Nic = rank(EEGwhite)-5;
end
if isempty(Nic)
    Nic = rank(EEGwhite)-5;
end

if nargin <4;
	sf = EEG.srate;
end
if isempty(sf)
    sf = EEG.srate;
end

%Checks that number of ICs is smaller than the rank of the data matrix
rankMat = rank(EEGwhite);
if Nic > rankMat
    error('Number of ICs for consideration (%d) is larger than than the rank of the data matrix (%d). Please revise.',Nic,rankMat);
end

if ~isempty(options.comps)
    if options.comps > rankMat 
        error('Number of ICs to plot (%d) is larger than than the rank of the data matrix (%d). Please revise.',options.comps,rankMat);
    elseif options.comps > Nic
        error('Number of ICs to plot (%d) is larger than than number of ICs for consideration (%d). Please revise.',options.comps,Nic);
    end
end

%%  Finding IC's by EDM

EEGwhite=double(EEGwhite);  

fprintf('Performing ICA using enhanced deflation method (EDM).\n');
[Wbest]=tesa_icaedm(EEGwhite,Nic); %Function to compute EDM
 
% Computing ICA time-courses and topographies
Sica=Wbest'*EEGwhite; %ICA time-courses
A=1/T*inMat*Sica'; % ICA topographies

% Reshapes 2D matrix to 3D matrix
S = reshape(Sica,size(Sica,1),size(EEG.data,2),size(EEG.data,3));
Sout=mean(S,3);

global redo
redo = true;

while redo == true;
    pre = mean(EEG.data,3);
    EEG1 = EEG;

    badComp = tesa_edmcompselect(EEG, S, A, 'comps',options.comps, 'tmsMuscle',options.tmsMuscle,'tmsMuscleThresh',options.tmsMuscleThresh,'tmsMuscleWin',options.tmsMuscleWin, 'tmsMuscleFeedback',options.tmsMuscleFeedback);

    % Correcting the data by removing artifactual components
    data = inMat-A(:,badComp)*Sica(badComp,:); %Corrected data
    data=data+Meanb; % Adding the mean back
    post = reshape(data,size(EEG.data,1),size(Sout,2),[]); % Reshape matrix in to eeglab format
    
    %Plot check
    f1 = figure('KeyPressFcn',@(obj,evt) 0);
    f1.Name = 'Press enter when ready to continue.';
    f1.NumberTitle = 'off';
    
    plotTimeX = [-100,500];
    
    subplot(1,2,1)
    plot(EEG.times,pre,'b'); grid on; hold on;
    plot([0 0], get(gca,'ylim'),'r--');
    set(gca,'Xlim', plotTimeX);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('Before correction','fontweight','bold');
    yScale = get(gca,'ylim');

    subplot(1,2,2)
    plot(EEG.times,mean(post,3),'b'); grid on; hold on;
    plot([0 0], yScale,'r--');
    set(gca,'Xlim', plotTimeX,'ylim',yScale);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('After correction','fontweight','bold');
    
    waitfor(gcf,'CurrentCharacter');
    curChar=uint8(get(gcf,'CurrentCharacter'));

    %Check if happy to continue
    ButtonName = questdlg('Are you satisfied with component removal?', 'Verify component removal', 'No-Redo', 'Yes', 'Cancel', 'Yes');
    switch ButtonName,
        case 'No-Redo',
            redo = true; EEG = EEG1;
        case 'Yes',
            redo = false; 
        case 'Cancel',
            redo = false; error('Script termninated by user.');
    end 

    close;

end

EEG.data = reshape(data,size(EEG.data,1),size(Sout,2),[]); % Reshape matrix in to eeglab format
EEG.edmBadComp = badComp;

fprintf('EDM correction complete.\n');

end




% pop_tesa_edm()    - find artefactual components automatically by using the 
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
%   >>  EEG = pop_tesa_edm( EEG ); % pop up window
%   >>  EEG = pop_tesa_edm( EEG, chanlocs, Nic, sf );
%   >>  EEG = pop_tesa_edm( EEG, chanlocs, Nic, sf,  'key1', value1...);
%
% Inputs:
%   EEG                - EEGLAB EEG structure
%   chanlocs           - Channel locations - taken from EEG.chanlocs
%   Nic                - [Integer] Number of independent components to look for
%                        Default = rank(EEG)-5 (to make sure the
%                        algorithm converges).
%   sf                 - Sampling Frequency in Hz - taken from EEG.srate
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
%                       Default: [11,30]
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
%  EEG = pop_tesa_edm( EEG, [], 30, 1000); % only look for 30 components, sampling rate is 1000 Hz
%  EEG = pop_tesa_edm( EEG, [], [], [], 'comps', 10, 'tmsMuscleThresh',10,'tmsMuscleWin',[11,50],'tmsMuscleFeedback','on'); % only plot the top 10 components, change the threshold for artefact detection to 10, change the window for comparions to 11-50 ms and return threshold values for each component in the command window.
% 
% See also:
%   pop_tesa_pcacompress, pop_tesa_fastica, tesa_edmpreprocess, tesa_edmcompselect, icaedm

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

 function [EEG com] = pop_tesa_edm(EEG, chanlocs, Nic, sf, varargin);
 
 
com = '';          

%check that data is present
if isempty(EEG.data);
    error('Data is empty');
end

rankMat=rank(reshape(EEG.data,size(EEG.data,1),[],1));
% pop up window
% -------------
if nargin < 3
   
   geometry = { 1 [1 0.5] [1 0.5] 1 [1 0.5] [1 0.5] [1 0.5] [1 0.5] 1 1 1};
    
    uilist =  {{'style', 'text', 'string', 'The number of ICs must be <= to the rank of the data matrix','fontweight','bold'}...
              {'style', 'text', 'string', 'Number of ICs to consider (blank for all)'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', 'Number of ICs to plot (blank for all)'} ...
              {'style', 'edit', 'string', ''} ...
              {}...
              {'style', 'text', 'string', 'Auto detect TMS-evoked muscle artefact','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '8'} ...
              {'style', 'text', 'string', 'Target times for artefact (ms): start, end'} ...              
              {'style', 'edit', 'string', '11, 30'} ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {}...
              {'style', 'text', 'string', 'Press enter when you have categorised each IC.','fontAngle','italic'}...
              {'style', 'text', 'string', 'Blue plots are auto-categorised as neural, red plots as artefacts.','fontAngle','italic'}...
              };  
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Enahnced deflation method (EDM) -- pop_tesa_edm()', 'helpcom', 'pophelp(''pop_tesa_edm'')');
    if isempty(result), return; end;
    
    %IC values and Sampling frequency
    if strcmp(result{1,1},'')
        Nic = [];
    else
        Nic = str2num(result{1,1});
    end
    comps = str2num(result{1,2});
    if result{1,3} == 1
        tmsMuscle = 'on';
    else 
        tmsMuscle = 'off';
    end
    tmsMuscleThresh = str2num(result{1,4});
    tmsMuscleWin = str2num(result{1,5});
    if result{1,6} == 1
        tmsMuscleFeedback = 'on';
    else
        tmsMuscleFeedback = 'off';
    end    

end
if nargin < 2
    chanlocs = [];
end
if nargin < 3 && ~exist('Nic','var')
    Nic = [];
end
if nargin < 4
    sf = [];
end

if nargin == 1
    EEG = tesa_edm(EEG, [], Nic, [], 'comps',comps,'tmsMuscle',tmsMuscle,'tmsMuscleThresh',tmsMuscleThresh,'tmsMuscleWin',tmsMuscleWin,'tmsMuscleFeedback',tmsMuscleFeedback);
    com = sprintf('%s = pop_tesa_edm( %s, [], %s, [], ''comps'', %s, ''tmsMuscle'', ''%s'', ''tmsMuscleThresh'', %s, ''tmsMuscleWin'', %s, ''tmsMuscleFeedback'', ''%s'' );', inputname(1), inputname(1), mat2str(Nic), mat2str(comps), tmsMuscle, mat2str(tmsMuscleThresh), mat2str(tmsMuscleWin), tmsMuscleFeedback );
elseif nargin > 1 && nargin < 5
    EEG = tesa_edm(EEG, chanlocs, Nic, sf);
    com = sprintf('%s = pop_tesa_edm( %s, %s, %s, %s );', inputname(1), inputname(1), mat2str(chanlocs),mat2str(Nic),mat2str(sf));
elseif nargin > 4
    EEG = tesa_edm(EEG, [], Nic, [], varargin{:});
    com = sprintf('%s = pop_tesa_edm( %s, [], %s, [], %s );', inputname(1), inputname(1), mat2str(Nic),vararg2str(varargin));
end

% pop_tesa_compselect() - sorts components following ICA by time course
%                   variance, applies a series of rules (determined by the
%                   user) to classify components as artifacts or neural
%                   activity, presents the components one-by-one for manual
%                   clarification and then removes components selected as
%                   artifacts from the data. During the clarification step,
%                   users can re-classify components using a drop down
%                   menu. When satisfied, press enter (or any other button) to 
%                   continue to the next component. Components are classified as one of the
%                   following: neural, TMS-evoked muscle, eye, muscle,
%                   electrode noise, sensory or other. If a component is selected
%                   as an artifact using a rule, it is plotted as red, if
%                   neural it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaBadComp.
% 
%                   Users are encouraged to report the settings used with
%                   this function in any publications.
% 
%                   Note that the component selection rules are intended
%                   to work as a guide only, and each component must be
%                   manually checked and re-classified if necessary. The rule 
%                   thresholds can be adapted by the user and may differ 
%                   between datasets. There is no guarantee for the success
%                   of these rules in accurately classifying components.
%                   For further information on selecting components, please
%                   read:
%                   
%                   Rogasch NC et al. (2014) Removing artefacts from TMS-EEG 
%                   recordings using independent component analysis: Importance 
%                   for assessing prefrontal and motor cortex network
%                   properties. NeuroImage, 101:425-439               
%
% Usage:
%   >>  EEG = pop_tesa_compselect( EEG ); % pop up window
%   >>  EEG = pop_tesa_compselect( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'compCheck','str' - 'on' | 'off'. Turns on or off the component plots
%                   for checking automated component selection. On is highly recommended. 
%                   Default: 'on'
%   'comps',int     - int is an integer describing the number of components
%                   to perform selection on (e.g. first 10 components).
%                   Leave empty for all components. 
%                   Default: []
%   'figSize','str' - 'small' | 'medium' | 'large'; Determines the size of
%                   the figures that are plotted displaying the information
%                   on the components.
%                   Default: 'small'
%   'plotTimeX',[start,end] - vector with integers for plotting the
%                   component time course (in ms).
%                   Default: [-200,500]
%   'plotFreqX',[low,high] - vector with integers for plotting the
%                   frequency distribution of components (in Hz)
%                   Default: [1,100]
% 
% Optional input pairs for detecting artifact components (varargin):
%   
%   TMS-evoked muscle activity
%   This type of artifact is detected by comparing the mean absolute amplitude
%   of the component time course within a target window ('tmsMuscleWin') 
%   and the mean absolute amplitude across the entire component time course. A threshold is
%   set by the user for detection (e.g. 8 means the mean absolute amplitude in the target window
%   is 8 times larger than the mean absolute amplitude across the entire time course). 
%   Components are stored under tmsMuscle.
%   'tmsMuscle','str' - 'on' | 'off'. A string which turns on TMS-evoked muscle 
%                   activity detection. 
%                   Default: 'on'
%   'tmsMuscleThresh',int - An integer determining the threshold for
%                   detecting components representing TMS-evoked muscle activity.
%                   Default: 8
%   'tmsMuscleWin',[start,end] - a vector describing the target window for
%                   TMS-evoked muscle activity (in ms).
%                   Default: [11,30]
%   'tmsMuscleFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of TMS-evoked muscle threshold value for each component
%                   in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Eye blinks
%   This type of artifact is detected by comparing the mean absolute z score
%   (calculated on the component topography weights) of two electrodes close 
%   to the eyes (e.g. 'Fp1' and 'Fp2'). A threhold is set by the user for detection 
%   (e.g. 2.5 - the mean z score of the two electrodes needs to be larger than
%   2.5 for detection). Components are stored under eyes.
%   'blink','str' - 'on' | 'off'. A string which turns on blink detection. 
%                   Default: 'on'
%   'blinkThresh',int - An integer determining the threshold for
%                   detecting components representing eye blinks.
%                   Default: 2.5
%   'blinkElecs',{'str','str'} - a cell with the two electrodes for blink
%                   detection. The two electrodes closest to the eyes are
%                   best for this purpose.
%                   Default: {'Fp1','Fp2'}
%   'blinkFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of blink detection threshold value for each component in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Lateral eye movement
%   This type of artifact is detected by comparing the z scores
%   (calculated on the component topography weights) of two electrodes on either side
%   of the forehead (e.g. 'F7','F8'). The z score must be positive for one and 
%   negative for the other. A threhold is set by the user for detection 
%   (e.g. 2 means the z score in one electrode must be greater than 2 and less than
%   -2 in the other electrode). Components are stored under eyes.
%   'move','str' - 'on' | 'off'. A string which turns on movement detection. 
%                   Default: 'on'
%   'moveThresh',int - An integer determining the threshold for
%                   detecting components representing lateral eye movement.
%                   Default: 2
%   'moveElecs',{'str','str'} - a cell with the two electrodes for movement
%                   detection. The two electrodes on either side of the
%                   forehead are best for this purpose.
%                   Default: {'F7','F8'}
%   'moveFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of movement detection threshold value for each
%                   component in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Persistent muscle activity
%   This type of artifact is detected by comparing the mean power of the 
%   component time course frequency distribution between a target window
%   and the mean across all frequencies(calculated using an FFT across all trials). 
%   A threhold is set by the user for detection (e.g. 0.6 means the high frequency 
%   power is 60% of the total power). 
%   Components are stored under muscle.
%   'muscle','str' - 'on' | 'off'. A string which turns on muscle detection. 
%                   Default: 'on'
%   'muscleThresh',int - An integer determining the threshold for
%                   detecting components representing muscle activity.
%                   Default: 0.6
%   'muscleFreqWin',[int,int] - a vector with the frequencies for the target w
%                   window(in Hz).
%                   Default: [31,100] - Note the default higher value will
%                   adjust to match the maximum plot frequency defined by
%                   'plotFreqX'.
%   'muslceFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of muscle detection threshold value for each component
%                   in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Electrode noise
%   This type of artifact is detected by comparing z scores in individual 
%   electrodes (calculated on the component topography weights). A threhold is 
%   set by the user for detection (e.g. 4 means one or more electrodes has
%   an absolute z score of at least 4). Components are stored in elecNoise.
%   'elecNoise','str' - 'on' | 'off'. A string which turns on electrode noise detection. 
%                   Default: 'on'
%   'elecNoiseThresh',int - An integer determining the threshold for
%                   detecting components representing electrode noise.
%                   Default: 4
%   'elecNoiseFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of electrode noise detection threshold value for each component
%                   in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Note that a sensory category is also available in the dropdown menu.
%   There are currently no specific rules suitable for accurately detecting
%   these components and are therefore classified entirely at the user's discretion.
% 
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_compselect( EEG ); % default use
%   EEG = pop_tesa_compselect( EEG, 'comps', 10, 'figSize', 'medium','plotTimeX',[-600,600], 'plotFreqX', [2,45] ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots.
%   EEG = pop_tesa_compselect( EEG, 'elecNoise','off','blinkThresh',3,'blinkElecs',{'AF3','AF4'},'blinkFeedback','on'); % turn off electrode noise detection, change threshold for blinks to 3, change electrodes used to AF3 and AF4 and turn on the feedback of blink threhsolds for  individual components in the command window.
% 
% See also:
%   pop_tesa_fastica, tesa_sortcomps 

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

 function [EEG com] = pop_tesa_compselect( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
            
    geometry = {[0.5 1] 1 [1 0.5 1 0.5] [1 0.5 1 0.5] 1 ... %Figure
                [0.45 1] 1 [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1 ... %TMS-evoked muscle and Eye Blink
                [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1  ... %Lateral eye movement and Muscle
                [1 0.5 1 0.5] [1 0.5 1.5] [1 0.5 1.5] 1 1 ... %Electrode noise
                };
            
    uilist = {{'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Independent component selection','fontweight','bold'} ... 
              {'style', 'text', 'string', 'Figure settings','fontweight','bold'} ...  
              {'style', 'text', 'string', 'Components to plot (blank for all)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Figure size'} ...
              {'style', 'popupmenu', 'string', 'small|medium|large', 'tag', 'interp' } ...
              {'style', 'text', 'string', 'Time course limits (ms): start, end'} ...
              {'style', 'edit', 'string', '-200, 500'} ...
              {'style', 'text', 'string', 'Frequency limits (Hz): start, end'} ...
              {'style', 'edit', 'string', '1, 100'}...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Settings for auto-detecting artifact components','fontweight','bold'} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'TMS-evoked muscle activity','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Eye blinks','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...                            
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '8'} ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '2.5'} ...
              {'style', 'text', 'string', 'Target times (ms): start, end'} ...              
              {'style', 'edit', 'string', '11, 30'} ...
              {'style', 'text', 'string', 'Electrodes for detection'} ...
              {'style', 'edit', 'string', 'Fp1, Fp2'} ...              
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Lateral eye movements','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...              
              {'style', 'text', 'string', 'Persistant muscle activity','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '2'} ...              
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '0.6'} ...
              {'style', 'text', 'string', 'Electrodes for detection'} ...
              {'style', 'edit', 'string', 'F7, F8'} ...
              {'style', 'text', 'string', 'Target frequencies (Hz): start, end'} ...
              {'style', 'edit', 'string', '30, 100'} ...              
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', ''} ...                            
              {'style', 'text', 'string', 'Electrode noise','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Check auto-detect results','fontweight','bold'} ...
              {'style', 'popupmenu', 'string', 'on|off', 'tag', 'check' } ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '4'} ...
              {'style', 'text', 'string', '  For further information on settings for artifact','fontAngle','italic'} ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', '  component selection, see help or the user manual','fontAngle','italic'} ...
              {}...
              {'style', 'text', 'string', 'Press enter when you have categorised each IC. Blue plots are auto-categorised as neural, red plots as artifacts.','fontAngle','italic'}...
              };
          
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find artifact components -- pop_tesa_compselect()', 'helpcom', 'pophelp(''pop_tesa_compselect'')');
    if isempty(result), return; end;
    
    %Extract data
    comps = str2num(result{1,1});
    if result{1,2} == 1
        figSize = 'small';
    elseif result{1,2} == 2
       	figSize = 'medium';
    elseif result{1,2} == 3
       	figSize = 'large';
    end
    plotTimeX = str2num(result{1,3});
    plotFreqX = str2num(result{1,4});
    
    if result{1,5} == 1
        tmsMuscle = 'on';
    else
        tmsMuscle = 'off';
    end
    if result{1,6} == 1
        blink = 'on';
    else
        blink = 'off';
    end
    tmsMuscleThresh = str2num(result{1,7});
    blinkThresh = str2num(result{1,8});
    tmsMuscleWin = str2num(result{1,9});
    blinkElecs = strtrim(strsplit(strrep(result{1,10},',','')));
    blinkString = ['''',result{1,10},''''];
    blinkString = strrep(blinkString,',',''',''');
    blinkString = strrep(blinkString,' ','');
    if result{1,11} == 1
        tmsMuscleFeedback = 'on';
    else
        tmsMuscleFeedback = 'off';
    end
    if result{1,12} == 1
        blinkFeedback = 'on';
    else
        blinkFeedback = 'off';
    end
    
    if result{1,13} == 1
        move = 'on';
    else
        move = 'off';
    end
    if result{1,14} == 1
        muscle = 'on';
    else
        muscle = 'off';
    end
    moveThresh = str2num(result{1,15});
    muscleThresh = str2num(result{1,16});
    moveElecs = strtrim(strsplit(strrep(result{1,17},',','')));
    moveString = ['''',result{1,17},''''];
    moveString = strrep(moveString,',',''',''');
    moveString = strrep(moveString,' ','');
    muscleFreqWin = str2num(result{1,18});
    if result{1,19} == 1
        moveFeedback = 'on';
    else
        moveFeedback = 'off';
    end
    if result{1,20} == 1
        muscleFeedback = 'on';
    else
        muscleFeedback = 'off';
    end
    
    if result{1,21} == 1
        elecNoise = 'on';
    else
        elecNoise = 'off';
    end
    if result{1,22} == 1
        compCheck = 'on';
    else
        compCheck = 'off';
    end
    elecNoiseThresh = str2num(result{1,23});
    if result{1,24} == 1
        elecNoiseFeedback = 'on';
    else
        elecNoiseFeedback = 'off';
    end
    
end

%Run script from input
if nargin <2
    EEG = tesa_compselect(EEG,'compCheck',compCheck,'comps',comps,'figSize',figSize,'plotTimeX',plotTimeX,'plotFreqX',plotFreqX,'tmsMuscle',tmsMuscle,'tmsMuscleThresh',tmsMuscleThresh,'tmsMuscleWin',tmsMuscleWin,'tmsMuscleFeedback',tmsMuscleFeedback,'blink',blink,'blinkThresh',blinkThresh,'blinkElecs',blinkElecs,'blinkFeedback',blinkFeedback,'move',move,'moveThresh',moveThresh,'moveElecs',moveElecs,'moveFeedback',moveFeedback,'muscle',muscle,'muscleThresh',muscleThresh,'muscleFreqWin',muscleFreqWin,'muscleFeedback',muscleFeedback,'elecNoise',elecNoise,'elecNoiseThresh',elecNoiseThresh,'elecNoiseFeedback',elecNoiseFeedback);
    com = sprintf('%s = pop_tesa_compselect( %s,''compCheck'',''%s'',''comps'',%s,''figSize'',''%s'',''plotTimeX'',%s,''plotFreqX'',%s,''tmsMuscle'',''%s'',''tmsMuscleThresh'',%s,''tmsMuscleWin'',%s,''tmsMuscleFeedback'',''%s'',''blink'',''%s'',''blinkThresh'',%s,''blinkElecs'',{%s},''blinkFeedback'',''%s'',''move'',''%s'',''moveThresh'',%s,''moveElecs'',{%s},''moveFeedback'',''%s'',''muscle'',''%s'',''muscleThresh'',%s,''muscleFreqWin'',%s,''muscleFeedback'',''%s'',''elecNoise'',''%s'',''elecNoiseThresh'',%s,''elecNoiseFeedback'',''%s'' );', inputname(1), inputname(1),compCheck,mat2str(comps),figSize,mat2str(plotTimeX),mat2str(plotFreqX),tmsMuscle,mat2str(tmsMuscleThresh),mat2str(tmsMuscleWin),tmsMuscleFeedback,blink,mat2str(blinkThresh),blinkString,blinkFeedback,move,mat2str(moveThresh),moveString,moveFeedback,muscle,mat2str(muscleThresh),mat2str(muscleFreqWin),muscleFeedback,elecNoise,mat2str(elecNoiseThresh),elecNoiseFeedback );
elseif nargin > 2
    EEG = tesa_compselect(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_compselect( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

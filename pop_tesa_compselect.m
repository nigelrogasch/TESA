% pop_tesa_compselect() - applies a series of rules (determined by the
%                   user) to classify components as artifacts, presents the 
%                   components one-by-one for manual
%                   clarification and then removes components selected as
%                   artifacts from the data. During the clarification step,
%                   users can re-classify components using a drop down
%                   menu. When satisfied, press the 'Next' button to 
%                   continue to the next component. Components are classified as one of the
%                   following: keep, reject, reject - TMS-evoked muscle, reject - blink, 
%                   reject - eye movement, reject - muscle, reject - 
%                   electrode noise, or reject - sensory . If a component is selected
%                   as an artifact using a rule, it is plotted as red, if
%                   not it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaCompClassTesaX.
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
%                   and
% 
%                   Rogasch NC et al. (2014) Analysing concurrent transcranial magnetic stimulation and
%                   electroencephalographic data: A review and introduction to the open-source
%                   TESA software. NeuroImage, 147:934–951                             
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
%   'remove','str'  - 'on' | 'off'. Turns on or off removal of components
%                   after classification
%                   Default: 'on'
%   'saveWeights','str' - 'on' | 'off'. Saves the winv (topoplot weights)
%                   and icaact (time course weights) that the
%                   classification was performed on. Note this will
%                   increase the storage size of the files.
%                   Default: 'off'
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
%   'freqScale','str' - 'raw' | 'log' | 'log10' | 'db'
%                   y-axis scaling for the frequency plots
%                   Default: 'log'
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
%                   in a figure.
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
%                   of blink detection threshold value for each component in a figure.
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
%                   component in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Persistent muscle activity
%   This type of artifact is detected by calculating the slope of the line 
%   of best fit between the log EEG power and the log frequency. See the 
%   following manuscript for details:
%
%   Fitzgibbon SP et al. (2016) Automatic determination of EMG-contaminated 
%   components and validation of independent component analysis using EEG during
%   pharmacologic paralysis. Clinical Neurophysiology 127:1781-1793
% 
%   Components are stored under muscle.
%   'muscle','str' - 'on' | 'off'. A string which turns on muscle detection. 
%                   Default: 'on'
%   'muscleThresh',int - An integer determining the threshold for
%                   detecting components representing muscle activity.
%                   Default: -0.31
%   'muscleFreqIn',[int,int] - a vector with the frequency limits for
%                   including in the slope analysis. Leave empty to include
%                   all of the frequency spectra. [7,70] is recommended.
%                   Default: [7,70]
%   'muscleFreqEx',[int,int] - a vector with frequencies to exclude (e.g
%                   due to line noise. Leave empty to include all of the
%                   frequency spectra.
%                   Example: [48,52] % Exclude 50 Hz line noise
%                   Default: [48,52];
%   'muslceFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of muscle detection threshold value for each component
%                   in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Note that the legacy TESA heuristic for detecting persitent muscle
%   activity can still be called from the command line. The two following
%   arguments must be called with the tesa_compselect function:
%   'muscleLegacy','str' - 'on' | 'off'. Call this key/value pair to turn
%                   on TESA legacy heuristic (e.g. comparing the mean
%                   values between a target window and the remaining
%                   window using a ratio)
%                   Default: 'off'
%   'muscleFreqWin',[int,int] - a vector with the frequencies for the target w
%                   window(in Hz).
%                   Example: [30,100]
%                   Default: []
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
%                   in a figure.
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
%   EEG = pop_tesa_compselect( EEG, 'figSize', 'medium','plotTimeX',[-600,600], 'plotFreqX', [2,45] ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots.
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

% Change log
% 13.6.2018 - removed sort components step (added directly to tesa_fastica)
% 13.6.2018 - removed plotting function and replaced with tesa_compplot
% 13.6.2018 - changed the way component classifications are stored. Now stored in EEG.icaCompClassTesaX
% 13.6.2018 - feedback on classification thresholds is now provided as a
%               figure as opposed to on the command line
% 13.6.2018 - changed the classification categories to 'keep' and 'reject' categories
% 14.6.2018 - included 'remove' input to determine whether to classify
%               components ('off'), or classify and remove components ('on')
% 10.4.2019 - Added an option to choose y-axis frequency scale
% 11.4.2019 - Changed the heuristic rule for detecting persistent muscle
%               activity to the slope of the line of best fit on the
%               log-log of the frequency spectra data (see Fitzgibbon et al
%               2016 Clin Neurophysiology). Note that the legacy TESA
%               muscle heuristic can still be called from the command line
% 11.4.2019 - Included a call to save the IC weights and time series with
%               the classification selections 

 function [EEG com] = pop_tesa_compselect( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
            
    remComps = {'on','off'};
    
    % Frequency scales
    freqScaleNames = {'raw','log','log10','db'};
    
    geometry = {[0.5 1] 1 [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1 ... %Figure
                [0.45 1] 1 [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1 ... %TMS-evoked muscle and Eye Blink
                [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1  ... %Lateral eye movement and Muscle
                [1 0.5 1 0.5] [1 0.5 1 0.5] [1 0.5 1 0.5] 1 1 1 1 ... %Electrode noise
                };
            
    uilist = {{'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Independent component selection','fontweight','bold'} ... 
              {'style', 'text', 'string', 'Figure settings','fontweight','bold'} ...  
              {'style', 'text', 'string', 'Remove components'} ...
              {'style', 'popupmenu', 'string', remComps, 'tag', 'interp'} ...
              {'style', 'text', 'string', 'Figure size'} ...
              {'style', 'popupmenu', 'string', 'small|medium|large', 'tag', 'interp', 'value', 2 } ...
              {'style', 'text', 'string', 'Time course limits (ms): start, end'} ...
              {'style', 'edit', 'string', '-200, 500'} ...
              {'style', 'text', 'string', 'Frequency limits (Hz): start, end'} ...
              {'style', 'edit', 'string', '1, 100'}...
              {'style', 'text', 'string', 'Frequency y-axis scale'} ...
              {'style', 'popupmenu', 'string', freqScaleNames, 'tag', 'freqNames','Value',2} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ''} ...
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
              {'style', 'text', 'string', 'Persistent muscle activity','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '2'} ...              
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '-0.31'} ...
              {'style', 'text', 'string', 'Electrodes for detection'} ...
              {'style', 'edit', 'string', 'F7, F8'} ...
              {'style', 'text', 'string', 'Frequencies to include (Hz): start, end'} ...
              {'style', 'edit', 'string', '7, 70'} ...              
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Frequencies to exclude (Hz): start, end'} ...
              {'style', 'edit', 'string', '48, 52'} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', ''} ...                            
              {'style', 'text', 'string', 'Electrode noise','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Check auto-detect results','fontweight','bold'} ...
              {'style', 'popupmenu', 'string', 'on|off', 'tag', 'check' } ...
              {'style', 'text', 'string', 'Threshold for detection'} ...
              {'style', 'edit', 'string', '4'} ...
              {'style', 'text', 'string', 'Save IC weights?','fontweight','bold'} ...
              {'style', 'popupmenu', 'string', 'on|off', 'tag', 'weights', 'value', 2 } ...
              {'style', 'text', 'string', 'Threshold feedback'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', ' For further information on settings for artifact component selection, see help or the user manual','fontAngle','italic'} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Press next when you have categorised each IC. Blue plots are auto-categorised as keep, red plots as reject.','fontAngle','italic'}...
              };
          
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find artifact components -- pop_tesa_compselect()', 'helpcom', 'pophelp(''pop_tesa_compselect'')');
    if isempty(result), return; end;
    
    %Extract data
    removeComps = remComps{result{1,1}};
    if result{1,2} == 1
        figSize = 'small';
    elseif result{1,2} == 2
       	figSize = 'medium';
    elseif result{1,2} == 3
       	figSize = 'large';
    end
    plotTimeX = str2num(result{1,3});
    plotFreqX = str2num(result{1,4});
    freqScale = freqScaleNames{result{1,5}};
    
    if result{1,6} == 1
        tmsMuscle = 'on';
    else
        tmsMuscle = 'off';
    end
    if result{1,7} == 1
        blink = 'on';
    else
        blink = 'off';
    end
    tmsMuscleThresh = str2num(result{1,8});
    blinkThresh = str2num(result{1,9});
    tmsMuscleWin = str2num(result{1,10});
    blinkElecs = strtrim(strsplit(strrep(result{1,11},',','')));
    blinkString = ['''',result{1,11},''''];
    blinkString = strrep(blinkString,',',''',''');
    blinkString = strrep(blinkString,' ','');
    if result{1,12} == 1
        tmsMuscleFeedback = 'on';
    else
        tmsMuscleFeedback = 'off';
    end
    if result{1,13} == 1
        blinkFeedback = 'on';
    else
        blinkFeedback = 'off';
    end
    
    if result{1,14} == 1
        move = 'on';
    else
        move = 'off';
    end
    if result{1,15} == 1
        muscle = 'on';
    else
        muscle = 'off';
    end
    moveThresh = str2num(result{1,16});
    muscleThresh = str2num(result{1,17});
    moveElecs = strtrim(strsplit(strrep(result{1,18},',','')));
    moveString = ['''',result{1,18},''''];
    moveString = strrep(moveString,',',''',''');
    moveString = strrep(moveString,' ','');
    muscleFreqIn = str2num(result{1,19});
    if result{1,20} == 1
        moveFeedback = 'on';
    else
        moveFeedback = 'off';
    end
    muscleFreqEx = str2num(result{1,21});
    if result{1,22} == 1
        muscleFeedback = 'on';
    else
        muscleFeedback = 'off';
    end
    
    if result{1,23} == 1
        elecNoise = 'on';
    else
        elecNoise = 'off';
    end
    if result{1,24} == 1
        compCheck = 'on';
    else
        compCheck = 'off';
    end
    elecNoiseThresh = str2num(result{1,25});
    if result{1,26} == 1
        saveWeights = 'on';
    else
        saveWeights = 'off';
    end
    if result{1,27} == 1
        elecNoiseFeedback = 'on';
    else
        elecNoiseFeedback = 'off';
    end
    
end

%Run script from input
if nargin <2
    EEG = tesa_compselect(EEG,'compCheck',compCheck,'remove',removeComps,'saveWeights',saveWeights,'figSize',figSize,'plotTimeX',plotTimeX,'plotFreqX',plotFreqX,'freqScale',freqScale,'tmsMuscle',tmsMuscle,'tmsMuscleThresh',tmsMuscleThresh,'tmsMuscleWin',tmsMuscleWin,'tmsMuscleFeedback',tmsMuscleFeedback,'blink',blink,'blinkThresh',blinkThresh,'blinkElecs',blinkElecs,'blinkFeedback',blinkFeedback,'move',move,'moveThresh',moveThresh,'moveElecs',moveElecs,'moveFeedback',moveFeedback,'muscle',muscle,'muscleThresh',muscleThresh,'muscleFreqIn',muscleFreqIn,'muscleFreqEx',muscleFreqEx,'muscleFeedback',muscleFeedback,'elecNoise',elecNoise,'elecNoiseThresh',elecNoiseThresh,'elecNoiseFeedback',elecNoiseFeedback);
    com = sprintf('%s = pop_tesa_compselect( %s,''compCheck'',''%s'',''remove'',''%s'',''saveWeights'',''%s'',''figSize'',''%s'',''plotTimeX'',%s,''plotFreqX'',%s,''freqScale'',''%s'',''tmsMuscle'',''%s'',''tmsMuscleThresh'',%s,''tmsMuscleWin'',%s,''tmsMuscleFeedback'',''%s'',''blink'',''%s'',''blinkThresh'',%s,''blinkElecs'',{%s},''blinkFeedback'',''%s'',''move'',''%s'',''moveThresh'',%s,''moveElecs'',{%s},''moveFeedback'',''%s'',''muscle'',''%s'',''muscleThresh'',%s,''muscleFreqIn'',%s,''muscleFreqEx'',%s,''muscleFeedback'',''%s'',''elecNoise'',''%s'',''elecNoiseThresh'',%s,''elecNoiseFeedback'',''%s'' );', inputname(1), inputname(1),compCheck,removeComps,saveWeights,figSize,mat2str(plotTimeX),mat2str(plotFreqX),freqScale,tmsMuscle,mat2str(tmsMuscleThresh),mat2str(tmsMuscleWin),tmsMuscleFeedback,blink,mat2str(blinkThresh),blinkString,blinkFeedback,move,mat2str(moveThresh),moveString,moveFeedback,muscle,mat2str(muscleThresh),mat2str(muscleFreqIn),mat2str(muscleFreqEx),muscleFeedback,elecNoise,mat2str(elecNoiseThresh),elecNoiseFeedback );
elseif nargin > 1
    EEG = tesa_compselect(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_compselect( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

% pop_tesa_compplot() - Plots independent components following ICA and allows users to 
%                   classify components for removal using a dropdown menu.
%                   Each plot includes a butterfly plot of the EEG data with
%                   and without the current IC, the average time course of
%                   the component, the trial-by-trial amplitude of the
%                   component, the topography of the component, and the
%                   power spectrum of the component. Infomration on the %
%                   mean variance (of total variance) accounted for by each
%                   component is also provided. Users can interactively
%                   alter the x and y axes of the plots by selecting
%                   'Update plots'. Note that three figures sizes are
%                   available: 'small','medium' and 'large'.
%                   When satisfied with classification, press the 'Next' button to 
%                   continue to the next component. Users can scroll backwards and forward
%                   through components by selecting the 'Next IC' and then pressing
%                   'Next'. Users can also skip to the 'Finish' from this dropdown.
%                   Components are classified as one of the following: 
%                   keep, reject, reject - TMS-evoked muscle, reject - blink, 
%                   reject - eye movement, reject - muscle, reject - 
%                   electrode noise, or reject - sensory . If a component classification
%                   algorithm has been used to classify components in TESA, this
%                   can be used as an input to tesa_compplot to automatically
%                   mark components for rejection. If a component is selected
%                   for rejection, it is plotted as red on subsequent plots, if
%                   not it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaCompClass.
% 
%                   Note that tesa_compplot is automatically called by
%                   tesa_compselect if 'compCheck' is 'on'.
% 
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
%                   TESA software NeuroImage, 147:934–951             
%
%
% Usage:
%   >>  EEG = pop_tesa_compplot( EEG ); % popup window
%   >>  EEG = pop_tesa_compplot( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'compClassInput',EEGLAB structure field - location in EEGLAB structure
%                   where output from TESA components classification
%                   algorithms are stored (e.g. tesa_compselect). Note that the number
%                   of components that have been classified must match the
%                   the number of components currently stored in the data
%                   (i.e. this won't work if you have already removed the
%                   components in tesa_compselect).
%                   Example: EEG.icaCompClass.TESA1
%                   Default: empty
%   'saveWeights','str' - 'on' | 'off'. Saves the winv (topoplot weights)
%                   and icaact (time course weights) that the
%                   classification was performed on. Note this will
%                   increase the storage size of the files.
%                   Default: 'off'
%   'figSize','str' - 'small' | 'medium' | 'large'; Determines the size of
%                   the figures that are plotted displaying the information
%                   on the components.
%                   Default: 'medium'
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
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_compplot( EEG ); % popup window
%   EEG = pop_tesa_compplot( EEG, 'compClassInput', EEG.icaCompClassTesa1 ); % Loads automated component classifications performed by tesa_compselect
%   EEG = pop_tesa_compplot( EEG, 'figSize', 'small','plotTimeX',[-600,600], 'plotFreqX', [2,45], 'freqScale','db' ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots, and the scale of the y-axis on the frequency plots.
% 
% See also:
%   tesa_compplot, tesa_compselect, tesa_sortcomps 

% Copyright (C) 2018  Nigel Rogasch, Monash University,
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

 function [EEG com] = pop_tesa_compplot( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    % Look for component classifications in EEGLAB structure
    eegFields = fieldnames(EEG);
    fieldLog = strncmp('icaCompClass',eegFields,12);
    if sum(fieldLog)>0
        if ~isempty(EEG.icaCompClass)
            classIn = fieldnames(EEG.icaCompClass);
        else
            classIn = [];
        end
    else
        classIn = [];
    end
    emptyIn = {''};
    if isempty(classIn)
        inputName = emptyIn;
    else
        inputName = [emptyIn;classIn];
    end
    
    % Figure input
    figInput = {'small','medium','large'};
    
    % Frequency scales
    freqScaleNames = {'raw','log','log10','db'};
    
    geometry = {[0.2 1] [1 0.5] 1 1 1 [1 0.5] [1 0.5] [1 0.5] [1 0.5] [1 0.5]... %Figure
                };
            
    uilist = {{'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Plot independent components','fontweight','bold'} ...   
              {'style', 'text', 'string', 'Load IC classification'} ...
              {'style', 'popupmenu', 'string', inputName} ...
              {'style', 'text', 'string', 'To generate IC classification, run tesa_compselect','fontAngle','italic'} ...
              {'style', 'text', 'string', ''} ...
              {'style', 'text', 'string', 'Figure settings','fontweight','bold'} ...
              {'style', 'text', 'string', 'Figure size'} ...
              {'style', 'popupmenu', 'string', figInput, 'tag', 'interp','Value',2} ...
              {'style', 'text', 'string', 'Time course limits (ms): start, end'} ...
              {'style', 'edit', 'string', '-200, 500'} ...
              {'style', 'text', 'string', 'Frequency limits (Hz): start, end'} ...
              {'style', 'edit', 'string', '1, 100'}...
              {'style', 'text', 'string', 'Frequency y-axis scale'} ...
              {'style', 'popupmenu', 'string', freqScaleNames, 'tag', 'freqNames','Value',2} ...
              {'style', 'text', 'string', 'Save IC weights?'} ...
              {'style', 'popupmenu', 'string', 'on|off', 'tag', 'weights', 'value', 2 } ...
              };
          
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find artifact components -- pop_tesa_compplot()', 'helpcom', 'pophelp(''pop_tesa_compplot'')');
    if isempty(result), return; end;
    
    %Extract data
    compClassInput = inputName{result{1}};
    figSize = figInput{result{2}};
    plotTimeX = str2num(result{3});
    plotFreqX = str2num(result{4});
    freqScale = freqScaleNames{result{5}};
    if result{6} == 1
        saveWeights = 'on';
    else
        saveWeights = 'off';
    end
    
end

%Run script from input
if nargin <2
    if strcmp(compClassInput,'')
        EEG = tesa_compplot(EEG,'figSize',figSize,'plotTimeX',plotTimeX,'plotFreqX',plotFreqX,'freqScale',freqScale, 'saveWeights',saveWeights);
        com = sprintf('%s = pop_tesa_compplot( %s,''figSize'',''%s'',''plotTimeX'',%s,''plotFreqX'',%s, ''freqScale'',''%s'', ''saveWeights'',''%s'' );', inputname(1), inputname(1),figSize,mat2str(plotTimeX),mat2str(plotFreqX),freqScale, saveWeights);
    else
        EEG = tesa_compplot(EEG,'compClassInput',EEG.icaCompClass.(compClassInput),'figSize',figSize,'plotTimeX',plotTimeX,'plotFreqX',plotFreqX,'freqScale',freqScale, 'saveWeights',saveWeights);
        com = sprintf('%s = pop_tesa_compplot( %s,''compClassInput'',EEG.icaCompClass.%s,''figSize'',''%s'',''plotTimeX'',%s,''plotFreqX'',%s, ''freqScale'',''%s'', ''saveWeights'',''%s'' );', inputname(1), inputname(1),compClassInput,figSize,mat2str(plotTimeX),mat2str(plotFreqX),freqScale,saveWeights);
    end
elseif nargin > 1
    EEG = tesa_compplot(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_compplot( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

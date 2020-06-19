% pop_tesa_sspsir()-  Uses SSP_SIR method to suppress TMS-evoked muscle artifacts [1] OR
%                     control data [2] from the provided EEGLAB dataset
%
%                     [1] Mutanen, T. P., Kukkonen, M., Nieminen, J. O., Stenroos, M., Sarvas, J.,
%                     & Ilmoniemi, R. J. (2016). Recovering TMS-evoked EEG responses masked by
%                     muscle artifacts. Neuroimage, 139, 157-166.
%
%                     [2] Biabani, M, Fornito, A, Mutanen, T. P., Morrow, J, & Rogasch, N. C.(2019).
%                     Characterizing and minimizing the contribution of sensory inputs to TMS-evoked
%                     potentials.Brain Stimulation, 12(6):1537-1552.
%
%
%
% Usage:
%  >>  [EEG] = pop_tesa_sspsir(EEG); % pop up window; run SSP-SIR using default values
%  >>  [EEG] = pop_tesa_sspsir(EEG, 'key1',value1... );% run SSP-SIR using customised
%                              inputs
%
% Optional input pairs (varargin)
%
%
% 'leadfieldIn',      - Lead field matrix : channel x vertices - If not provided, a simple spherical                     
%                      lead field for the EEGLAB structure will be constructed.
%
% 'leadfieldChansFile'- Leadfield channel matrix (channelName x 1). If provided, sorts the leadfield
%                     channel order to match EEGLAB data (if required).
%
%
% 'artScale',         - Choose the method to estimate the artefacts in the PCA :
%                     'automatic' (default): Uses a sliding window to scale the signal relative to 
%                     the amplitude.
%                     'manual': Only uses the data within a window provided by the user in timeRange
%                     'manualConstant': Works like 'manual' but projects out the artifact dimensions 
%                     estimated from timeRange across the whole time domain. Useful for  the validation 
%                     of SSPSIR.
%                     'control': The only option to remove control data. Uses the  data from the 
%                     control data as provided in EEG_control.
%                     NOTE: 
%                     Setting artScale to 'automatic' only cleans the detected time window of
%                     muscle artifacts.
%                     Setting the artScale to 'manual' only cleans the specified time window of muscle 
%                     artifacts. Setting the artScale to 'manualConstant' estimates muscle artifacts 
%                     from the given time window (e.g. [0 50]), but cleans the dataset indentically                     
%                     across the whole time domain. Option 'manualConstant' is particularly useful for 
%                     controlling how much the neuronal signals of interest might be attenuated by 
%                     SSPSIR.
%
% 'timeRange',        - Required for artScale :'manual' and 'manualConstant'. Vector with start and 
%                     end times in ms of window containing TMS-evoked muscle response.
%                     - Optional for artScale :'control'. Vector with start and end times in ms of 
%                     window containing control responses.
%                     Note that multiple time windows can be defined.
%                     Example: [5,50] - one window
%                     Example: [5,50;100,120] - multiple windows
%                     Multiple windows are useful e.g. when comparing pre- and post datasets and 
%                     attempting to clean the datasets identically.
%                     In such a case the datasets should be concatenated prior to TESA_SSPSIR and 
%                     the time windows should correspond the artifactual window of both of the 
%                     compared datasets.
%
%
% 'PC',               - The number of artifact PCs to be removed. If not provided plots the PC 
%                     dimentions to choose from (default).
%                     [N]: N = The number of principal components to remove
%                     {'data', [N]} : Uses the provided EEG data to find the number of components
%                     which explain N% of the variations.
%                     Example:  {'data', [90]} - removes the components which explain more than 90%
%                     of the variance in the EEG.data
%
%
% 'M',                - The truncation dimension used in the SIR step. If not provided estimates it 
%                     from the EEG data as> rank(data) - PC
%
%
% 'EEG_control'       - EEG data (*.set) from control condition. Required for artScale = 'control'.
%                     [fileName,filePath] Example: ['/Users/myPC/Desktop/controlResponse.set/']
%                     NOTE: EEG_control is expected to be provided in an evoked form.
%
% 
% Outputs:
% EEG                 - EEGLAB EEG structure ( Output of SSPSIR applied on single trials )
%                     NOTE: 
%                     EEG.meanTrials  is the output of SSPSIR applied to the average of all trials
%                     saved on a new field  
%
% Examples:
%  >> [EEG] = pop_tesa_sspsir( EEG ); % default use
%  >> [EEG] = pop_tesa_sspsir( EEG, 'artScale', 'manual','timeRange',[5,50], 'PC',...  
%     {'data', [90]} ); Suppresses muscle artefacts by removing the data components explaining more 
%     than 90% of variance in the time winodw of 5-50ms 
%  >> [EEG] = pop_tesa_sspsir( EEG, 'artScale', 'control','PC',  [5], 'EEG_control',...
%     ['/Users/myPC/Desktop/controlResponse.set/']); % Suppresses control data by removing the first
%     5 principal components of controlResponse.Data from EEG.data 
%
%
% Copyright (C) 2019
% Mana Biabani Mana.Biabanimoghadam@monash.edu, Monash University,
% Tuomas Mutanen tuomas.mutanen@aalto.fi, Aalto University
% Nigel Rogasch,nigel.rogasch@adelaide.edu.au University of Adelaide,
% Jukka Sarvas Jukka.Sarvas@aalto.fi, Aalto University
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
% 30/09/2019

function [EEG, com] = pop_tesa_sspsir( EEG, varargin )

% Check that data is present
if isempty(EEG.data)
    error('Data is empty');
end


PC_select = [ 'if get(gcbo, ''value'') == 1,' ...
    '   set(findobj(gcbf, ''tag'', ''PC''), ''enable'', ''off'');' ...
    'elseif get(gcbo, ''value'') == 2,' ...
    '   set(findobj(gcbf, ''tag'', ''PC''), ''enable'', ''on'', ''string'', ''number of PCs'');'...
    'elseif get(gcbo, ''value'') == 3,' ...
    '   set(findobj(gcbf, ''tag'', ''PC''), ''enable'', ''on'', ''string'', ''% of variance'');'...
    'end;' ];


art_select = [ 'if get(gcbo, ''value'') == 1,' ...
    ' set(findobj(gcbf, ''tag'', ''ctrl''), ''enable'', ''off'');' ...
    ' set(findobj(gcbf, ''tag'', ''txttr''), ''enable'', ''off'');' ...
    ' set(findobj(gcbf, ''tag'', ''tr''),''string'', '''');' ...
    ' set(findobj(gcbf, ''tag'',''pushCtrl''),''enable'', ''off'');' ...
    'elseif get(gcbo, ''value'') == 2,' ...
    ' set(findobj(gcbf, ''tag'', ''ctrl''), ''enable'', ''off'');'...
    ' set(findobj(gcbf, ''tag'', ''txttr''), ''enable'', ''on'');' ...
    ' set(findobj(gcbf, ''tag'', ''tr''),''string'', ''5,50'');' ...
    ' set(findobj(gcbf, ''tag'',''pushCtrl''),''enable'', ''off'');' ...
    'elseif get(gcbo, ''value'') == 3,' ...
    ' set(findobj(gcbf, ''tag'', ''ctrl''), ''enable'', ''off'');'...
    ' set(findobj(gcbf, ''tag'', ''txttr''), ''enable'', ''on'');' ...
    ' set(findobj(gcbf, ''tag'', ''tr''),''string'', ''5,50'');' ...
    ' set(findobj(gcbf, ''tag'',''pushCtrl''),''enable'', ''off'');' ...
    'elseif get(gcbo, ''value'') == 4,' ...
    ' set(findobj(gcbf, ''tag'', ''ctrl''), ''enable'', ''on'');'...
    ' set(findobj(gcbf, ''tag'', ''txttr''), ''enable'', ''on'');' ...
    ' set(findobj(gcbf, ''tag'', ''tr''),''string'', '''');' ...
    ' set(findobj(gcbf, ''tag'',''pushCtrl''),''enable'', ''on'');' ...
    'end;' ];

%fileIn = [ '[filename, filepath] = uigetfile(''*.mat'');' ...
 %   'if filename(1) ~=0,' ...
 %   'set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
 %   'end;' ...
 %   'clear filename filepath tagtest;' ];

fileIn = [ '[filename, filepath] = uigetfile(''*'');' ...
    'if filename(1) ~=0,' ...
    'set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
    'end;' ...
    'clear filename filepath tagtest;' ];

% Pop up window
if nargin < 2
    geometry = {1 [1 0.75 0.25] [1 0.75 0.25] [1 1] [1 1] [1 1] [1 1] [1 1] [1 0.75 0.25] };
    
    uilist = {{'style', 'text', 'string', 'Run SSPSIR','fontweight','bold'} ...
        {'style', 'text', 'string', 'leadfieldIn (optional)'} ...
        {'style', 'edit', 'string', '', 'tag', 'lfFile' } ...
        {'style', 'pushbutton', 'string', '...', 'callback', [ 'tagtest = ''lfFile'';' fileIn ] } ...
        {'style', 'text', 'string', 'leadfieldChansFile (optional)'} ...
        {'style', 'edit', 'string', '', 'tag', 'lfChanFile' } ...
        {'style', 'pushbutton', 'string', '...', 'callback',  [ 'tagtest = ''lfChanFile'';' fileIn ]  }...
        {'style', 'text', 'string', 'artScale'} ...
        {'style', 'popupmenu', 'string', 'automatic|manual|manualConstant|control','value' 1 'callback' art_select }...
        {'style', 'text', 'string', 'timeRange (ms): start,end','enable', 'off', 'tag', 'txttr'} ...
        {'style', 'edit', 'string', '','tag', 'tr'} ...
        {'style', 'text', 'string', 'PC'} ...
        {'style', 'popupmenu', 'string', 'plot|N|data','value' 1 'callback' PC_select}...
        {'style' 'text' 'string' 'number of PCs', 'enable', 'off','tag' 'PC'}...
        {'style' 'edit' 'string' ''} ...
        {'style', 'text', 'string', 'M (optional)'} ...
        {'style', 'edit', 'string', ''} ...
        {'style', 'text', 'string', 'EEG_control (*.set)', 'enable', 'off', 'tag', 'ctrl'} ...
        {'style', 'edit', 'string', '', 'tag', 'ctrlFile'}...
        {'style', 'pushbutton', 'string', '...', 'callback',[ 'tagtest = ''ctrlFile'';' fileIn ],'tag','pushCtrl', 'enable', 'off' }};
    
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Run SSPSIR -- pop_tesa_sspsir()', 'helpcom', 'pophelp(''pop_tesa_sspsir'')');
    
    options = struct('leadfieldIn',[],'leadfieldChansFile',[],'artScale',[],'timeRange',[],'PC',[],'M',[],'EEG_control',[]);
    
    if isempty(result), return; end;
    
    args = {};
    
    % Extract data
    if ~isempty(result{1,1}), args = {'leadfieldIn', result{1,1}}; end %struct2array(load(result{1,1}))}; end;
    if ~isempty(result{1,2}), args = {args{:}, 'leadfieldChansFile', result{1,2}}; end % struct2array(load(result{1,2}))}; end;
    switch  result{1,3}
        case 1,  args = {args{:}, 'artScale','automatic'};
        case 2,  args = {args{:}, 'artScale','manual'};
        case 3,  args = {args{:}, 'artScale','manualConstant'};
        case 4,  args = {args{:}, 'artScale','control'};
    end;
    if ~isempty(result{1,4}), args = {args{:}, 'timeRange', str2num(result{1,4})}; end;
    switch  result{1,5}
        case 1,  args = {args{:},'PC', []};
        case 2,  args = {args{:},'PC',  str2num(result{1,6})};
        case 3,  args = {args{:},'PC',{'data', str2num(result{1,6})}};
    end;
    if ~isempty(result{1,7}), args = {args{:}, 'M',  str2num(result{1,7})}; end;
    if ~isempty(result{1,8}), args = {args{:}, 'EEG_control', result{1,8}}; end;
    
    % Run tesa_sspsir 
    [EEG] = tesa_sspsir(EEG,args{:});   
    com = sprintf('%s = pop_tesa_sspsir( %s, %s);', inputname(1), inputname(1), vararg2str(args) );

elseif nargin > 2
    
    [EEG] = tesa_sspsir(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_sspsir( %s);', inputname(1), inputname(1), vararg2str(varargin)  );
   
end
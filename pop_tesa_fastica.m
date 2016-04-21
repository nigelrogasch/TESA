% pop_tesa_fastica()    - runs FastICA on data using some common settings applied
%                   for TMS-EEG data analysis. tesa_fastica uses EEGLAB pop_runica function. 
%                   See the publications listed below for further details. 
%                   A stabilization option is also included which can help 
%                   if data are not converging. 
% 
%                   Korhonen, Hernandez-Pavon et al (2011) Removal of 
%                   large muscle artifacts from transcranial magnetic
%                   stimulation-evoked EEG by independent component
%                   analysis. Med Biol Eng Compt, 49:397-407.
% 
%                   Rogasch et al (2014) Removing artefacts from TMS-EEG 
%                   recordings using independent component analysis: 
%                   Importance for assessing prefrontal and motor cortex 
%                   network properties. NeuroImage, 101: 429-435.
% 
%                   Note that this script requires the FastICA algorithm to
%                   be included in the matlab path. The package can be
%                   downloaded from:
%                   http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml                   
%                   
% Usage:
%   >>  EEG = pop_tesa_fastica( EEG ); % pop up window
%   >>  EEG = pop_tesa_fastica( EEG, 'key1',value1... ); 
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs
%   'approach','str' - 'symm' | 'defl'. Symmetric or deflation approach for
%                   FastICA. The symmetric approach is more reliable and therefore 
%                   highly recommended. See Korhonen et al for details.
%                   Default: 'symm'
%   'g','str'       - 'tanh' | 'gauss' | 'pow3' | 'skew'. Contrast function
%                   for FastICA. Either 'tanh' or 'gauss' perform equally
%                   well for TMS-EEG analysis. See Korhonen et al for
%                   details.
%                   Default: 'tanh'
%   'stabilization','str' - 'on' | 'off'. Controls whether FastICA uses
%                   stabilized version which detects 'strokes' (i.e. when
%                   algorithm gets stuck between 2 points and won't
%                   converge) and halves the learning rate.
%                   Default: 'off'
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = pop_tesa_fastica( EEG ); %default use
%   EEG = pop_tesa_fastica( EEG, 'g','gauss','stabilization','on' ); % Uses the gauss contrast function and turns on the stabilized FastICA version to aid with convergence.
% 
% See also:
%   eegfiltnew

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

 function [EEG com] = pop_tesa_fastica( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
       
    geometry = {1 [1 1] [1 1] [1 1]};

    uilist = {{'style', 'text', 'string', 'Run FastICA','fontweight','bold'} ...
              {'style', 'text', 'string', 'Approach'} ...
              {'style', 'popupmenu', 'string', 'symmetric|deflation', 'tag', 'interp' } ...
              {'style', 'text', 'string', 'Contrast function (g)'} ...
              {'style', 'popupmenu', 'string', 'tanh|gauss|pow3|skew', 'tag', 'interp' } ...
              {'style', 'text', 'string', 'stabilization'} ...
              {'style', 'popupmenu', 'string', 'off|on', 'tag', 'interp' }};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Run FastICA -- pop_tesa_fastica()', 'helpcom', 'pophelp(''pop_tesa_fastica'')');
    if isempty(result), return; end;
    
    %Extract data 
    if result{1,1} == 1
        approach = 'symm';
    elseif result{1,1} == 2
       	approach = 'defl';
    end
    if result{1,2} == 1
        g = 'tanh';
    elseif result{1,2} == 2
       	g = 'gauss';
    elseif result{1,2} == 3
       	g = 'pow3';
    elseif result{1,2} == 4
       	g = 'skew';
    end
    if result{1,3} == 1
        stabilization = 'off';
    elseif result{1,3} == 2
       	stabilization = 'on';
    end
    
end

%Run script from input
if nargin <2
    EEG = tesa_fastica(EEG,'approach',approach,'g',g,'stabilization',stabilization);
    com = sprintf('%s = pop_tesa_fastica( %s, ''approach'', ''%s'', ''g'', ''%s'', ''stabilization'', ''%s'' );', inputname(1), inputname(1), approach, g, stabilization );
elseif nargin > 2
    EEG = tesa_fastica(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_fastica( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

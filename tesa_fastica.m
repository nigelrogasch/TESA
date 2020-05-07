% tesa_fastica()    - runs FastICA on data using some common settings applied
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
%   >>  EEG = tesa_fastica( EEG ); %default use
%   >>  EEG = tesa_fastica( EEG, 'key1',value1... ); 
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
%   EEG = tesa_fastica( EEG ); %default use
%   EEG = tesa_fastica( EEG, 'g','gauss','stabilization','on' ); % Uses the gauss contrast function and turns on the stabilized FastICA version to aid with convergence.
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

% Change log
% 13.6.2018 - Included rank and sort of components with tesa_sortcomps

function EEG = tesa_fastica( EEG, varargin )

if nargin < 1
	error('Not enough input arguments.');
end

%Check inputs
%define defaults
options = struct('approach','symm','g','tanh','stabilization','off');

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

% Checks inputs
if ~(strcmp(options.approach,'symm') || strcmp(options.approach,'defl'))
    error('Input for ''approach'' must be either ''symm'' or ''defl''. Input type ''symm'' is recommended for TMS-EEG. See help for further details.');
elseif ~(strcmp(options.g,'tanh') || strcmp(options.g,'gauss') || strcmp(options.g,'pow3') || strcmp(options.g,'skew'))
    error('Input for ''g'' must be either ''tanh'', ''gauss'',''pow3'' or ''skew''. Either ''tanh'' or ''gauss'' are recommended for TMS-EEG. See help for further details.');
elseif ~(strcmp(options.stabilization,'on') || strcmp(options.stabilization,'off'))
    error('Input for ''stabilization'' must be either ''on'' or ''off''.')
end

% Save random number generator settings and date/time
tmpfield = fieldnames(EEG);
tmplog = strncmp('fastica',tmpfield,7);
tmpfastica = tmpfield(tmplog);
if isempty(tmpfastica)
    EEG.fastica1.rng = rng;
    EEG.fastica1.date = datetime('now');
else
    tmpname = ['fastica',num2str(length(tmpfastica)+1)];
    EEG.(tmpname).rng = rng;
    EEG.(tmpname).date = datetime('now');
end

% Run FastICA using EEGLAB pop_runica
EEG = pop_runica(EEG,'icatype','fastica', 'approach', options.approach, 'g', options.g,'stabilization',options.stabilization);

% Ranks and sorts components
EEG = tesa_sortcomps(EEG);

%display message
fprintf('FastICA performed on data using ''%s'' approach and ''%s'' contrast function.\n',options.approach,options.g);
    

end

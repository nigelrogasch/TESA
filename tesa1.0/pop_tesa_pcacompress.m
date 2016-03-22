% pop_tesa_pcacompress()    - Compresses data to n-dimensions as advocated in
%                           the following papers:
%
%                           Korhonen, hernandez-pavon et al (2011) Removal of large muscle artifacts 
%                           from transcranial magnetic stimulation-evoked EEG 
%                           by independent component analysis. Med Biol Eng
%                           Compt, 49:397-407.
%
%                           Hernandez-Pavon et al (2012) Uncovering neural  
%                           independent components from highly artifactual 
%                           TMS-evoked EEG data. J Neurosci Meth,
%                           209:144-57
% 
% Usage:
%   >>  EEG = pop_tesa_pcacompress( EEG ); % pop up window
%   >>  EEG = pop_tesa_pcacompress( EEG , 'key1', value1...);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
% 
% Optional input pairs:
%   'compVal',int       - int = number of dimensions to compress data
%                       Default = 30.
%   'plot','str'        - 'on' | 'off'. Turns on/off plot summarising the
%                       variance explained by principal components.
%                       Default = 'on'
% 
% Outputs:
%   EEG                 - EEGLAB EEG structure
%
% Examples:
%   EEG = pop_tesa_pcacompress( EEG );
%   EEG = pop_tesa_pcacompress( EEG, 'compVal', 30 ); %compress to top 30 dimensions
%   EEG = pop_tesa_pcacompress( EEG, 'plot','off' ); %turns off summary plot
% 
% See also:
%   pop_tesa_pca_supress, pop_tesa_edm 

% Copyright (C) 2016  Nigel Rogasch & Julio C. Hernandez-Pavon
% Monash University and Aalto University
% nigel.rogasch@monash.edu; julio.hpavon@gmail.com
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

 function [EEG com] = pop_tesa_pcacompress( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {[1 0.3] [1 0.3]};
    
    uilist = {{'style', 'text', 'string', 'Number of dimesions for compression'} ...
              {'style', 'edit', 'string', '30'}...
              {'style', 'text', 'string', 'Plot dimension ranks?'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              };
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Compress data dimensionality-- pop_tesa_pcacompress()', 'helpcom', 'pophelp(''pop_tesa_pcacompress'')');
    if isempty(result), return; end;
    
    %extract compression value
    compVal = str2num(result{1,1});
    if result{1,2} == 1
        plot = 'on';
    else
        plot ='off';
    end
    
end

%remove data
if nargin < 2
    EEG = tesa_pcacompress( EEG, 'compVal', compVal, 'plot', plot );
    com = sprintf('%s = pop_tesa_pcacompress( %s, ''compVal'', %s, ''plot'',''%s'' );', inputname(1), inputname(1), mat2str(compVal), plot);
else
    EEG = tesa_pcacompress( EEG, varargin{:} );
    com = sprintf('%s = pop_tesa_pcacompress( %s, %s );', inputname(1), inputname(1), vararg2str(varargin));
end

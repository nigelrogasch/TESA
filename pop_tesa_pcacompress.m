% pop_tesa_pcacompress()  - Compresses data to n-dimensions as advocated in
%                           the following papers:
%
%                           Korhonen et al (2011) Removal of large muscle artifacts 
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
%   >>  EEG = tesa_sortcomps( EEG ); % pop-window
%   >>  EEG = tesa_sortcomps( EEG , compVal);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   compVal             - [Integer] Number of dimensions to compress data
%                           Default = 25.
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% See also:
%   SAMPLE, EEGLAB 
%
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

 function [EEG com] = pop_tesa_pcacompress( EEG, compVal )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {[1 0.5]};
    
    uilist = {{'style', 'text', 'string', 'Number of dimesions for compression'} ...
              {'style', 'edit', 'string', '25'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Compress data dimensionality-- pop_tesa_pcacompress()', 'helpcom', 'pophelp(''tesa_pcacompress'')');
      
    %extract compression value
    compVal = str2num(result{1,1});
    
end

%remove data
EEG = tesa_pcacompress( EEG, compVal );

% return the string command

com = sprintf('%s = pop_tesa_pcacompress( %s, %s );', inputname(1), inputname(1), compVal);

end

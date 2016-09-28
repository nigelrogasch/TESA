% tesa_pcacompress()    - Compresses data to n-dimensions as advocated in
%                           the following papers:
%
%                           Korhonen, Hernandez-pavon et al (2011) Removal of large muscle artifacts 
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
%   >>  EEG = tesa_pcacompress( EEG );
%   >>  EEG = tesa_pcacompress( EEG , 'key1', value1...);
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
%   EEG = tesa_pcacompress( EEG );
%   EEG = tesa_pcacompress( EEG, 'compVal', 30 ); %compress to top 30 dimensions
%   EEG = tesa_pcacompress( EEG, 'plot','off' ); %turns off summary plot
% 
% See also:
%   tesa_pca_supress, tesa_edm 

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

function [EEG,inMat,outMat] = tesa_pcacompress( EEG, varargin )

if nargin < 1
	error('Not enough input arguments.');
end

%define defaults
options = struct('compVal',30,'plot','on');

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

%Check compVal input is correct
if ischar(options.compVal)
    error('''CompVal'' must be an integer, not a string. e.g. ''compVal'', 30.')
end
    

%Check plot input is correct
if ~(strcmp(options.plot,'on') || strcmp(options.plot,'off'))
    error('''Plot must be either ''on'' or ''off'', e.g. ''plot'',''on''');
end

%Reshapes 3D matrix to 2D
inMat=reshape(EEG.data,size(EEG.data,1),[],1);

%Checks that number of dimensions is larger than compression value
covarianceMatrix = cov(inMat', 1);
[E, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
rankMat = sum (diag (D) > rankTolerance);

if rankMat <= options.compVal
    error('Dimension of data (%d) is lower than the compression value (%d). Function terminated.', rankMat, options.compVal);
end


%Runs singular value decomposition
[U,S,V]=svd(inMat*inMat');
d=diag(S);

%Compresses data by truncating to 'compVal' dimensions 
C=U(:,1:options.compVal)*U(:,1:options.compVal)';
outMat=C*inMat;

%Reshapes 2D matrix to 3D matrix
EEG.data = reshape(outMat,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));

%Figures before and after compression
if strcmp(options.plot,'on')
    h1=figure; subplot(1,2,1), bar(d);grid;title('Data before compression'); 
    xlabel ('Dimensions'); ylabel ('Amplitude (A.U)');
    subplot(1,2,2), bar(d(1:options.compVal));grid;title('Data after compression'); 
    xlabel ('Dimensions'); ylabel ('Amplitude (A.U)');
end

fprintf('Data compressed to %d dimensions\n', options.compVal);
        
end

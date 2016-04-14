% tesa_icaedm( )               - This function computes ICA by using the 
%                           enhanced deflation method (EDM)  
%                           as advocated in the following papers
%
%                           Korhonen, Hernandez-Pavon et al (2011) 
%                           Removal of large muscle artifacts 
%                           from transcranial magnetic stimulation-evoked EEG 
%                           by independent component analysis. Med Biol Eng
%                           Compt, 49:397-407.
%

% Usage:
%   >>  [Wbest]=tesa_icaedm(EEGwhite,Nic); %Function to compute EDM;
%
% Inputs:
%   EEGwhite            - Whitened data
%   Nic                 - [Integer] Number of independent components to look for
% 
% Outputs:
%   Wbest               - Best Weight vectors 
%
% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com

function [W]=tesa_icaedm(Xwhite,N)
% Xwhite = Whitened Data matrix
% N= Number of independent components

itercount=10; % Number of iterations in finding a candidate weight vector w 
candcount=5;  % Number of candidates for choosing the best w 
x=Xwhite;
M=size(x,1);
T=size(x,2);
P=eye(M);
W=[];
Neg=[];
for j=1:N
    nege=0;
    for jj=1:candcount
        w=randn(M,1);
        w=P*w;
        w=w/norm(w);
        count=0;
        while count<itercount
            t=w'*x; % a row vector
                g=tanh(t);
                gp=1-g.^2;
            Eg=1/T*sum(x.*(ones(M,1)*g),2);
            wplus=Eg-1/T*sum(gp,2)*w;
            w=P*wplus;
            w=w/norm(w);
            count=count+1;
        end
        s=w'*x;
          Jneg=abs(1/T*sum(log(cosh(s)),2)-0.374567207491438);%G(t)=log(cosh(t))
        if Jneg>nege
            nege=Jneg;
            wbest=w;
        end
    end
    w=wbest;
    P=P-w*w';
    W=[W,w];
    Neg=[Neg;nege];
end

% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This software is the implementation of the Ariadne algorithm by Nasso et al.
% Copyright (C) 2012 Sara Nasso
%
% Ariadne is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at
% your option) any later version.
%
% Ariadne is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mProphetOutput falseHitsmProphet]=filter_mProphetOutput(filemProphetOutput,FDRthreshold)
% it filters mProphet output results, user definable is the FDR threshold,
% suggested values are 0.01 or 0.05

global settings

mProphetOutput=readTxt(filemProphetOutput);
idxCol=strmatch('decoy',mProphetOutput(1,:),'exact');
idxsDecoys=find(ismember(mProphetOutput(1:end,idxCol),{'TRUE'})==1);
idxCol=strmatch('peak_group_rank',mProphetOutput(1,:),'exact');
idxsNotBestPeakGroup=find(ismember(mProphetOutput(1:end,idxCol),{'1'})==0);
idxCol=strmatch('m_score',mProphetOutput(1,:),'exact');
m_score=(mProphetOutput(2:end,idxCol));
if settings.mProphetFilter
    boolFDR=zeros(size(m_score,1)+1,1);
    for i=1:size(m_score,1)
        tmpScore=str2num(m_score{i});
        if (~isempty(tmpScore)) & ~isequal(tmpScore,'NA') & ~isequal(tmpScore,'NaN') & ~isequal(tmpScore,'infinite') & tmpScore<FDRthreshold % ok is 1
            boolFDR(i+1)=1;
        end
    end
elseif ~settings.mProphetFilter
    boolFDR = ones(size(m_score,1)+1,1);
end
idxsFDR=find(boolFDR==0);
idxsToTrash=setdiff(union([idxsDecoys;idxsNotBestPeakGroup],idxsFDR),1); % row 1 is the header

falseHitsmProphet=mProphetOutput(idxsToTrash,:); % saving false hits according to mProphet in a separate variable
mProphetOutput(idxsToTrash,:)=[]; % getting rid of false hits

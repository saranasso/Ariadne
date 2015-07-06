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
function [idx_Data, Q3s, rts, trans]=getPrecursorData(srm, precursor, Q3s)
% it returns the index position of the current precursor and its data from
% the SRM data structure
global settings 

[min_dist idx_min]=min(abs([srm.data.Q1]-precursor));
if min_dist<settings.ionsResolution
    idx_Data=idx_min;
else
    disp(['Precursor ' num2str(precursor) ' not in ' srm.metadata.Filename])
    return
end
rts=[];
trans=[];

for i=1:numel(Q3s)
    bool=any(abs(Q3s(i)-srm.data(idx_Data).Q3)<settings.ionsResolution);
    idx(i)=find(abs(Q3s(i)-srm.data(idx_Data).Q3)<settings.ionsResolution);
end

if bool
    if settings.filterTrans || ~settings.mProphetFilter
        rts=srm.data(idx_Data).rts;
        trans=srm.data(idx_Data).trans(idx);
    elseif ~settings.filterTrans
        rts=srm.data(idx_Data).rts;
        trans=srm.data(idx_Data).trans;
    end
elseif ~ bool
    disp(['Q3s from trans list do not match to Q3s from data for precursor: ' num2str(precursor)])
    return
end	

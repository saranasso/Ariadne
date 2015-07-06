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
function [data currentRts currentTrans]=getPeakData(rts,trans,PExt,idx_peak,idx)
% it selects data from a certain peak and gives back also the occurrences
% data

if idx_peak~=0
    if numel(idx_peak)==1
        peak_in_rt=PExt(idx_peak,1);
        peak_fi_rt=PExt(idx_peak,2);
    else
        peak_in_rt=PExt(idx_peak(idx),1);
        peak_fi_rt=PExt(idx_peak(idx),2);
        disp('warning: more than one index')
    end
    idxs_data_to_keep=find(rts>peak_in_rt & rts<peak_fi_rt);
    currentTrans=trans{idx}(idxs_data_to_keep);
    currentRts=rts(idxs_data_to_keep);
    data=hist2data(currentTrans,currentRts);
else
    idxs_data_to_keep=1:1:numel(trans{idx});
    currentTrans=trans{idx}(idxs_data_to_keep);
    currentRts=rts(idxs_data_to_keep);
    data=hist2data(currentTrans,currentRts);
end

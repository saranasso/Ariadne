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
function [endo ref rts_endo rts_ref] = getPeakGroupWindow(rti,rtf,rts,trans)
% it cuts out the ROI around the selected peak group

[rti_endo rti_ref rtf_endo rtf_ref]=getRange(rti,rtf,rts);
idxs_endo=find(rts.endo>rti_endo & rts.endo<rtf_endo);
idxs_ref=find(rts.ref>rti_ref & rts.ref<rtf_ref);
N=round(rtf-rti); %+1
idxs_common=intersect(idxs_endo,idxs_ref);
idxs_endo=idxs_common;
idxs_ref=idxs_common;

for i=1:size(trans.ref,2)
    endo{i}=trans.endo{i}(idxs_endo)';
    ref{i}=trans.ref{i}(idxs_ref)';
    rts_endo{i}=rts.endo(idxs_endo)';
    rts_ref{i}=rts.ref(idxs_ref)';
end	

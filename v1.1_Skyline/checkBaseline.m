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
function [endo, ref, rts_endo,rts_ref, rtsFilt, trans] = checkBaseline(endo,ref,rts_endo,rts_ref,trans,rts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added for QTRAP data
baseline=0;
for i=1:size(trans.ref,2)
    idxs_endo=find(endo{i}>baseline);
    idxs_ref=find(ref{i}>baseline);
    endo{i}=endo{i}(idxs_endo);
    ref{i}=ref{i}(idxs_ref);
    rts_endo{i}=rts_endo{i}(idxs_endo);
    rts_ref{i}=rts_ref{i}(idxs_ref);
    idxs_endo=find(trans.endo{i}>baseline);
    idxs_ref=find(trans.ref{i}>baseline);
    rtsFilt.endo{i}=rts.endo(idxs_endo)';
    rtsFilt.ref{i}=rts.ref(idxs_ref)';
    trans.endo{i}=trans.endo{i}(idxs_endo)';
    trans.ref{i}=trans.ref{i}(idxs_ref)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

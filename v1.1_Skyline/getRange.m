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
function [rti_endo rti_ref rtf_endo rtf_ref]=getRange(rti,rtf,rts)
% it defines the range on retention times for data analysis

% obvious hp: sorted vectors for retention times
% conservative version to get identical ret times vector
if rti<max(rts.endo(1),rts.ref(1))
    rti_endo=max(rts.endo(1),rts.ref(1));
    rti_ref=max(rts.endo(1),rts.ref(1));
else
    rti_endo=rti;
    rti_ref=rti;
end

if rtf>min(rts.endo(end),rts.ref(end))
    rtf_endo=min(rts.endo(end),rts.ref(end));
    rtf_ref=min(rts.endo(end),rts.ref(end));
else
    rtf_endo=rtf;
    rtf_ref=rtf;
end

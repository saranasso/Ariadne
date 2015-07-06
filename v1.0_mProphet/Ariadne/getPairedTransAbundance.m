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
function   abundance_estimate=getPairedTransAbundance(rts,trans,PExt,idx_peak,idx,estimatedPeakWidth)
% statistically assesses peak boundaries at 1% and estimates trans
% abundance

[data.endo currentRts.endo currentTrans.endo]=getPeakData(rts.endo,trans.endo,PExt.endo,idx_peak.endo,idx);
[data.ref currentRts.ref currentTrans.ref]=getPeakData(rts.ref,trans.ref,PExt.ref,idx_peak.ref,idx);

[currentTrans,currentRts]=synchTrans(currentTrans,currentRts);

dataMerged = cat(1,data.endo,data.ref);
rts_bounds=quantile(dataMerged,[.01 .99]);
subplot 211
abundance_estimate.endo=getAbundanceEstimate(rts_bounds,currentRts.endo,currentTrans.endo,dataMerged);
subplot 212
abundance_estimate.ref=getAbundanceEstimate(rts_bounds,currentRts.ref,currentTrans.ref,dataMerged);

	

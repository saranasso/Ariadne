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
function [endo, ref, cut] = intCut (endo, ref, method,input)
% reducing the total number of counts for speeding up the GMM estimate. It
% doesn't affect the relative intensities for the GMM (areas ratio invariant to this transformation)

switch method
    
    case 'maxCountsPerTrans'
        
        maxCounts = input * size(endo,2)*2; % both endo and ref
        
        sumCounts = 0;
        for j = 1:size(endo,2);
            sumCounts = sumCounts + sum(endo{j})+ sum(ref{j});
        end
        
        cut = sumCounts/maxCounts;
        
    case 'minCountIntensity'
        
        cut = input;
        
    case 'estimateConservative'
        
        [medianCount.endo] = getMedianCount(endo);
        [medianCount.ref] = getMedianCount(ref);
        cut = min([medianCount.endo, medianCount.ref]);
               
end



for j = 1:size(endo,2);
    ref{j} = ref{j}./cut;
    endo{j} = endo{j}./cut;
end	

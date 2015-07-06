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
function data=hist2data(hist_counts,rts)
% reverts histograms to data

tot_counts=(sum(ceil(hist_counts(hist_counts>1))));
data=ones(tot_counts,1);
counts_done=0;
if nargin==1
    rts=(1:1:length(hist_counts));
end
for i=1:length(hist_counts)
        if hist_counts(i)>1
            counts=ceil(hist_counts(i));
            data(counts_done+1:counts_done+counts)=rts(i)*ones(counts,1);
            counts_done=counts_done+counts;
        end
end
	

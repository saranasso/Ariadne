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
function [retTime idxPG gm  X minDist]=pickPeakGroup(Peaklist,PExt,prior,endo,rts,pepLibraryInfo)
% it fits the GMM model to the data and finds the putative retention time
% that is associated to the best peak group for the monitored target peptide 

NlogLFitted=[];
X=[];
dist=[];
if any(prior.PComponents==0)
   prior.PComponents(find(prior.PComponents==0)) = 0.00001;
end
gm = gmdistribution(prior.mu,prior.Sigma,prior.PComponents);
gmmFitted = cell(1,size(Peaklist,1));
X = cell(1,size(Peaklist,1));
dist = [];
for k=1:size(Peaklist,1)
    rti=PExt(k,1);
    rtf=PExt(k,2);
    [gmmFitted{k}, X{k}, dist(k)] = fitPeakGroupsGMM(rti,rtf,gm,prior,endo,rts,pepLibraryInfo);  
end

[minDist, idxPG]=min(dist);

retTime=Peaklist(idxPG,1);

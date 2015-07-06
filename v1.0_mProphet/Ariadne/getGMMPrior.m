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
function prior=getGMMPrior(pepLibraryInfo,settings,numOfTrans, method)
% it defines the GMM prior using the "a priori" expected values of its parameters

prior.mu=[];
prior.Sigma=[];
prior.PComponents=[];
switch method
    
    case 'endo'
        for i=1:numOfTrans
            prior.mu = cat(1,[pepLibraryInfo.RT_time pepLibraryInfo.Q3s(i)],prior.mu);
            prior.Sigma = cat(3,[settings.expPeakWidth/6 0; 0 1],prior.Sigma); % divided by 6 because 99,7% lay in 3*sigma per each tail=6*sigma
            prior.PComponents = cat(1,pepLibraryInfo.relInt(i),prior.PComponents); % no need for normalizing the rel int as real mix props
        end
        
    case 'ref'
        for i=1:numOfTrans
            prior.mu = cat(1,[pepLibraryInfo.RT_time pepLibraryInfo.Q3s(i)],prior.mu);
            prior.Sigma = cat(3,[settings.expPeakWidth/6 0; 0 1],prior.Sigma); % divided by 6 because 99,7% lay in 3*sigma per each tail=6*sigma
            prior.PComponents = cat(1,pepLibraryInfo.relIntRef(i),prior.PComponents); % no need for normalizing the rel int as real mix props
        end
        
    case 'merge' 
        for i=1:numOfTrans
            prior.mu = cat(1,[pepLibraryInfo.RT_time pepLibraryInfo.Q3s(i)],prior.mu);
            prior.Sigma = cat(3,[settings.expPeakWidth/6 0; 0 1],prior.Sigma); % divided by 6 because 99,7% lay in 3*sigma per each tail=6*sigma
            prior.PComponents = cat(1,pepLibraryInfo.relIntRef(i)*pepLibraryInfo.relInt(i),prior.PComponents); % no need for normalizing the rel int as real mix props
        end
        
end	

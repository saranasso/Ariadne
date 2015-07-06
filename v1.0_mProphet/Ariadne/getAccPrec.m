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
function  [RE sd qcd]=getAccPrec(estimatedValues,trueValues)

estimatedValues(isinf(estimatedValues)) = nan;

for i=1:size(trueValues,1)
    RE(i,1)=nanmedian(abs(estimatedValues(:,i)-trueValues(i)))/trueValues(i)*100;
    sd(i,1)=mad(estimatedValues(:,i),1);
%     qcd_mod(i,1)=((iqr(estimatedValues(:,i))/nanmedian(estimatedValues(:,i))))*100;
%     cv(i,1)=(sd(i,1)/nanmedian(estimatedValues(:,i)))*100; 
    qcd(i,1)=((iqr(estimatedValues(:,i)))/(prctile(estimatedValues(:,i),75)+prctile(estimatedValues(:,i),25)))*100; %QCD

end



	

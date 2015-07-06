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

%% Update pepLibrary if you want to quantify additional peps/transitions you may have added to your transition list. 
%  Do not forget relative intensities! After running this, run again analyze, it
%  will attempt to quantify those without a quantification value

[pepLibrarySRM]=updatePepLibFromTransList(pepLibrarySRM,transList)
for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).dataFile=association_file(:,2);
end

if settings.mProphet && ~settings.mProphetFilter
    [ pepLibrarySRM pepInListNotIn_mProphet]=getmProphetResults(pepLibrarySRM,transList,mProphetOutput,falseHitsmProphet)
end

if settings.mProphet && settings.mProphetFilter
    [pepLibrarySRM pepInListNotIn_mProphet]=matchTransListTomProphetOut(transList,mProphetOutput,falseHitsmProphet)
end
pepLibrarySRM(find([pepLibrarySRM(:).hasRef]==0))=[]
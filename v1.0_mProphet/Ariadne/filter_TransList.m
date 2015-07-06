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
function [transList decoysTransList]=filter_TransList(fileTransList)
% it combines all trans lists and filters out decoys from input transition list

transList=[];
for i=1:numel(fileTransList)
    eval(['[transList_tmp header]=readTxt(fileTransList{' num2str(i) '});'])
    transList=cat(1,transList,transList_tmp(2:end,:));
end
transList=cat(1,transList_tmp(1,:),transList);

idxColDecoy=strmatch('decoy',transList(1,:),'exact');
idxsDecoysTransList=find(ismember(transList(1:end,idxColDecoy),{'1'})==1);
decoysTransList=transList(idxsDecoysTransList,:); % saving decoys in another variable
transList(idxsDecoysTransList,:)=[]; % getting rid of decoys
	

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

global settings

try
    association_dil_file = association_file; % necessary to cope with different versions
end
molarities=sort(unique([association_dil_file{:,1}]),'descend');
molarities(molarities>1)=round(molarities(molarities>1));

quants= {'Protein','Peptide','Charge','Retention T','FileName','Ratio','Absolute Quant'};
for i=1:numel(pepLibrarySRM)
    if pepLibrarySRM(i).hasRef
        for j = 1:numel(pepLibrarySRM(i).ratio)
            if ~isnan(pepLibrarySRM(i).ratio(j))
                
                if settings.absoluteQuantification
                    grpID = association_file{strmatch(pepLibrarySRM(i).dataFile{j}, association_file(:,2),'exact'),1};
                    grpIdx = find(ismember(molarities,grpID)==1);
                    currentAbsQuant = pepLibrarySRM(i).estimatedAbundance(grpIdx);
                    if isnan(currentAbsQuant) && pepLibrarySRM(i).ratio(j)==0
                        continue
                    end
                else
                    currentAbsQuant = 'NA';
                end
                quants=cat(1,quants,{pepLibrarySRM(i).protein,pepLibrarySRM(i).sequence,pepLibrarySRM(i).charge,pepLibrarySRM(i).estRetTime(j),pepLibrarySRM(i).dataFile{j},pepLibrarySRM(i).ratio(j),currentAbsQuant});
            else
                continue
            end
        end
    end
end

export(cell2dataset(quants),'File',[pwd filesep 'AriadneQuants.txt'],'Delimiter','\t')

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
function [pepLibrarySRM]=getCalibration(pepLibrarySRM,association_file,idx_library)
% it estimates the calibration curve for each peptide

global settings
i=idx_library;
[idx_bool idx_ht]=ismember([association_file(:,2)],[pepLibrarySRM(i).dataFile]);
molarities_files=[cell2mat(association_file(idx_bool,1)) idx_ht(idx_bool)];
molarities=sort(unique([association_file{:,1}]),'descend');
standard_ratios=[];
k=1;
for j=1:numel(molarities)
    idx_mol=molarities_files(:,1)==molarities(j);
    if any(idx_mol)
        idxs_cc_dataFiles = intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio],2)); % if some data files were not quantified then the variable is truncated at a shorter position
        if strcmp(settings.constantTarget, 'light')
            standard_ratios(k,1)=nanmedian([pepLibrarySRM(i).ratio(molarities_files(idxs_cc_dataFiles,2))].^(-1)); % heavy to light needed here cause light is constant and unknown while heavy is variable and known
        else strcmp(settings.constantTarget, 'heavy')
            standard_ratios(k,1)=nanmedian([pepLibrarySRM(i).ratio(molarities_files(idxs_cc_dataFiles,2))]); % heavy to light needed here cause light is constant and unknown while heavy is variable and known
        end
        conc_tmp(k,1)=molarities(j);
        k=k+1;
    else
        standard_ratios(k,1)=NaN;
    end
    clear idx_mol
end
conc_tmp(isnan(standard_ratios))=[];
standard_ratios(isnan(standard_ratios))=[];
conc_tmp(isinf(standard_ratios))=[];
standard_ratios(isinf(standard_ratios))=[];
conc_tmp(standard_ratios==0)=[];
standard_ratios(standard_ratios==0)=[];
if numel(standard_ratios)>1
    [pepLibrarySRM(i).regression]=getCalibrationCurve(standard_ratios,conc_tmp);
else
    return
end
[r pepLibrarySRM(i).R2_pValue.Ariadne] = corr(standard_ratios,conc_tmp,'tail','gt','rows','pairwise');
pepLibrarySRM(i).R2.Ariadne=r^2;



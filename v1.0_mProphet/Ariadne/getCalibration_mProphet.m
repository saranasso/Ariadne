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
function [pepLibrarySRM]=getCalibration_mProphet(pepLibrarySRM,association_file,idx_library)
% it estimates the calibration curve for each peptide

global settings
i=idx_library;
[idx_bool idx_ht]=ismember([association_file(:,2)],[pepLibrarySRM(i).dataFile]);
molarities_files=[cell2mat(association_file(idx_bool,1)) idx_ht(idx_bool)];
molarities=sort(unique([association_file{:,1}]),'descend');

mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};
standard_ratios=[];
k=1;
for j=1:numel(molarities)
    idx_mol=molarities_files(:,1)==molarities(j);
    if any(idx_mol)
        for numMethod=1:numel(mProphet_methods)
            if isfield(pepLibrarySRM(i).ratio_mProphet, mProphet_methods{numMethod})
                idxs_cc_dataFiles=intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio_mProphet],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                idxs_cc_dataFiles = intersect(idxs_cc_dataFiles,1:1:size(pepLibrarySRM(i).mScore,2));
                idxs_cc_dataFiles = idxs_cc_dataFiles(find(pepLibrarySRM(i).mScore(idxs_cc_dataFiles)<settings.FDRthreshold));
                if strcmp(settings.constantTarget, 'light')
                    eval(['standard_ratios.' mProphet_methods{numMethod}  '(k,1)=nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod} ' ].^(-1));']); % heavy to light needed here cause light is constant and unknown while heavy is variable and known
                else strcmp(settings.constantTarget, 'heavy')
                    eval(['standard_ratios.' mProphet_methods{numMethod}  '(k,1)=nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod} ' ]);']); % light to heavy  here cause heavy is constant and unknown while light is variable and known
                end
            end
        end
        conc_tmp(k,1)=molarities(j);
        k=k+1;
    end
    clear idx_mol
end

options={'tail','gt','rows','pairwise'};
for numMethod=1:numel(mProphet_methods)
    if isfield(standard_ratios, mProphet_methods{numMethod})
        eval(['standard_ratios_tmp=standard_ratios.' mProphet_methods{numMethod} ';'])
    else
        continue
    end
    conc_2xtmp=conc_tmp;
    subplot(2,2,numMethod)
    conc_2xtmp(isnan(standard_ratios_tmp))=[];
    standard_ratios_tmp(isnan(standard_ratios_tmp))=[];
    conc_2xtmp(isinf(standard_ratios_tmp))=[];
    standard_ratios_tmp(isinf(standard_ratios_tmp))=[];
    conc_2xtmp(standard_ratios_tmp==1)=[];
    standard_ratios_tmp(standard_ratios_tmp==1)=[];
    conc_2xtmp(standard_ratios_tmp==0)=[];
    standard_ratios_tmp(standard_ratios_tmp==0)=[];
    if numel(standard_ratios_tmp)>1
        eval(['[pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod} ']=getCalibrationCurve(standard_ratios_tmp,conc_2xtmp);'])
    else
        return
    end
    title(['mProphet method ' mProphet_methods{numMethod}])
    eval(['[r pepLibrarySRM(i).R2_pValue.' mProphet_methods{numMethod} '] = corr(standard_ratios_tmp,conc_2xtmp,options{1},options{2},options{3},options{4});'])
    eval(['pepLibrarySRM(i).R2.' mProphet_methods{numMethod} '=r^2;'])
    clear conc_2xtmp
end

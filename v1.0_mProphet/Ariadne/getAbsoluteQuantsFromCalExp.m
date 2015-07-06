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
function [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExp(pepLibrarySRM,association_dil_file,pepLibrarySRMCal)

global settings
mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};

for i=1: numel(pepLibrarySRM)
    %     if pepLibrarySRM(i).hasRef
    [idx_bool idx_ht]=ismember([association_dil_file(:,2)],[pepLibrarySRM(i).dataFile]);
    molarities_files=[cell2mat(association_dil_file(idx_bool,1)) idx_ht(idx_bool)];
    molarities=sort(unique([association_dil_file{:,1}]),'descend');
    for j=1:numel(molarities)
        idx_mol=molarities_files(:,1)==molarities(j);
        idx_lib_cal =  intersect(find(ismember({pepLibrarySRMCal(:).sequence},pepLibrarySRM(i).sequence)==1),find([pepLibrarySRMCal(:).charge]==pepLibrarySRM(i).charge));
        
        if numel(idx_lib_cal)>1
            disp('vidimu ch''ami''i''fa')
            idx_lib_cal = idx_lib_cal(1); % it shouldn't happen but, if it does, I expect them to be the same peptide analysed more times
        elseif isempty(idx_lib_cal)
            disp(['Peptide ' pepLibrarySRM(i).sequence ' was not analysed in the calibration experiment'])
            pepLibrarySRM(i).estimatedAbundance(j)=NaN;
            if settings.mProphet
                for numMethod=1:numel(mProphet_methods)
                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                end
            end
        else
            
            if ~isempty(pepLibrarySRMCal(idx_lib_cal).regression) && (nanmedian([pepLibrarySRM(i).ratio(molarities_files(idx_mol,2))]).^(-1))~=0 && ~isinf((nanmedian([pepLibrarySRM(i).ratio(molarities_files(idx_mol,2))]).^(-1)))
                
                pepLibrarySRM(i).estimatedAbundance(j)=pepLibrarySRMCal(idx_lib_cal).regression.m*(nanmedian([pepLibrarySRM(i).ratio(molarities_files(idx_mol,2))]).^(-1))+pepLibrarySRMCal(idx_lib_cal).regression.q;
 
                if pepLibrarySRM(i).estimatedAbundance(j)<0
                    pepLibrarySRM(i).estimatedAbundance(j)=NaN;
                    disp('Abundance smaller than zero in Ariadne!!!')
                elseif isinf(pepLibrarySRM(i).estimatedAbundance(j))
                    pepLibrarySRM(i).estimatedAbundance(j) = NaN;
                end
            else
                pepLibrarySRM(i).estimatedAbundance(j)=NaN;
            end
            
            if settings.mProphet
                if ~isempty(pepLibrarySRMCal(idx_lib_cal).regression_mProphet)
                    for numMethod=1:numel(mProphet_methods)
                        %                         molarities_files(idx_mol,2)
                        if isfield(pepLibrarySRMCal(idx_lib_cal).regression_mProphet, mProphet_methods{numMethod})
                            if ~isempty(eval(['pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod}]))
                                idxs_cc_dataFiles=intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio_mProphet],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                                idxs_cc_dataFiles = intersect(idxs_cc_dataFiles,1:1:size(pepLibrarySRM(i).mScore,2));
                                idxs_cc_dataFiles = idxs_cc_dataFiles(find(pepLibrarySRM(i).mScore(idxs_cc_dataFiles)<settings.FDRthreshold));
                                
                                if ~isempty(idxs_cc_dataFiles) && eval(['(nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  ']).^(-1))~=0']) && ~isinf(eval(['(nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  ']).^(-1))']))
                                    eval([ 'pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.m*(nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  ']).^(-1)) + pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.q;'])
                                else
                                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                                end
                                if eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)<0'])
                                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                                end
                            else
                                eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                            end
                        else
                            eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                        end
                    end
                elseif isempty(pepLibrarySRMCal(idx_lib_cal).regression_mProphet)
                    for numMethod=1:numel(mProphet_methods)
                        eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                    end
                end
            end
            
        end
        clear idx_mol
    end
    %     end
end

for i=1:size(pepLibrarySRM,2)
    i
    estimatedAbundances.Ariadne(i,:)=pepLibrarySRM(i).estimatedAbundance; % after calibration values
    if settings.mProphet
        estimatedAbundances.totalxic(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_totalxic; % after calibration values
        estimatedAbundances.maxapex(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_maxapex; % after calibration values
        estimatedAbundances.apexsum(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum; % after calibration values
        estimatedAbundances.apexsum_outlier(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum_outlier; % after calibration values
    end
end

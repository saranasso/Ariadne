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
function [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExpWithSameBackground(pepLibrarySRM,association_file,pepLibrarySRMCal)
% it computes the absolute quantifications starting from the calibration
% curves estimated in a previous experiment and when you spike in light reference peptides you might have in the background.

global settings
% mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier';'Skyline_mProphet'};
mProphet_methods={'Skyline_mProphet'};
if isfield(pepLibrarySRM,'backCorrectedAbEst')
    pepLibrarySRM = rmfield(pepLibrarySRM,'backCorrectedAbEst')
    pepLibrarySRM = rmfield(pepLibrarySRM,'estimatedAbundance')
    for i=1:numel(mProphet_methods)
        eval(['pepLibrarySRM = rmfield(pepLibrarySRM,''estimatedAbundance_mProphet_' mProphet_methods{i} ''')'])
        eval(['pepLibrarySRM = rmfield(pepLibrarySRM,''backCorrectedAbEst_' mProphet_methods{i} ''')'])
    end
end
for i=1: numel(pepLibrarySRM)
    
    [idx_bool idx_ht]=ismember([association_file(:,2)],[pepLibrarySRM(i).dataFile]);
    groups_files=[cell2mat(association_file(idx_bool,1)) idx_ht(idx_bool)];
    groups=sort(unique([association_file{:,1}]),'descend');
    for j=1:numel(groups)
        idx_mol=groups_files(:,1)==groups(j);
        idx_lib_cal =  intersect(find(ismember({pepLibrarySRMCal(:).sequence},pepLibrarySRM(i).sequence)==1),find([pepLibrarySRMCal(:).charge]==pepLibrarySRM(i).charge));
        
        if numel(idx_lib_cal)>1
            disp(['Current peptide shows multiple calibrations: first one is going to be used'])
            disp(['Here the sequence ' pepLibrarySRM(i).sequence ' and the charge state ' num2str(pepLibrarySRM(i).charge)])
            idx_lib_cal = idx_lib_cal(1); % it shouldn't happen but, if it does, I expect them to be the same peptide (with same charge state!) analysed more times
        elseif isempty(idx_lib_cal)
            disp(['Current peptide was not analysed in the calibration experiment'])
            disp(['Here the sequence ' pepLibrarySRM(i).sequence ' and the charge state ' num2str(pepLibrarySRM(i).charge)])
            pepLibrarySRM(i).estimatedAbundance(j)=NaN;
            if settings.mProphet
                for numMethod=1:numel(mProphet_methods)
                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                end
            end
        else
            
            if ~isempty(pepLibrarySRMCal(idx_lib_cal).regression)
                if strcmp(settings.constantTarget, 'light')
                    idxs_cc_dataFiles=intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                    pepLibrarySRM(i).estimatedAbundance(j) = pepLibrarySRMCal(idx_lib_cal).regression.m*nanmedian(([pepLibrarySRM(i).ratio(groups_files(idx_mol,2))].^(-1)))+pepLibrarySRMCal(idx_lib_cal).regression.q;
                    pepLibrarySRM(i).backgroundCorrectionFactor(j) = ((groups(j) - pepLibrarySRMCal(idx_lib_cal).regression.q)/pepLibrarySRMCal(idx_lib_cal).regression.m)/(nanmedian(([pepLibrarySRM(i).ratio(groups_files(idx_mol,2))])).^(-1));

                    for k = (setxor(j,1:1:numel(groups)))
                           idx_mol_correction=groups_files(:,1)==groups((k));
                            value = [pepLibrarySRMCal(idx_lib_cal).regression.m*((nanmedian(([pepLibrarySRM(i).ratio(groups_files(idx_mol_correction,2))])).^(-1))*pepLibrarySRM(i).backgroundCorrectionFactor(j))+pepLibrarySRMCal(idx_lib_cal).regression.q];
                            if isinf(value) || value<=0 || pepLibrarySRM(i).backgroundCorrectionFactor(j)==0 %|| k ~= min(setxor(j,1:1:numel(groups)))
                                value = NaN;
                            end
                            if isfield(pepLibrarySRM(i),'backCorrectedAbEst')
                                if k<=numel(pepLibrarySRM(i).backCorrectedAbEst)
                                    pepLibrarySRM(i).backCorrectedAbEst{k} =[pepLibrarySRM(i).backCorrectedAbEst{k}, value];
                                else
                                    pepLibrarySRM(i).backCorrectedAbEst{k} =value;
                                end
                            else
                                pepLibrarySRM(i).backCorrectedAbEst{k} = value;
                            end
                            clear idx_mol_correction
                    end
                    
                      
                else strcmp(settings.constantTarget, 'heavy')
                    pepLibrarySRM(i).estimatedAbundance(j)=pepLibrarySRMCal(idx_lib_cal).regression.m*nanmedian([pepLibrarySRM(i).ratio(groups_files(idx_mol,2))])+pepLibrarySRMCal(idx_lib_cal).regression.q;
                end
                if pepLibrarySRM(i).estimatedAbundance(j)<0
                    pepLibrarySRM(i).estimatedAbundance(j)=NaN;
                    disp('Abundance smaller than zero in Ariadne!')
                elseif isinf(pepLibrarySRM(i).estimatedAbundance(j))
                    pepLibrarySRM(i).estimatedAbundance(j) = NaN;
                end
            else
                pepLibrarySRM(i).estimatedAbundance(j)=NaN;
            end
            
            if settings.mProphet
                if ~isempty(pepLibrarySRMCal(idx_lib_cal).regression_mProphet)
                    for numMethod=1:numel(mProphet_methods)
                        if isfield(pepLibrarySRMCal(idx_lib_cal).regression_mProphet, mProphet_methods{numMethod})
                            if ~isempty(eval(['pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod}]))
                                idxs_cc_dataFiles=intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio_mProphet],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                                idxs_cc_dataFiles = intersect(idxs_cc_dataFiles,1:1:size(pepLibrarySRM(i).mScore,2));
                                idxs_cc_dataFiles = idxs_cc_dataFiles(find(pepLibrarySRM(i).mScore(idxs_cc_dataFiles)<settings.FDRthreshold));
                                
                                if isfield(pepLibrarySRM(i).ratio_mProphet, mProphet_methods{numMethod}) && ~isempty(idxs_cc_dataFiles) %&& boolVar
                                    if strcmp(settings.constantTarget, 'light')
                                        eval([ 'pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.m*nanmedian(([pepLibrarySRM(i).ratio_mProphet(groups_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  '].^(-1))) + pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.q;'])
                                        eval([ 'pepLibrarySRM(i).backgroundCorrectionFactor_' mProphet_methods{numMethod} '(j) = ((groups(j) - pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.q)/pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.m)/(nanmedian(([pepLibrarySRM(i).ratio_mProphet(groups_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  '])).^(-1)) ;'])
                                        
                                                                                
                                        for k = (setxor(j,1:1:numel(groups)))
                                            idx_mol_correction=groups_files(:,1)==groups((k));
                                            idx_mol_correction=intersect(find(idx_mol_correction==1),1:1:size([pepLibrarySRM(i).ratio_mProphet],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                                            idx_mol_correction = intersect(idx_mol_correction,1:1:size(pepLibrarySRM(i).mScore,2));
                                            idx_mol_correction = idx_mol_correction(find(pepLibrarySRM(i).mScore(idx_mol_correction)<settings.FDRthreshold));
                                            
                                            
                                            value = eval(['pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.m*((nanmedian(([pepLibrarySRM(i).ratio_mProphet(groups_files(idx_mol_correction,2)).' mProphet_methods{numMethod}  '])).^(-1))*pepLibrarySRM(i).backgroundCorrectionFactor_' mProphet_methods{numMethod} '(j))+pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.q;']);

                                            if isinf(value) || value<=0 ||  eval([ 'pepLibrarySRM(i).backgroundCorrectionFactor_' mProphet_methods{numMethod} '(j)==0'])%|| k ~= min(setxor(j,1:1:numel(groups)))
                                                value = NaN;
                                            end
                                            if isfield(pepLibrarySRM(i),['backCorrectedAbEst_' mProphet_methods{numMethod}])
                                                if k<=numel(eval(['pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod}]))
                                                    eval(['pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} '{k} = [pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} '{k}, value];'])
                                                else
                                                    eval(['pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} '{k} = value;'])
                                                end
                                            else
                                                eval(['pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} '{k} = value;'])
                                            end
                                            clear idx_mol_correction
                                        end
                                        
                                        
                                        pause(0.01)
                                    else strcmp(settings.constantTarget, 'heavy')
                                        eval([ 'pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.m*nanmedian([pepLibrarySRM(i).ratio_mProphet(groups_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  ']) + pepLibrarySRMCal(idx_lib_cal).regression_mProphet.' mProphet_methods{numMethod} '.q;'])
                                    end
                                else
                                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                                end
                                boolVar = 0;
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
end

for i=1: numel(pepLibrarySRM)
    for j=1:numel(groups)
            if isfield(pepLibrarySRM(i),'backCorrectedAbEst') && ~isempty(pepLibrarySRM(i).backCorrectedAbEst)
                pepLibrarySRM(i).estimatedAbundance(j) = nanmean(pepLibrarySRM(i).backCorrectedAbEst{j});
            else
                pepLibrarySRM(i).estimatedAbundance(j) = nan;
            end

    end
end

for i=1: numel(pepLibrarySRM)
    for j=1:numel(groups)
        for numMethod=1:numel(mProphet_methods)
            
                if isfield(pepLibrarySRM(i),['backCorrectedAbEst_' mProphet_methods{numMethod}]) && eval(['~isempty(pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} ')'])
                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j) = nanmean(pepLibrarySRM(i).backCorrectedAbEst_' mProphet_methods{numMethod} '{j});'])
                else
                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j) = nan;'])
                end

            
        end
    end
end



for i=1:size(pepLibrarySRM,2) % absolute quantifications
    i
    estimatedAbundances.Ariadne(i,:)=pepLibrarySRM(i).estimatedAbundance;
    if settings.mProphet
        
%         estimatedAbundances.totalxic(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_totalxic;
%         estimatedAbundances.maxapex(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_maxapex;
%         estimatedAbundances.apexsum(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum;
%         estimatedAbundances.apexsum_outlier(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum_outlier;
%         estimatedAbundances.Skyline_mProphet(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_Skyline_mProphet; % after calibration values
        for numMethod=1:numel(mProphet_methods)
            eval(['estimatedAbundances.' mProphet_methods{numMethod} '(i,:) = pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} ';'])
        end
    end
end

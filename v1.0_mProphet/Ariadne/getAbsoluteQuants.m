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
function [estimatedAbundances pepLibrarySRM]=getAbsoluteQuants(pepLibrarySRM,association_file)
% it computes absolute abundances

global settings
mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};

for i=1: numel(pepLibrarySRM)
    [idx_bool idx_ht]=ismember([association_file(:,2)],[pepLibrarySRM(i).dataFile]);
    molarities_files=[cell2mat(association_file(idx_bool,1)) idx_ht(idx_bool)];
    molarities=sort(unique([association_file{:,1}]),'descend');
    for j=1:numel(molarities)
        idx_mol=molarities_files(:,1)==molarities(j);
        
        if ~isempty(pepLibrarySRM(i).regression)
            %             try
            idxs_cc_dataFiles = intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio],2)); % if some data files were not quantified then the variable is truncated at a shorter position
            if strcmp(settings.constantTarget, 'light')
                pepLibrarySRM(i).estimatedAbundance(j)=pepLibrarySRM(i).regression.m*nanmedian([pepLibrarySRM(i).ratio(molarities_files(idxs_cc_dataFiles,2))].^(-1))+pepLibrarySRM(i).regression.q;
            else strcmp(settings.constantTarget, 'heavy')
                pepLibrarySRM(i).estimatedAbundance(j)=pepLibrarySRM(i).regression.m*nanmedian([pepLibrarySRM(i).ratio(molarities_files(idxs_cc_dataFiles,2))])+pepLibrarySRM(i).regression.q;
            end
            %             catch
            %                 pepLibrarySRM(i).estimatedAbundance(j)=NaN;
            %             end
            if pepLibrarySRM(i).estimatedAbundance(j)<0
                pepLibrarySRM(i).estimatedAbundance(j)=NaN;
                disp('Abundance smaller than zero in Ariadne!!!')
            end
        else
            pepLibrarySRM(i).estimatedAbundance(j)=NaN;
        end
        
        if settings.mProphet
            if ~isempty(pepLibrarySRM(i).regression_mProphet)
                for numMethod=1:numel(mProphet_methods)
                    %                 try
                    if isfield(pepLibrarySRM(i).regression_mProphet, mProphet_methods{numMethod})
                        if ~isempty(eval(['pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod}]))
                            idxs_cc_dataFiles = intersect(find(idx_mol==1),1:1:size([pepLibrarySRM(i).ratio_mProphet],2)); % if some data files were not quantified then the variable is truncated at a shorter position
                            idxs_cc_dataFiles = intersect(idxs_cc_dataFiles,1:1:size(pepLibrarySRM(i).mScore,2));
                            idxs_cc_dataFiles = idxs_cc_dataFiles(find(pepLibrarySRM(i).mScore(idxs_cc_dataFiles)<settings.FDRthreshold));
                            if strcmp(settings.constantTarget, 'light')
                                eval([ 'pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod} '.m*nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  '].^(-1))+pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod} '.q;'])
                            else strcmp(settings.constantTarget, 'heavy')
                                eval([ 'pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod} '.m*nanmedian([pepLibrarySRM(i).ratio_mProphet(molarities_files(idxs_cc_dataFiles,2)).' mProphet_methods{numMethod}  '])+pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod} '.q;'])
                            end
                            if eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)<0'])
                                eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                            end
                        else
                            eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                        end
                    else
                        %                 catch
                        %                         eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                        %                     end
                        eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                    end
                end
            elseif isempty(pepLibrarySRM(i).regression_mProphet)
                for numMethod=1:numel(mProphet_methods)
                    eval(['pepLibrarySRM(i).estimatedAbundance_mProphet_' mProphet_methods{numMethod} '(j)=NaN;'])
                end
            end
        end
        clear idx_mol
    end
end

for i=1:size(pepLibrarySRM,2) % after calibration values
    estimatedAbundances.Ariadne(i,:)=pepLibrarySRM(i).estimatedAbundance; 
    estimatedAbundances.totalxic(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_totalxic; 
    estimatedAbundances.maxapex(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_maxapex; 
    estimatedAbundances.apexsum(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum; 
    estimatedAbundances.apexsum_outlier(i,:)=pepLibrarySRM(i).estimatedAbundance_mProphet_apexsum_outlier;
end

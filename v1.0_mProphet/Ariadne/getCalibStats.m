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
function [median_R2 median_RMSE median_localDynRange]=getCalibStats(pepLibrarySRM,method)
% it computes stats associated to the calibration curves

switch method
    
    case 'Ariadne'
        for i=1:numel(pepLibrarySRM)
            if ~isempty(pepLibrarySRM(i).regression)
                [localDynamicRange(i), RMSE(i), R2(i), minEstimate(i), maxEstimate(i)] = getCalibInfo4Stats(pepLibrarySRM(i).regression);
            end
        end
        median_localDynRange=nanmedian(localDynamicRange);
        median_RMSE=nanmedian(RMSE);
        median_R2=nanmedian(R2);
        
    case 'mProphet'
        mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};
        for numMethod=1:numel(mProphet_methods)
            for i=1:numel(pepLibrarySRM)
                if isfield(pepLibrarySRM(i).regression_mProphet, mProphet_methods{numMethod})
                    if ~isempty(eval(['pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod}]))
                        [localDynamicRange(i), RMSE(i), R2(i), minEstimate(i), maxEstimate(i)] = getCalibInfo4Stats(eval(['pepLibrarySRM(i).regression_mProphet.' mProphet_methods{numMethod}]));
                    end
                end
            end
            eval(['median_localDynRange.' mProphet_methods{numMethod} '=nanmedian(localDynamicRange);']);
            eval(['median_RMSE.' mProphet_methods{numMethod} '=nanmedian(RMSE);']);
            eval(['median_R2.' mProphet_methods{numMethod} '=nanmedian(R2);']);
        end
end	

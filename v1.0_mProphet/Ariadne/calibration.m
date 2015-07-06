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

%% it estimates calibration curves and some quality metrics

pepLibrarySRM(find([pepLibrarySRM(:).hasRef]==0))=[];

%% possible concentrations
try
    association_file=association_dil_file;
end
molarities=sort(unique([association_file{:,1}]),'descend');
molarities(molarities>1)=round(molarities(molarities>1));

%% estimating calibration curves
for i=1:numel(pepLibrarySRM)
    
    h=figure;
    disp('Ariadne')
    [pepLibrarySRM]=getCalibration(pepLibrarySRM,association_file,i);
    title([{['Peptide # ' num2str(i) '  ' pepLibrarySRM(i).sequence]};{'Ariadne'}])
    attachPlotToFile(h,'regressions')
    
    if settings.mProphet
        h=figure;
        disp('mProphet')
        [pepLibrarySRM]=getCalibration_mProphet(pepLibrarySRM,association_file,i);
%         title([{['Peptide # ' num2str(i) '  ' pepLibrarySRM(i).sequence]};{'mProphet'}])
        attachPlotToFile(h,'regressions')
    end
    close all
end
pepLibrarySRMCal = pepLibrarySRM;
save('pepLibCalExp','pepLibrarySRMCal')

%% dynamic Range
[median_R2_Ariadne median_RMSE_Ariadne median_localDynRange_Ariadne]=getCalibStats(pepLibrarySRM,'Ariadne')
[median_R2_mProphet median_RMSE_mProphet median_localDynRange_mProphet]=getCalibStats(pepLibrarySRM,'mProphet')

%% diary off
diary off









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

% this scripts allows evaluating performance in presence of a dilution
% series and a calibration experiment.
clear
clc

dataDir = uigetdir(pwd, 'Select the Ariadne code directory');
addpath dataDir

% prompt = {['Is there any endogenous peptide in the background that might interfere with your spiked-in reference?']};
% dlg_title = 'Input boolean value';
% num_lines = 1;
% def = {'yes = 1, no =0'};
% interferingBackground = getBooleanInput(prompt,dlg_title,num_lines,def);
interferingBackground = 0; % see considerations in the paper (Performance assessment section)

[resultFile pathResultFile]=uigetfile('*.mat', 'Select the result file .mat','MultiSelect', 'Off');
resultFileFilePath=[pathResultFile resultFile];
load(resultFileFilePath)
pepLibrarySRM(find([pepLibrarySRM(:).hasRef]==0))=[];

cd(pathResultFile)
[settings pathsettings]=uigetfile('*.m', 'Select the  file configSettings.m for this dataset','MultiSelect', 'Off');
settingsFilePath=[pathsettings settings];
run(settingsFilePath)

compare2Skyline 


[calExp pathCalExp]=uigetfile('*.mat', 'Select the calibration experiment file pepLibCalExp.mat','MultiSelect', 'Off');
calExpFilePath=[pathCalExp calExp];
load(calExpFilePath)

cd(pathResultFile)
run(settingsFilePath)

mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier';'Skyline_mProphet'};
if settings.multFactor~= 1
    for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).ratio = pepLibrarySRM(i).ratio.*settings.multFactor;
        for j =1:numel(pepLibrarySRM(i).ratio_mProphet)
            for k=1:numel(mProphet_methods)
                if isfield(pepLibrarySRM(i).ratio_mProphet(j), mProphet_methods{k})
                    eval([ 'pepLibrarySRM(i).ratio_mProphet(j).' mProphet_methods{k} '= pepLibrarySRM(i).ratio_mProphet(j).' mProphet_methods{k} '*settings.multFactor;'])
                end
            end
        end
    end
end
association_file=association_dil_file

cd(pathResultFile)

if ~interferingBackground
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExp(pepLibrarySRM,association_dil_file,pepLibrarySRMCal)
else
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExpWithSameBackground(pepLibrarySRM,association_file,pepLibrarySRMCal)
end


pause off
getResultsPaperWithSkyline
pause on
break
stats










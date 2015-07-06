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

% this script imports the Skyline results for comparison aims

fileSkylineOutput=dir('*.tsv')
[SkylineOutputFile pathSkylineOutputFile]=uigetfile('*.tsv', 'Select the Skyline output file','MultiSelect', 'Off');
SkylineOutputFilePath=[pathSkylineOutputFile SkylineOutputFile];
[SkylineOutput falseHitsSkyline]=filter_SkylineOutput(SkylineOutputFilePath,settings.FDRthreshold);
pepLibrarySRM = rmfield(pepLibrarySRM, 'mScore'); 
[ pepLibrarySRM pepInListNotIn_Skyline]=getSkylineResults(pepLibrarySRM,transList,SkylineOutput,falseHitsSkyline);

if settings.calibration
    [settings pathsettings]=uigetfile('*.m', 'Select the  file configSettings.m for this dataset','MultiSelect', 'Off');
    settingsFilePath=[pathsettings settings];
    run(settingsFilePath)
    calibration
    save('pepLibCalExpSky.mat')
    break
end
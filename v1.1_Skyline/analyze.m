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
for k=1:numel(pepLibrarySRM)
    
    if isfield(pepLibrarySRM,'analyzed')
        if ~pepLibrarySRM(k).analyzed
            boolAnalyzeIt = 1;
        else
            boolAnalyzeIt = 0;
        end
    else
        boolAnalyzeIt = 1;
    end
    
    if boolAnalyzeIt
        pepLibrarySRM(k).analyzed =1;
        pepLibrarySRM(k).sequence
        if pepLibrarySRM(k).hasRef
            precursor.endo=pepLibrarySRM(k).precursor;
            
            if ~settings.mProphetFilter
                dataFiles=association_file(:,2);
            else
                dataFiles=pepLibrarySRM(k).dataFile;
            end
            
            for i=1:numel(pepLibrarySRM(k).dataFile)
                disp(['k is ' num2str(k) ' and i is ' num2str(i)])
                try
                    eval(['srm=' genvarname(dataFiles{i}) ';']);
                catch
                    disp('DataFile null')
                    pepLibrarySRM(k).ratio(i)=NaN;
                    continue
                end
                %% data retrieval -begin
                try % getPrecursorData where convenient checkings are executed
                    [idxInData, Q3s.endo, rts.endo, trans.endo]=getPrecursorData(srm, precursor.endo, unique(sort(pepLibrarySRM(k).Q3s)));
                catch
                    disp('Q3 from pepLibrarySRM and those from srm.data differ more than settings.ionsResolution')
                    pause(1)
                    pepLibrarySRM(k).ratio(i)=NaN;
                    continue
                end
                
                precursor.ref=pepLibrarySRM(k).precursorRef;
                try
                    [idxInData_ref Q3s.ref, rts.ref, trans.ref]=getPrecursorData(srm, precursor.ref, unique(sort(pepLibrarySRM(k).Q3sRef)));
                catch
                    disp('Q3 from pepLibrarySRM and those from srm.data differ more than settings.ionsResolution')
                    pause(1)
                    continue
                end
                %% data retrieval -end
                %% pep abundances and ratio computation

                [pepLibrarySRM]=ariadneSthread(precursor, Q3s, rts,trans,pepLibrarySRM,k,i);
      
            end
            close all
            java.lang.System.gc()
        elseif ~pepLibrarySRM(k).hasRef
            pepLibrarySRM(k).ratio=nan(1,numel(pepLibrarySRM(k).dataFile));
        end
        clear Q3s rts trans idxInData idxInData_ref
        pause(0.0001)
        save(['synch_lowess_sg_quantili_001099_merge_minCorr' num2str(settings.minXcorr) '_1PD_' num2str(settings.expPeakWidth) '_FDR_' num2str(settings.FDRthreshold) '.mat'])
        pause(0.0001)
    else
        continue
    end
end

if settings.calibration
    calibration
    save('pepLibCalExp.mat')
end

try
    association_file = association_dil_file; % backward compatibility
end

if ~settings.calibration && settings.absoluteQuantification
    try
        [estimatedAbundances pepLibrarySRM] = getAbsQuantsFromCalExp(pepLibrarySRM, association_file, pepLibrarySRMCal)
    catch ME
        report = getReport(ME)
        [calExp pathCalExp]=uigetfile('*.mat', 'Select the  file pepLibCalExp.mat','MultiSelect', 'Off');
        calExpFilePath=[pathCalExp calExp];
        load(calExpFilePath)
        [estimatedAbundances pepLibrarySRM] = getAbsQuantsFromCalExp(pepLibrarySRM, association_file, pepLibrarySRMCal)
    end
end

if settings.dilution_series && settings.absoluteQuantification
    % test performance
    pause off
    assessQuantQuality
    pause on
end

save(['analysisResults_minTransCorr' num2str(settings.minXcorr)  '.mat'])

%% export relative and absolute abundances


% matlabpool close

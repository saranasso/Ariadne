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

clear
close all
clc
try
    matlabpool open
catch
    disp('Warning: could not run the parallel computing option!')
end

%% adding Ariadne path to Matlab path
dataDir = uigetdir(pwd, 'Select the Ariadne code directory');
addpath dataDir
 
%% config
[settings pathsettings]=uigetfile('*.m', 'Select the  file configSettings.m for this dataset','MultiSelect', 'Off');
settingsFilePath=[pathsettings settings];
run(settingsFilePath)

%% loading calibration experiment results (pepLibCalExp.mat) when needed

if settings.absoluteQuantification && ~settings.calibration
    [calExp pathCalExp]=uigetfile('*.mat', 'Select the  file pepLibCalExp.mat','MultiSelect', 'Off');
    calExpFilePath=[pathCalExp calExp];
    load(calExpFilePath)
end

run(settingsFilePath)

%% importing mProphet results
if settings.mProphet
    filemProphetOutput=dir('*mProphet_all_peakgroups.xls')
    [mProphetOutputFile pathmProphetOutputFile]=uigetfile('*.xls', 'Select the mProphet output file','MultiSelect', 'Off');
    mProphetOutputFilePath=[pathmProphetOutputFile mProphetOutputFile];
    [mProphetOutput falseHitsmProphet]=filter_mProphetOutput(mProphetOutputFilePath,settings.FDRthreshold); % removing mProphet results without metadata. FDR filtering takes place if, and only if, settings.mProphetFilter=1
end


%% importing transition lists
[fileTransList pathTransList]=uigetfile('*.xls', 'Select the transitions list file','MultiSelect', 'On');
if iscell(fileTransList)
    for i=1:numel(fileTransList)
        fileTransList_tmp{i}=[pathTransList fileTransList{i}];
    end
else
    fileTransList_tmp{1}=[pathTransList fileTransList];
end
fileTransList=fileTransList_tmp;

%% Data import or load
mzXMLFilesPath_tmp=[];
choice=questdlg('Are you importing new data in Matlab?','Data importing...');

if isequal(choice,'Yes')
    dataDir = uigetdir(pwd, 'Select the data files directory');
    %% importing the data
    cd(dataDir)
    fileMzXML=struct2cell(dir('*.mzXML'));
    mzXMLFilesPath=fileMzXML(1,:);
    for i=1:numel(mzXMLFilesPath)
        fileNames{i}=strtok(mzXMLFilesPath{i},'.');
    end
    association_file = getAssociation(fileNames);
    for i=1:numel(mzXMLFilesPath)
        disp([ 'i is ' num2str(i) ' out of ' num2str(numel(mzXMLFilesPath))])
        filePath=[pwd '/' mzXMLFilesPath{i}];
        dataFromXML{i}=mzxmlread(filePath);
        pause(0.01)
        dataFromXML{i}=dataReFormat(filePath,dataFromXML{i});
        pause(0.01)
    end
    for i=1:numel(mzXMLFilesPath)
        filePath=[pwd '/' mzXMLFilesPath{i}];
        varName_tmp=strtok(mzXMLFilesPath{i},'.');
        eval([genvarname(varName_tmp) '=dataFromXML{i};'])
        dataFromXML{i}=[];
    end
        % association file
    mzXMLFilesPath=who('x*')';
    for i=1:numel(mzXMLFilesPath)
        mzXMLFilesPath{i}=mzXMLFilesPath{i}(2:end);
    end
    association_file = getAssociation(mzXMLFilesPath);
    save([dataDir filesep 'dataRead.mat'])
else
    flag='No';
    i=0;
    while strmatch('No', flag)
        i=i+1;
        [file{i} path{i}]=uigetfile('*.mat', 'Select the matlab file dataRead.mat storing the data','MultiSelect', 'On');
        flag=questdlg('Are you done?','Data importing...');
    end
    for i=1:numel(file)
        currentDataFile=[path{i} file{i}];
        if isstr(currentDataFile)
            load(currentDataFile)
        end
    end
    try
        association_file =association_dil_file; % backward compatibility
    end
%     % association file
%     mzXMLFilesPath=who('x*')';
%     for i=1:numel(mzXMLFilesPath)
%         mzXMLFilesPath{i}=mzXMLFilesPath{i}(2:end);
%     end
%     association_file = getAssociation(mzXMLFilesPath);
end

run(settingsFilePath)

%%  reading and filtering trans lists
[transList decoysTransList]=filter_TransList(fileTransList);

%% retrieving useful metadata to access data only from the trans list
if ~settings.mProphet && ~settings.mProphetFilter
    [pepLibrarySRM]=getPepLibFromTransList(transList)
    for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).dataFile=association_file(:,2);
    end
end

%% reading and filtering the mProphet output results
% matching peptides in trans list to true hits from mProphet results and retrieving useful metadata to access data
if settings.mProphet && settings.mProphetFilter
    [pepLibrarySRM pepInListNotIn_mProphet]=matchTransListTomProphetOut(transList,mProphetOutput,falseHitsmProphet)
end

%% retrieving useful metadata to access data from the trans list and from mProphet output
if settings.mProphet && ~settings.mProphetFilter
    [pepLibrarySRM]=getPepLibFromTransList(transList)
    for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).dataFile=association_file(:,2);
    end
    [ pepLibrarySRM pepInListNotIn_mProphet]=getmProphetResults(pepLibrarySRM,transList,mProphetOutput,falseHitsmProphet)
end

pepLibrarySRM(find([pepLibrarySRM(:).hasRef]==0))=[]

save ready4analysis

%% start analysis
analyze

%% export results to text file
exportResults2Txt

save done

return

%% ADD - ON

%% Update pepLibrary if you want to quantify additional peps/transitions you may have added to your transition list. 
%  Do not forget relative intensities! After this, run again analysis (i.e., analyze script), it
%  will attempt to quantify those without a quantification value
% Note: you only need to evaluate these rows of code

[pepLibrarySRM]=updatePepLibFromTransList(pepLibrarySRM,transList)
for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).dataFile=association_file(:,2);
end

if settings.mProphet && ~settings.mProphetFilter
    [ pepLibrarySRM pepInListNotIn_mProphet]=getmProphetResults(pepLibrarySRM,transList,mProphetOutput,falseHitsmProphet)
end

if settings.mProphet && settings.mProphetFilter
    [pepLibrarySRM pepInListNotIn_mProphet]=matchTransListTomProphetOut(transList,mProphetOutput,falseHitsmProphet)
end
pepLibrarySRM(find([pepLibrarySRM(:).hasRef]==0))=[]



	

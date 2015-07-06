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
function [retTime, maxApex, estimatedPeakWidth]=peakGroupSelection(rts_endo,endo,rts_ref,ref,pepLibraryInfo)
% it selects the best peak group for the monitored peptide target

global settings

plotBool='false';
numOfTrans=size(ref,2);

%% reducing data points number for improving  speed
% % maxCountsPerTrans = 100000; % enough for an EM estimate of the GMM - risk to be biased to intense peakgroups
% % [endo, ref, cut] = intCut (endo, ref, 'maxCountsPerTrans',maxCountsPerTrans);
% % [endo, ref, cut] = intCut (endo, ref, 'minCountIntensity',settings.minExpectPeakHeight);
[endo, ref, cut] = intCut (endo, ref, 'estimateConservative',settings.noiseLevel);
currentVars=whos;
if ~ismember({'cut'},{currentVars.name}) % if the above line is commented
    cut=1;
end

%% getting needed data

rti_min=min(rts_endo{1}(1),rts_ref{1}(1));
rtf_max=max(rts_endo{1}(end),rts_ref{1}(end));
newRetTimes=linspace(rti_min,rtf_max,(rtf_max-rti_min))';
transSumData.endo=[];
transSumData.ref=[];
transMergeEndoRef=[];
for i=1:numOfTrans
    try
        transSumData.endo=cat(1,transSumData.endo,hist2data(endo{i},rts_endo{i}));
        transSumData.ref=cat(1,transSumData.ref,hist2data(ref{i},rts_ref{i}));
        transMergeEndoRef{i}=hist(cat(1,hist2data(endo{i},rts_endo{i}),hist2data(ref{i},rts_ref{i})),newRetTimes)';
    catch 
        disp(['No data for trans #' num2str(i)]) % it should almost never happen...
    end
end

mergedTransSumData=cat(1,transSumData.endo,transSumData.ref);
mergedTrans=hist(mergedTransSumData,newRetTimes)';


%% peak picking on merged trans
mergedTrans = mergedTrans.*cut; % scales back to "original" the counts
[Peaklist, PFWHH, PExt] =  mspeaks(newRetTimes,mergedTrans,'HeightFilter',settings.minExpectPeakHeight,'OverSegmentationFilter', settings.expPeakWidth,'ShowPlot', plotBool,'Style','extline');
expPeakWidth =  settings.expPeakWidth;
while isempty(Peaklist) &&  expPeakWidth>1
    expPeakWidth=expPeakWidth/2
    [Peaklist, PFWHH, PExt] =  mspeaks(newRetTimes,mergedTrans,'HeightFilter',settings.minExpectPeakHeight,'OverSegmentationFilter', expPeakWidth,'ShowPlot', plotBool,'Style','extline');
end

if isempty(Peaklist)
    retTime = [];
    maxApex = [];
    estimatedPeakWidth = [];
    return
end

%% peakgroup model 
prior=getGMMPrior(pepLibraryInfo,settings,numOfTrans,'merge');
[retTime idxPG gm.merge  X minDist.merge]=pickPeakGroup(Peaklist,PExt,prior,transMergeEndoRef,newRetTimes,pepLibraryInfo);
prior=getGMMPrior(pepLibraryInfo,settings,numOfTrans,'endo');
[retTime_endo idxPG_endo gm.endo  X minDist.endo]=pickPeakGroup(Peaklist,PExt,prior,endo,rts_endo{1},pepLibraryInfo);
prior=getGMMPrior(pepLibraryInfo,settings,numOfTrans,'ref');
[retTime_ref idxPG_ref gm.ref  X minDist.ref]=pickPeakGroup(Peaklist,PExt,prior,ref,rts_ref{1},pepLibraryInfo);

idx=mode([idxPG,idxPG_endo,idxPG_ref]);
retTime=Peaklist(idx,1);
rti=PExt(idx,1);
rtf=PExt(idx,2);
rts.endo=rts_endo{1};
rts.ref=rts_ref{1};
[rti_endo rti_ref rtf_endo rtf_ref]=getRange(rti,rtf,rts);
idxs_endo=find(rts.endo>rti_endo & rts.endo<rtf_endo);
idxs_ref=find(rts.ref>rti_ref & rts.ref<rtf_ref);

for i=1:numOfTrans
    maxApex.endo(i)=max(endo{i}(idxs_endo))*cut;
    maxApex.ref(i)=max(ref{i}(idxs_ref))*cut;
end
gm = gm.merge;
estimatedPeakWidth=max([min(gm.Sigma(1,1,:)), settings.expPeakWidth]);

% wDistance = 
	

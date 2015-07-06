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
function [ pepLibrary pepInListNotIn_Skyline]=getSkylineResults(pepLibrary,transList,SkylineOutput,falseHitsSkyline)
% this function finds common peptides between trans list and Skyline output
% and creates the peptides library

global settings 

%%  finding common peptides between trans list and Skyline output
pepSeqIdxTransList=strmatch('stripped_sequence',transList(1,:),'exact');
for i=2:size(transList,1)
    pepSeqFromTransList{i-1}=[transList{i,pepSeqIdxTransList} transList{i,strmatch('prec_z',transList(1,:),'exact')}];
end
pepSeqIdx_SkylineOut=strmatch('PeptideSequence',SkylineOutput(1,:),'exact');
for i=2:size(SkylineOutput,1)
    pepSeqFromSkylineOutput{i-1}=[SkylineOutput{i,pepSeqIdx_SkylineOut} SkylineOutput{i,strmatch('PrecursorCharge',SkylineOutput(1,:),'exact')}];
end
pepInListNotIn_Skyline=setdiff(pepSeqFromTransList,pepSeqFromSkylineOutput);
for i=1:size(pepLibrary,2)
    pepSeqInPepLib{i}=[pepLibrary(i).sequence num2str(pepLibrary(i).charge)];
end
% let's try to see why they're not in the Skyline output results
if ~isempty(pepInListNotIn_Skyline)
    for i=1:size(pepInListNotIn_Skyline,1)
        currPep=pepInListNotIn_Skyline{i};
        if ~isempty(falseHitsSkyline(:,strmatch(currPep,SkylineOutput(1,:),'exact')))
            disp(['Peptide ' currPep ' was removed from Skyline results as false hit' ])
        elseif isempty(falseHitsSkyline(:,strmatch(currPep,SkylineOutput(1,:),'exact')))
            disp(['Peptide ' currPep ' never appeared in Skyline results' ])
        end
    end
end
for i=1:numel(pepLibrary)
    pepSeqFromPepLibrary{i}=[pepLibrary(i).sequence, num2str(pepLibrary(i).charge)];
end
[matchedPeps idxTransList idxProphOut]=intersect(pepSeqFromPepLibrary,pepSeqFromSkylineOutput);


for i=1:size(matchedPeps,2)
    currPep=[matchedPeps{i}];
    idx_SkylineOut = strmatch(currPep,pepSeqFromSkylineOutput,'exact')+1;
    idx_pepLibrary = strmatch(currPep,pepSeqInPepLib,'exact');
    
    for j=1:numel(idx_SkylineOut)
        %% info changing among different dataFiles: Tr (ret time) best peak group, datafile and m & d score
        idx=strmatch( strtok(SkylineOutput(idx_SkylineOut(j),strmatch('FileName',SkylineOutput(1,:),'exact')),'.'),pepLibrary(idx_pepLibrary).dataFile);
        if ~isempty(str2num(SkylineOutput{idx_SkylineOut(j),strmatch('annotation_QValue',SkylineOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).mScore(idx)=str2num(SkylineOutput{idx_SkylineOut(j),strmatch('annotation_QValue',SkylineOutput(1,:),'exact')});
        end
        if ~isempty(str2num(SkylineOutput{idx_SkylineOut(j),strmatch('RatioLightToHeavy',SkylineOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).ratio_mProphet(idx).Skyline_mProphet=str2num(SkylineOutput{idx_SkylineOut(j),strmatch('RatioLightToHeavy',SkylineOutput(1,:),'exact')});
        end
    end
end
pepLibrary(find([pepLibrary(:).hasRef]==0))=[];
	

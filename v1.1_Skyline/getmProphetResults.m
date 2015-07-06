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
function [ pepLibrary pepInListNotIn_mProphet]=getmProphetResults(pepLibrary,transList,mProphetOutput,falseHitsmProphet)
% this function finds common peptides between trans list and mProphet output
% and creates the peptides library

global settings 

%%  finding common peptides between trans list and mProphet output
pepSeqIdxTransList=strmatch('stripped_sequence',transList(1,:),'exact');
for i=2:size(transList,1)
    pepSeqFromTransList{i-1}=[transList{i,pepSeqIdxTransList} transList{i,strmatch('prec_z',transList(1,:),'exact')}];
end
pepSeqIdx_mProphetOut=strmatch('transition_group_pepseq',mProphetOutput(1,:),'exact');
for i=2:size(mProphetOutput,1)
    pepSeqFrommProphetOutput{i-1}=[mProphetOutput{i,pepSeqIdx_mProphetOut} mProphetOutput{i,strmatch('transition_group_charge',mProphetOutput(1,:),'exact')}];
end
pepInListNotIn_mProphet=setdiff(pepSeqFromTransList,pepSeqFrommProphetOutput);
for i=1:size(pepLibrary,2)
    pepSeqInPepLib{i}=[pepLibrary(i).sequence num2str(pepLibrary(i).charge)];
end
% let's try to see why they're not in the mProphet output results
if ~isempty(pepInListNotIn_mProphet)
    for i=1:size(pepInListNotIn_mProphet,1)
        currPep=pepInListNotIn_mProphet{i};
        if ~isempty(falseHitsmProphet(:,strmatch(currPep,mProphetOutput(1,:),'exact')))
            disp(['Peptide ' currPep ' was removed from mProphet results as false hit' ])
        elseif isempty(falseHitsmProphet(:,strmatch(currPep,mProphetOutput(1,:),'exact')))
            disp(['Peptide ' currPep ' never appeared in mProphet results' ])
        end
    end
end
for i=1:numel(pepLibrary)
    pepSeqFromPepLibrary{i}=[pepLibrary(i).sequence, num2str(pepLibrary(i).charge)];
end
[matchedPeps idxTransList idxProphOut]=intersect(pepSeqFromPepLibrary,pepSeqFrommProphetOutput);


for i=1:size(matchedPeps,2)
    currPep=[matchedPeps{i}];
    idx_mProphetOut = strmatch(currPep,pepSeqFrommProphetOutput,'exact')+1;
    idx_pepLibrary = strmatch(currPep,pepSeqInPepLib,'exact');
    
    for j=1:numel(idx_mProphetOut)
        %% info changing among different dataFiles: Tr (ret time) best peak group, datafile and m & d score
        idx=strmatch( mProphetOutput(idx_mProphetOut(j),strmatch('file_name',mProphetOutput(1,:),'exact')),pepLibrary(idx_pepLibrary).dataFile);
        if strcmp(settings.timeUnit,'sec')
            pepLibrary(idx_pepLibrary).exp_RT_BestPeakGroup(idx)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('Tr',mProphetOutput(1,:),'exact')});
        else
            pepLibrary(idx_pepLibrary).exp_RT_BestPeakGroup(idx)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('Tr_min',mProphetOutput(1,:),'exact')})*60;
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('m_score',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).mScore(idx)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('m_score',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('d_score',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).dScore(idx)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('d_score',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_totalxic',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).ratio_mProphet(idx).totalxic=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_totalxic',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_maxapex',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).ratio_mProphet(idx).maxapex=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_maxapex',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).ratio_mProphet(idx).apexsum=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum_outlier',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).ratio_mProphet(idx).apexsum_outlier=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum_outlier',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic_reference',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).abundance_mProphet(idx).total_xic_ref=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic_reference',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).abundance_mProphet(idx).total_xic_endo=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity_reference',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).abundance_mProphet(idx).max_apex_ref=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity_reference',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity',mProphetOutput(1,:),'exact')}))
            pepLibrary(idx_pepLibrary).abundance_mProphet(idx).max_apex_endo=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity',mProphetOutput(1,:),'exact')});
        end
    end
end
pepLibrary(find([pepLibrary(:).hasRef]==0))=[];
	

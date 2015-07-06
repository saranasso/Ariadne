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
function [ pepLibrary pepInListNotIn_mProphet]=matchTransListTomProphetOut(transList,mProphetOutput,falseHitsmProphet)
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
[matchedPeps idxTransList idxProphOut]=intersect(pepSeqFromTransList,pepSeqFrommProphetOutput);

%% creating the peptide library
allQ1s=transList(:,strmatch('Q1',transList(1,:),'exact'));
allQ1s=delEmptySpaces(allQ1s);

for i=1:size(matchedPeps,2)
    currPep=[matchedPeps{i}];
    idx_mProphetOut=strmatch(currPep,pepSeqFrommProphetOutput,'exact')+1;
    idxTransList=strmatch(currPep,pepSeqFromTransList,'exact')+1;
    pepLibrary(i).sequence=currPep(1:end-1);     
    pepLibrary(i).analyzed = 0;
    %% common info for all files info from mProphet output
    pepLibrary(i).protein=mProphetOutput(idx_mProphetOut(1),strmatch('protein',mProphetOutput(1,:),'exact'));
    pepLibrary(i).charge=str2num(mProphetOutput{idx_mProphetOut(1),strmatch('transition_group_charge',mProphetOutput(1,:),'exact')});
    if strcmp(settings.timeUnit,'sec')
        pepLibrary(i).RT_time=str2double(transList(idxTransList(1),strmatch('AverageMeasuredRetentionTime',transList(1,:),'exact')));
        if isempty(pepLibrary(i).RT_time)
            pepLibrary(i).RT_time=str2double(transList(idxTransList(1),strmatch('PredictedRetentionTime',transList(1,:),'exact')));
        end
    else
        pepLibrary(i).RT_time=str2double(transList(idxTransList(1),strmatch('AverageMeasuredRetentionTime',transList(1,:),'exact')))*60;
        if isempty(pepLibrary(i).RT_time)
            pepLibrary(i).RT_time=str2double(transList(idxTransList(1),strmatch('PredictedRetentionTime',transList(1,:),'exact')))*60;
        end
    end
    
    Q1s_tmp=transList(idxTransList,strmatch('Q1',transList(1,:),'exact'));
    Q1s_tmp=delEmptySpaces(Q1s_tmp);
    Q1s=unique(Q1s_tmp);
    cond=numel(Q1s)==2;
    cond1=0;
    if cond
        if ~isempty(strfind(pepLibrary(i).sequence,'K'))
            cond1=abs(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(i).charge-settings.deltaLabelK)<settings.massShiftTolerance;
            if ~cond1
                disp(['2 different Q1s for peptide ' currPep ' but label mass shift not respected: actual one in charge 1 would be ' num2str(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(i).charge) ' instead of ' num2str(settings.deltaLabelK)])
                pause(1)
            end
        elseif ~isempty(strfind(pepLibrary(i).sequence,'R'))
            cond1=abs(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(i).charge-settings.deltaLabelR)<settings.massShiftTolerance;
            if ~cond1
                disp(['2 different Q1s for peptide ' currPep ' but label mass shift not respected: actual one in charge 1 would be ' num2str(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(i).charge) ' instead of ' num2str(settings.deltaLabelR)])
                pause(1)
            end
        end
    end
    pepLibrary(i).hasRef=cond1;
   
    if pepLibrary(i).hasRef
        if ~settings.invertedLabel %heavy is the ref
            precursor(i)=(Q1s(1));
            precursorRef(i)=(Q1s(2));
        else
            precursor(i)=(Q1s(2));
            precursorRef(i)=(Q1s(1));
        end
        Q3sRef_tmp=(transList(strmatch((precursorRef(i)),allQ1s),strmatch('Q3',transList(1,:))));
        [Q3sRef_tmp idx_sort_Q3sRef]=unique(delEmptySpaces(Q3sRef_tmp));
        relIntRef_tmp_all=(transList(strmatch((precursorRef(i)),allQ1s),strmatch('relative_intensity',transList(1,:))));
        relIntRef_tmp=relIntRef_tmp_all(idx_sort_Q3sRef);
        for j=1:length(Q3sRef_tmp)
            pepLibrary(i).Q3sRef(j)=str2num(Q3sRef_tmp{j});
            pepLibrary(i).relIntRef(j)=str2num(relIntRef_tmp{j});
        end
        pepLibrary(i).RT_kit_pep=0;
    else
        precursor(i)=(Q1s(1));
        precursorRef(i)={-1};
        if strmatch('RT', pepLibrary(i).protein)
            pepLibrary(i).RT_kit_pep=1;
        else
            pepLibrary(i).RT_kit_pep=0;
        end
    end
    
    Q3s_tmp=(transList(strmatch(precursor(i),allQ1s),strmatch('Q3',transList(1,:))));
    [Q3s_tmp idx_sort_Q3]=unique(delEmptySpaces(Q3s_tmp));
    relInt_tmp_all=(transList(strmatch(precursor(i),allQ1s),strmatch('relative_intensity',transList(1,:))));
    relInt_tmp=relInt_tmp_all(idx_sort_Q3);
    for j=1:length(Q3s_tmp)
        pepLibrary(i).Q3s(j)=str2num(Q3s_tmp{j});
        pepLibrary(i).relInt(j)=str2num(relInt_tmp{j});
    end
    pepLibrary(i).precursor=str2num(precursor{i});
    if ~isequal(precursorRef{i},-1)
        pepLibrary(i).precursorRef=str2num(precursorRef{i});
    else
        pepLibrary(i).precursorRef=-1;
    end
    
    for j=1:numel(idx_mProphetOut)
        %% variables changing among different dataFiles: Tr (ret time) best peak group, datafile and m & d score
        pepLibrary(i).dataFile(j)=mProphetOutput(idx_mProphetOut(j),strmatch('file_name',mProphetOutput(1,:),'exact'));
        if strcmp(settings.timeUnit,'sec')
            pepLibrary(i).exp_RT_BestPeakGroup(j)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('Tr',mProphetOutput(1,:),'exact')});
        else
            pepLibrary(i).exp_RT_BestPeakGroup(j)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('Tr_min',mProphetOutput(1,:),'exact')})*60;
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('m_score',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).mScore(j)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('m_score',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('d_score',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).dScore(j)=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('d_score',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_totalxic',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).ratio_mProphet(j).totalxic=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_totalxic',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_maxapex',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).ratio_mProphet(j).maxapex=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_maxapex',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).ratio_mProphet(j).apexsum=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum_outlier',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).ratio_mProphet(j).apexsum_outlier=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('light_heavy_ratio_apexsum_outlier',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic_reference',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).abundance_mProphet(j).total_xic_ref=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic_reference',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).abundance_mProphet(j).total_xic_endo=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('total_xic',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity_reference',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).abundance_mProphet(j).max_apex_ref=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity_reference',mProphetOutput(1,:),'exact')});
        end
        if ~isempty(str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity',mProphetOutput(1,:),'exact')}))
            pepLibrary(i).abundance_mProphet(j).max_apex_endo=str2num(mProphetOutput{idx_mProphetOut(j),strmatch('max_apex_intensity',mProphetOutput(1,:),'exact')});
        end
    end
end
	

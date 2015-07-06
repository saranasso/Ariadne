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
function [pepLibrary]=updatePepLibFromTransList(pepLibrary,transList)
% this function finds common peptides between trans list and mProphet output
% and creates the peptides library

global settings

%%  finding common peptides between trans list and mProphet output
pepSeqIdxTransList=strmatch('stripped_sequence',transList(1,:),'exact');

for i=2:size(transList,1)
    pepSeqFromTransListAll{i-1}=[transList{i,pepSeqIdxTransList} transList{i,strmatch('prec_z',transList(1,:),'exact')}];
end
pepSeqFromTransList=unique(pepSeqFromTransListAll);
for i=1:size(pepLibrary,2)
    pepSeqInPepLib{i}=[pepLibrary(i).sequence num2str(pepLibrary(i).charge)];
end
%% creating the peptides library
allQ1s=transList(:,strmatch('Q1',transList(1,:),'exact'));
allQ1s=delEmptySpaces(allQ1s);

for i=1:size(pepSeqFromTransList,2)
    currPep=pepSeqFromTransList{i};
    idxPepLibrary = strmatch(currPep,pepSeqInPepLib);
    if isempty(idxPepLibrary);
        idxPepLibrary = size(pepLibrary,2) +1;
    end
    idxTransList=strmatch(currPep,pepSeqFromTransListAll)+1; % +1 is offset header in transList
    pepLibrary(idxPepLibrary).sequence=currPep(1:end-1);
    
    %% common info for all files info from mProphet output
    pepLibrary(idxPepLibrary).protein=transList(idxTransList(1),strmatch('protein_name',transList(1,:),'exact'));
    pepLibrary(idxPepLibrary).charge=str2num(transList{idxTransList(1),strmatch('prec_z',transList(1,:),'exact')});
    if strcmp(settings.timeUnit,'sec')
        pepLibrary(idxPepLibrary).RT_time=str2double(transList(idxTransList(1),strmatch('AverageMeasuredRetentionTime',transList(1,:),'exact')));
        if isempty(pepLibrary(idxPepLibrary).RT_time)
            pepLibrary(idxPepLibrary).RT_time=str2double(transList(idxTransList(1),strmatch('PredictedRetentionTime',transList(1,:),'exact')));
        end
    else
        pepLibrary(idxPepLibrary).RT_time=str2double(transList(idxTransList(1),strmatch('AverageMeasuredRetentionTime',transList(1,:),'exact')))*60;
        if isempty(pepLibrary(idxPepLibrary).RT_time)
            pepLibrary(idxPepLibrary).RT_time=str2double(transList(idxTransList(1),strmatch('PredictedRetentionTime',transList(1,:),'exact')))*60;
        end
    end
    
    Q1s_tmp=transList(idxTransList,strmatch('Q1',transList(1,:),'exact'));
    Q1s_tmp=delEmptySpaces(Q1s_tmp);
    Q1s=unique(Q1s_tmp);
    cond=numel(Q1s)==2;
    cond1=0;
    if cond
        if ~isempty(strfind(pepLibrary(idxPepLibrary).sequence,'K'))
            cond1=abs(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(idxPepLibrary).charge-settings.deltaLabelK)<settings.massShiftTolerance;
            if ~cond1
                disp(['2 different Q1s for peptide ' currPep ' but label mass shift not respected: actual one in charge 1 would be ' num2str(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(idxPepLibrary).charge) ' instead of ' num2str(settings.deltaLabelK)])
                pause(1)
            end
        elseif ~isempty(strfind(pepLibrary(idxPepLibrary).sequence,'R'))
            cond1=abs(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(idxPepLibrary).charge-settings.deltaLabelR)<settings.massShiftTolerance;
            if ~cond1
                disp(['2 different Q1s for peptide ' currPep ' but label mass shift not respected: actual one in charge 1 would be ' num2str(abs(str2num(Q1s{1})-str2num(Q1s{2}))*pepLibrary(idxPepLibrary).charge) ' instead of ' num2str(settings.deltaLabelR)])
                pause(1)
            end
        end
    end
    pepLibrary(idxPepLibrary).hasRef=cond1;
    
    if pepLibrary(idxPepLibrary).hasRef
        if ~settings.invertedLabel %heavy is the ref
            precursor(i)=(Q1s(1));
            precursorRef(i)=(Q1s(2));
        else %light is the ref
            precursor(i)=(Q1s(2));
            precursorRef(i)=(Q1s(1));
        end
        Q3sRef_tmp=(transList(strmatch((precursorRef(i)),allQ1s),strmatch('Q3',transList(1,:))));
        [Q3sRef_tmp idx_sort_Q3sRef]=unique(delEmptySpaces(Q3sRef_tmp));
        relIntRef_tmp_all=(transList(strmatch((precursorRef(i)),allQ1s),strmatch('relative_intensity',transList(1,:))));
        relIntRef_tmp=relIntRef_tmp_all(idx_sort_Q3sRef);
        for j=1:length(Q3sRef_tmp)
            pepLibrary(idxPepLibrary).Q3sRef(j)=str2num(Q3sRef_tmp{j});
            pepLibrary(idxPepLibrary).relIntRef(j)=str2num(relIntRef_tmp{j});
        end
        pepLibrary(idxPepLibrary).RT_kit_pep=0;
    else
        precursor(i)=(Q1s(1));
        precursorRef(i)={-1};
        
        if strmatch('RT-Kit', pepLibrary(idxPepLibrary).protein)
            pepLibrary(idxPepLibrary).RT_kit_pep=1;
        else
            pepLibrary(idxPepLibrary).RT_kit_pep=0;
        end
    end
    
    Q3s_tmp=(transList(strmatch(precursor(i),allQ1s),strmatch('Q3',transList(1,:))));
    [Q3s_tmp idx_sort_Q3]=unique(delEmptySpaces(Q3s_tmp));
    relInt_tmp_all=(transList(strmatch(precursor(i),allQ1s),strmatch('relative_intensity',transList(1,:))));
    relInt_tmp=relInt_tmp_all(idx_sort_Q3);
    for j=1:length(Q3s_tmp)
        pepLibrary(idxPepLibrary).Q3s(j)=str2num(Q3s_tmp{j});
        pepLibrary(idxPepLibrary).relInt(j)=str2num(relInt_tmp{j});
    end
    pepLibrary(idxPepLibrary).precursor=str2num(precursor{i});
    if ~isequal(precursorRef{i},-1)
        pepLibrary(idxPepLibrary).precursorRef=str2num(precursorRef{i});
    else
        pepLibrary(idxPepLibrary).precursorRef=-1;
    end
    
end


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
function srm=dataReFormat(filePath,out)
% it maps SRM data from mzXML to Matlab struct variables

info =mzxmlinfo(filePath);
precursor=zeros(size(out.scan,1),1);
rt_q3=zeros(size(out.scan,1),1);
mz_q3=cell(size(out.scan,1),1);
int_q3=cell(size(out.scan,1),1);
for i=1:size(out.scan,1)
    precursor(i)=out.scan(i,1).precursorMz.value;
    mz_q3{i}=out.scan(i,1).peaks.mz(1:2:end);
    int_q3{i}=out.scan(i,1).peaks.mz(2:2:end);
    rt_q3(i)=str2num(out.scan(i,1).retentionTime(3:end-1));
end

diff_precursors=unique(precursor);

for k=1:size(diff_precursors,1)
    srm.data(k).Q1=diff_precursors(k);
    idx=find(precursor==srm.data(k).Q1);
    srm.data(k).Q3=unique(cell2mat(mz_q3(idx)));
    srm.data(k).Q3(srm.data(k).Q3==0)=[];
    numOfTrans=numel(srm.data(k).Q3);
    N=numel(idx);
    
    for tN=1:numOfTrans
        trans{tN}=[];
        rts{tN}=[];
    end
    for i=1:N
        cc=idx(i);
        [currentScanMzs idx_unique]=unique(mz_q3{cc});
        currentScanInts =int_q3{cc}(idx_unique);
        for j=1:length(currentScanInts)
            current_mz=currentScanMzs(j);
            transNumber=find(srm.data(k).Q3==current_mz);
            if ~isempty(transNumber)
                trans{transNumber}=cat(2,trans{transNumber},currentScanInts(j));
                rts{transNumber}=cat(2,rts{transNumber},rt_q3(cc));
            else
                disp(['debug! Can''t find data for current Q3 ' num2str(current_mz) ' Da'])
            end
        end
    end
    if ~issorted(rts)
        disp('rts is not sorted! it might affect downstream analysis')
    end
    srm.data(k).trans=trans;
    srm.data(k).rts=rts;
    clear trans rts transOLD rtOLD
    java.lang.System.gc()
end
srm.metadata=info;




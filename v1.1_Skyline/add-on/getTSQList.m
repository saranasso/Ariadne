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
function getTSQList(list0,fileName,header,idx_retTimesCol,idx_dwellTimeCol,schedWinWidth,thresholdThermo,thermoBool,polarity)

% set scheduling window
rt_lowerBound(:,1)=cell2mat(list0(:,idx_retTimesCol))-schedWinWidth/2;
rt_upperBound(:,1)=cell2mat(list0(:,idx_retTimesCol))+schedWinWidth/2;
%compute trans overlapping (counts/hist)
rts_monitored=[];
for j=1:length(rt_upperBound)
    rts_monitored=[rts_monitored, [rt_lowerBound(j):1:rt_upperBound(j)]];
end
hh=figure
hold on
hist(rts_monitored,max(rts_monitored)-min(rts_monitored))
title(['Number of overlapping trans with a scheduling window of '  num2str(schedWinWidth) ' min' ])
% attachPlotToFile(hh,[fileName(1:end-4) '_hist_overlapping_transs'])
print_file=[fileName(1:end-4) '_hist_overlapping_transs.ps'];
print(hh,'-dpsc','-cmyk','-append',print_file)
close(hh)
% set threshold signal intensity for peak picking in Thermo online ret times realignment
tmp=cell(size(list0,1),1);
tmp(:,1)={thresholdThermo};
thresholdThermoCol=tmp;
% set boolean for Thermo online ret times realignment
tmp(:,1)={thermoBool};
thermoBoolCol=tmp;
% polarity column
tmp(:,1)={polarity};
polarityCol=tmp;
% extract data of interest from list
idx_Q1=strmatch('Q1',header,'exact')
idx_Q3=strmatch('Q3',header,'exact')
idx_CE=strmatch('CE',header,'exact')
idx_transIdentifier=strmatch('transition_name',header,'exact') ;
list_TSQ=[list0(:,idx_Q1),list0(:,idx_Q3),list0(:,idx_CE),...
    num2cell(rt_lowerBound),num2cell(rt_upperBound),polarityCol,thresholdThermoCol,thermoBoolCol,list0(:,idx_transIdentifier)]

% write list to csv file
tmpDS=dataset(list_TSQ);
export(tmpDS,'File',[fileName(1:end-4) '-' inputname(1) '.csv'],'Delimiter',',','WriteVarNames',false)

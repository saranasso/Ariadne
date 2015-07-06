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
function y_adj=getBackadjParam(X,noisy,peakWidth)
% it performs background adjustment with a cross validation approach for
% good window size choice.

%% Use "cvpartition" and "crossval"
wins = [peakWidth/2,peakWidth, 2*peakWidth];

steps=wins;
sse = nan(numel(steps),numel(wins));
cp = cvpartition(numel(X),'k',round(numel(X)/ceil(0.1*numel(X))));

for j=1:length(wins)
    for i=1:j % step should rather be smaller or equal to win
        f = @(train,test) norm(test(:,2) - mymsbackadj(train,test(:,1),wins(j),steps(i)))^2;
        try
            sse(i,j) = sum(crossval(f,[X,noisy],'partition',cp));
        end
    end
end
if size(sse,1)==1
    [minsse j_min] = min(sse);
else
    [minsse i_min] = min(sse);
    [minsse_rc j_min]=min(minsse);
    step = steps(i_min(j_min));
end
span = wins(j_min);

if size(sse,1)~=1
    y_adj=mymsbackadj([X,noisy],X,span,step);
elseif size(sse,1)==1
    y_adj=mymsbackadj([X,noisy],X,span);
end






	

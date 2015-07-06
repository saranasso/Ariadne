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
function y_smooth=getSmoothParamLowess(X,y,newX)
% it smooths the trace using a cross validation approach for estimating the
% span (window) size

minSpan=3;
minSpanFract = minSpan/numel(X);
if minSpanFract < 0.25
    spans=(0.25:0.05:0.5);
elseif minSpanFract > 0.5
    spans=(minSpanFract:0.05:0.75);
elseif minSpanFract > 0.75
    spans=(minSpanFract:0.05:0.95);
elseif minSpanFract > 0.25 && minSpanFract < 0.5
    spans=(minSpanFract:0.05:0.5);
elseif minSpanFract > 1
   disp('Less than 3 data points for current peak group, please check it out!')  % it shouldn't happen
end

sse = nan(1,numel(spans));

cp = cvpartition(numel(X),'k',round(numel(X)/ceil(0.1*numel(X))));

for j=1:length(spans)
    
    f = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
    sse(j) = sum(crossval(f,[X,y],'partition',cp));
end

if size(sse,1)==1
    [minsse j_min] = min(sse);
else
    [minsse i_min] = min(sse);
    [minsse_rc j_min]=min(minsse);
    bestDegree = degrees(i_min(j_min));
end
span = spans(j_min);

if nargin==2 && size(sse,1)==1
    y_smooth=mylowess([X,y],X,span);
elseif nargin==3 && size(sse,1)==1
    y_smooth=mylowess([X,y],newX,span);
    X=newX;
end





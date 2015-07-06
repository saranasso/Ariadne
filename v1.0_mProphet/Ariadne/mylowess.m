%% Copyright

%Copyright (c) 2010, The MathWorks, Inc.
%All rights reserved.

%Redistribution and use in source and binary forms, with or without 
%modification, are permitted provided that the following conditions are 
%met:

 %   * Redistributions of source code must retain the above copyright 
 %     notice, this list of conditions and the following disclaimer.
 %   * Redistributions in binary form must reproduce the above copyright 
 %     notice, this list of conditions and the following disclaimer in 
 %     the documentation and/or other materials provided with the distribution
 %   * Neither the name of the The MathWorks, Inc. nor the names 
 %     of its contributors may be used to endorse or promote products derived 
 %     from this software without specific prior written permission.
      
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.

function ys=mylowess(xy,xs,span)
%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.

if nargin<3 || isempty(span)
    span = .3;
end

% Sort and get smoothed version of xy data
xy = sortrows(xy);
x1 = xy(:,1);
y1 = xy(:,2);
ys1 = smooth(x1,y1,span,'loess');

% Remove repeats so we can interpolate
t = diff(x1)==0;
x1(t)=[]; ys1(t) = [];

% Interpolate to evaluate this at the xs values
% ys = interp1(x1,ys1,xs,'linear',NaN);
ys = interp1(x1,ys1,xs,'spline',NaN);

% Some of the original points may have x values outside the range of the
% resampled data.  Those are now NaN because we could not interpolate them.
% Replace NaN by the closest smoothed value.  
if any(isnan(ys))
    ys(xs<x1(1)) = ys1(1);
    ys(xs>x1(end)) = ys1(end);
end
	


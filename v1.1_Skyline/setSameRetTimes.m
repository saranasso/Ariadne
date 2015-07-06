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
function [ref,endo,rts_endo,rts_ref] = setSameRetTimes(ref,endo,rts_endo,rts_ref)
% setting the same retention times to enable downstream data analysis

for i=1:size(ref,2)
    rti_min=min(rts_endo{i}(1),rts_ref{i}(1));
    rtf_max=max(rts_endo{i}(end),rts_ref{i}(end));
    newRetTimes=linspace(rti_min,rtf_max,(rtf_max-rti_min))';
    % Remove repeats so we can interpolate
    t = diff(rts_endo{i})==0;
    rts_endo{i}(t)=[]; endo{i}(t) = [];
    endo{i} = interp1(rts_endo{i}, endo{i},newRetTimes,'spline',0);
    t = diff(rts_ref{i})==0;
    rts_ref{i}(t)=[]; ref{i}(t) = [];
    ref{i} = interp1(rts_ref{i}, ref{i},newRetTimes,'spline',0);
    rts_endo{i}=newRetTimes;
    rts_ref{i}=newRetTimes;
end







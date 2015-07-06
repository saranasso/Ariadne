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
function [idx_peak meas]=getPeakIdx(peaklist,prediction)
% refining the position of the best peak group: selecting the closest peak to the prediction

if iscell(peaklist)
    if ~isempty(peaklist)
        for i=1:size(peaklist,2)
            for j=1:size([peaklist{i}],1)
                tmp(j,:)=abs(peaklist{i}(j,:)-prediction);
            end
            [meas(i,:) idx_peak(i,:)]=min(tmp,[],1);
            clear tmp
        end
        idx_peak=(idx_peak(:,1));
    else
        idx_peak(i,:)=0;
        meas(i,:)=0;
    end
else
    if ~isempty(peaklist)
        tocluster=[prediction;peaklist];
        distances=pdist(tocluster);
        interesting_distances=distances(1:(size(tocluster,1)-1));
        [minValue idx_peak]=min(interesting_distances);
        meas=0;
    else
        idx_peak=0;
        meas=0;
    end
end

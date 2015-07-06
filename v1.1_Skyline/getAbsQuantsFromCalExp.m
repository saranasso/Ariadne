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
function [estimatedAbundances pepLibrarySRM] = getAbsQuantsFromCalExp(pepLibrarySRM, association_file, pepLibrarySRMCal)
% it estimates absolute abundances

global settings
diary('Ariadne.log')

mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};
if settings.multFactor~= 1
    for i=1:numel(pepLibrarySRM)
        pepLibrarySRM(i).ratio = pepLibrarySRM(i).ratio.*settings.multFactor;
        for j =1:numel(pepLibrarySRM(i).ratio_mProphet)
            for k=1:numel(mProphet_methods)
                if isfield(pepLibrarySRM(i).ratio_mProphet(j), mProphet_methods{k})
                    eval([ 'pepLibrarySRM(i).ratio_mProphet(j).' mProphet_methods{k} '= pepLibrarySRM(i).ratio_mProphet(j).' mProphet_methods{k} '*settings.multFactor;'])
                end
            end
        end
    end
end

if ~settings.interferingBackground
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExp(pepLibrarySRM,association_file,pepLibrarySRMCal)
else
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExpWithSameBackground(pepLibrarySRM,association_file,pepLibrarySRMCal)
end

diary off

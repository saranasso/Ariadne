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

% it tries to rescue the peptide traces merging them altogether (separately for endogenous and modified peptide)

[sorted_corr idx_sorted_xCorr]=sort(corr_coeff(1,2,:));
data.endo=[];
data.ref=[];
for i=1:size(trans.ref,2)
    try
        data.endo=cat(1,data.endo,hist2data(endo{idx_sorted_xCorr(i)},rts_endo{idx_sorted_xCorr(i)}));
        data.ref=cat(1,data.ref,hist2data(ref{idx_sorted_xCorr(i)},rts_ref{idx_sorted_xCorr(i)}));
    catch % if for at least one of the trans no peak was detected simply return
        disp(['No peak detected for trans #' num2str(i)])
        return
    end
end
if ~isempty(data.ref) && ~isempty(data.endo)
    endo{1}=hist(data.endo,rts_endo{1})';
    ref{1}=hist(data.ref,rts_ref{1})';
else
    return
end

[endo(1) ref(1)] = transProcessing (endo(1), ref(1), rts_endo(1), rts_ref(1));

%     [xcorrs lags]=xcorr(ref{1},endo{1},'unbiased');
[xcorrs lags]=xcorr(ref{1},endo{1});
[maxXcorr idx_max]=max(xcorrs);
lags(idx_max)
[corr_coeff(:,:,1),pValue(:,:,1),corr_coeff_low(:,:,1),corr_coeff_up(:,:,1)]=corrcoef(ref{1},endo{1});
corr_coeff(1,2,1)

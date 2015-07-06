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
function [gmmFitted, X, dist] = fitPeakGroupsGMM(rti,rtf,gm,prior,endo,rts,pepLibraryInfo)
% it fits the GMM model and recognizes the best fitting one as the one
% closest in weighted euclidean distance to the a priori model parameters
% of relative intensities.

idxs_endo=find(rts>rti & rts<rtf);
for j=1:size(endo,2)
    endo_tmp{j}=endo{j}(idxs_endo)';
    rts_endo{j}=rts(idxs_endo)';
end

try
    tic
    X=[];
    for i=1:size(endo,2)
        rts_data=hist2data(endo_tmp{i},rts_endo{i});
        tmp_data=[rts_data  pepLibraryInfo.Q3s(i)*ones(size(rts_data,1),size(rts_data,2))];
        X=cat(1,tmp_data,X);
    end
    
    disp('gmm')
    
    gmmFitted=gmdistribution.fit(X,size(endo,2),'Start',prior); % non e' confrontabile
    [post,nLogL] = posterior(gmmFitted,X);
    toc
    %Use a function handle to compute a distance that weights each coordinate contribution differently.
    Wgts = [];
    for i=1:numel(gm.PComponents)
        numberOfPtsxTrs=numel(find(post(:,i))==1);
        Wgts = [Wgts, numberOfPtsxTrs.^(-1) ];
    end
    weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
    dist = pdist2(gm.PComponents,gmmFitted.PComponents, @(Xi,Xj) weuc(Xi,Xj,Wgts));
catch
    disp('no gmm')
    dist =NaN;
    gmmFitted = [];
    X=[];
end

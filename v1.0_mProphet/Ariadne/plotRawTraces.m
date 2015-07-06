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

h=figure;
cmap = lines(numel(trans.ref));
for i=1:size(trans.ref,2)
    subplot 211
    plot(rtsFilt.endo{i}, trans.endo{i},'-','Color',cmap(i,:))
    hold on
    plot(rts_endo{i}, endo{i},'g--')
    try
        if ~settings.mProphetFilter
            stem(retTime,max(maxApex.endo),'b*')
        end
        stem(expRetTime,max(maxApex.endo),'r*')
        stem(mRetTime,mMaxApex.endo,'k*')
    end
    if settings.mProphetFilter
        legend([{'full trans'}; {'selected sub-trans'}; {'Expected ret time'}; {'mProphet ret time'}])
    else
        legend([{'full trans'}; {'selected sub-trans'}; {'Ariadne ret time'}; {'Expected ret time'}; {'mProphet ret time'}])
    end
    subplot 212
    plot(rtsFilt.ref{i},trans.ref{i},'-','Color',cmap(i,:))
    hold on
    plot(rts_ref{i},ref{i},'g--')
    try
        if ~settings.mProphetFilter
            stem(retTime,max(maxApex.ref),'b*')
        end
        stem(expRetTime,max(maxApex.ref),'r*')
        stem(mRetTime,mMaxApex.ref,'k*')
    end
    if settings.mProphetFilter
        legend([{'full trans'}; {'selected sub-trans'}; {'Expected ret time'}; {'mProphet ret time'}])
    else
        legend([{'full trans'}; {'selected sub-trans'}; {'Ariadne ret time'}; {'Expected ret time'}; {'mProphet ret time'}])
    end
    
end
title({['k is ' num2str(idx) ' and i is ' num2str(sub_idx)];['Peptide Sequence: ' pepLibrary(idx).sequence];['File name: ' strrep(pepLibrary(idx).dataFile{sub_idx},'_',' ')]})
attachPlotToFile(h,'Ariadne_plots')
close all

java.lang.System.gc();
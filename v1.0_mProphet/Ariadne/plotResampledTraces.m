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

try
    h=figure;
    cmap = lines(numel(trans.ref));
    for i=1:size(trans.ref,2)
        subplot 211
        plot(rts_endo{i}, endo{i},'-','Color',cmap(i,:))
        hold on
        title([{'Resampled sub-trans for endogenous peptide'}])
        subplot 212
        plot(rts_ref{i},ref{i},'-','Color',cmap(i,:))
        hold on
        title([{'Resampled sub-trans for reference peptide'}])
    end
    title({['k is ' num2str(idx) ' and i is ' num2str(sub_idx)];['Peptide Sequence: ' pepLibrary(idx).sequence];['File name: ' strrep(pepLibrary(idx).dataFile{sub_idx},'_',' ')]})
    toc
    attachPlotToFile(h,'Ariadne_plots')
    
    close all
end
java.lang.System.gc();
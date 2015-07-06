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
function boxplots(commonQuants,molarities,methods,mProphet_methods)    

whichstatsOption = '{''numel''}';
for numMethod=1:numel(methods)
        eval(['[x g] = straightenX(commonQuants.' methods{numMethod} ',molarities);']);
        eval(['numberOfCommonEstimatedAbundances.' methods{numMethod} ' = grpstats(x,g,' whichstatsOption ');']);
        
    end
    
options={'extrememode','compress','colors'};
cmap = lines(numel(molarities));
% options={[],[]}
YdataLim=[log2(molarities(end))-0.5*abs(log2(molarities(end))),log2(molarities(1))+0.5*abs(log2(molarities(1)))];
if YdataLim(1)==0
    YdataLim(1)= -3;
    YdataLim(2)= 9.3;
end

YdataLim(2)= 10;
YdataLim(1)= -10;
% YdataLim(2)= 8;
% YdataLim(1)= -4;

for i=1:numel(mProphet_methods)
    h=figure;
    subplot 211
    boxplot(log2(commonQuants.Ariadne),molarities,options{3},cmap)
    xlabel('True values (fmoles)')
    ylabel('Log2 values')
    title('Estimated  Absolute Quants: Ariadne')
    ylim('manual')
    ylim([YdataLim(1) YdataLim(2)])
    hold on
    xl=xlim;
    for j=1:numel(molarities)
        plot([xl(1) xl(2)],[log2(molarities(j)) log2(molarities(j))],'--','Color',cmap((end-j+1),:))
    end
    hbp = findobj(gca,'Tag','Box');
    for j=1:length(hbp)
%         text(median(get(hbp(j),'XData')),(log2(molarities(1))+0.3*abs(log2(molarities(1)))),num2str(numberOfCommonEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
        text(median(get(hbp(j),'XData')),(log2(molarities(end))- 3),num2str(numberOfCommonEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
    end
    
    subplot 212
    eval(['boxplot(log2(commonQuants.' mProphet_methods{i} '),molarities,options{3},cmap)'])
    xlabel('True values (fmoles)')
    ylabel('Log2 values')
    title(['Estimated  Absolute Quants: mProphet method ' mProphet_methods{i}])
    ylim('manual')
    ylim([YdataLim(1) YdataLim(2)])
    hold on
    xl=xlim;
    for j=1:numel(molarities)
        plot([xl(1) xl(2)],[log2(molarities(j)) log2(molarities(j))],'--','Color',cmap(end-j+1,:))
    end
    hbp = findobj(gca,'Tag','Box');
    for j=1:length(hbp)
%         text(median(get(hbp(j),'XData')),(log2(molarities(1))+0.3*abs(log2(molarities(1)))),num2str(eval(['numberOfCommonEstimatedAbundances.' mProphet_methods{i} '(length(hbp)-j+1)'])),'FontSize',10)
        text(median(get(hbp(j),'XData')),(log2(molarities(end))- 3),num2str(eval(['numberOfCommonEstimatedAbundances.' mProphet_methods{i} '(length(hbp)-j+1)'])),'FontSize',10)
    end
    
%     pause
%     attachPlotToFile(h,'boxplots')
%     %       pause
%     close all
end
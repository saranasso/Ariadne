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

global settings

% [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExp(pepLibrarySRM,association_file,pepLibrarySRMCal)
[estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExpWithSameBackground(pepLibrarySRM,association_file,pepLibrarySRMCal)


groups=sort(unique([association_file{:,1}]),'descend');
groups(groups>1)=round(groups(groups>1));

diary on
diary('Ariadne.log')

[RE_A sd_A CV_A]=getAccPrec(estimatedAbundances.Ariadne,groups');


if settings.mProphet
    mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};
    for i=1:numel(mProphet_methods)
        disp(mProphet_methods{i})
        eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(estimatedAbundances.' mProphet_methods{i} ',transpose(groups));']);
        m_RE(i)=nanmedian((RE(i)));
        m_CV(i)=nanmedian((CV(i)));
    end
    RE(RE==0)=nan;
    CV(CV==0)=nan;

    
    estimatedAbundancesTmp=estimatedAbundances;
    whichstatsOption = '{''numel''}';
    methods=cat(1,{'Ariadne'},mProphet_methods);
    for numMethod=1:numel(methods)
        estimatedAbundancesTmp.Ariadne(isinf(estimatedAbundancesTmp.Ariadne))=nan;
        eval(['[x g] = straightenX(estimatedAbundancesTmp.' methods{numMethod} ',groups);']);
        eval(['numberOfEstimatedAbundances.' methods{numMethod} ' = grpstats(x,g,' whichstatsOption ');']);
        
    end
    %% plot acc-prec
    h = figure;
    % set(gca, 'XTick',sort(log10(groups)),'XTickLabel',sort(groups))
    hold on
    plot(log10(groups),log2(abs([RE_A';RE])),'o-')
    plot(log10(groups),log2([CV_A';CV]),'*:')
    set(gca, 'XTick',sort(log10(groups)),'XTickLabel',sort(groups))
    legend(methods, 'Location', 'southeast' )
    xlabel('Dilution groups (fmol)')
    ylabel('Log2(%)')
    title('Relative Error (solid line) & Coefficient of Variation (dotted line) for complex background data')
    attachPlotToFile(h,'accPrec')
    close(h)
    
    
    
    %% efficiency
    for i=1:numel(methods)
        methods{i}
        eval(['efficiency.' methods{i} ' = nanmedian(numberOfEstimatedAbundances.' methods{i} ')/numel(pepLibrarySRM)*100;']);
        eval(['ariadnesEfficiencyGain.' methods{i} '= (nanmedian(numberOfEstimatedAbundances.Ariadne)-nanmedian(numberOfEstimatedAbundances.' methods{i} '))/nanmedian(numberOfEstimatedAbundances.' methods{i} ')*100;'])
    end
    %% efficiency
    for i=1:numel(methods)
        methods{i}
        eval(['efficiency.' methods{i} ' = nanmedian(numberOfEstimatedAbundances.' methods{i} ')/numel(pepLibrarySRM)*100;']);
        eval(['ariadnesEfficiencyGain.' methods{i} '= (nanmedian(numberOfEstimatedAbundances.Ariadne)-nanmedian(numberOfEstimatedAbundances.' methods{i} '))/nanmedian(numberOfEstimatedAbundances.' methods{i} ')*100;'])
    end
    ariadnesEfficiencyGain
end

diary off
% return
%
%
% targetsNumber=numel([pepLibrarySRM(:).dataFile]);
%
% absquant.Ariadne=[pepLibrarySRM(:).estimatedAbundance];
% absquant.Ariadne=absquant.Ariadne(absquant.Ariadne>0);
% absquant.Ariadne=absquant.Ariadne(~isnan(absquant.Ariadne));
% absquant_efficiency.Ariadne= numel(absquant.Ariadne)/targetsNumber*100;
% % total number of targets
% % targetsNumber=numel([pepLibrarySRM(:).dataFile]);
% ratios.Ariadne=[pepLibrarySRM(:).ratio];
% ratios.Ariadne=ratios.Ariadne(ratios.Ariadne>0);
% ratios.Ariadne=ratios.Ariadne(~isnan(ratios.Ariadne));
% efficiency.Ariadne=numel(ratios.Ariadne)/targetsNumber*100;
% for i=1:numel(mProphet_methods)
%     eval(['ratios.' mProphet_methods{i} '=[];'])
% end
% for i=1:numel(mProphet_methods)
%     for j=1:numel(pepLibrarySRM)
%         try
%             eval(['ratios.' mProphet_methods{i} '=cat(2, ratios.' mProphet_methods{i} ',[pepLibrarySRM(j).ratio_mProphet.' mProphet_methods{i} ']);'])
%         end
%     end
% end
% for i=1:numel(mProphet_methods)
%     eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '(~isnan(ratios.' mProphet_methods{i} '));'])
%     eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '((ratios.' mProphet_methods{i} ')>0);'])
%     eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '((ratios.' mProphet_methods{i} ')~=1);'])
%     eval(['efficiency.' mProphet_methods{i} '=numel(ratios.' mProphet_methods{i} ')/targetsNumber*100;'])
%     eval(['ariadnesEfficiencyGain.' mProphet_methods{i} '=(numel(ratios.Ariadne)-numel(ratios.' mProphet_methods{i} '))/numel(ratios.' mProphet_methods{i} ')*100;'])
% end

%% boxplots

options={'extrememode','compress','colors'};
cmap = lines(numel(groups));
% options={[],[]}
YdataLim=[log2(groups(end))-2,log2(groups(1))+1];
% YdataLim=[log2(groups(end))-1,log2(groups(1))+1];

for i=1:numel(mProphet_methods)
    h=figure;
    subplot 211
    %     boxplot(log2(estimatedAbundances.Ariadne),groups,options{1},options{2},options{3},cmap)
    boxplot(log2(estimatedAbundances.Ariadne),groups,options{3},cmap)
    xlabel('True values (fmoles)')
    ylabel('Log2 values')
    title('Estimated  Absolute Quants: Ariadne')
    ylim('manual')
    ylim([YdataLim(1) YdataLim(2)])
    hold on
    xl=xlim;
    for j=1:numel(groups)
        plot([xl(1) xl(2)],[log2(groups(j)) log2(groups(j))],'--','Color',cmap((end-j+1),:))
    end
    hbp = findobj(gca,'Tag','Box');
    for j=1:length(hbp)
        text(median(get(hbp(j),'XData')),(log2(max(groups))+0.05*abs(log2(max(groups)))),num2str(numberOfEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
    end
    
    subplot 212
    %     eval(['boxplot(log2(estimatedAbundances.' mProphet_methods{i} '),groups,options{1},options{2},options{3},cmap)'])
    eval(['boxplot(log2(estimatedAbundances.' mProphet_methods{i} '),groups,options{3},cmap)'])
    xlabel('True values (fmoles)')
    ylabel('Log2 values')
    title(['Estimated  Absolute Quants: mProphet method ' mProphet_methods{i}])
    ylim('manual')
    ylim([YdataLim(1) YdataLim(2)])
    hold on
    xl=xlim;
    for j=1:numel(groups)
        plot([xl(1) xl(2)],[log2(groups(j)) log2(groups(j))],'--','Color',cmap(end-j+1,:))
    end
    hbp = findobj(gca,'Tag','Box');
    for j=1:length(hbp)
        text(median(get(hbp(j),'XData')),(log2(max(groups))+0.05*abs(log2(max(groups)))),num2str(eval(['numberOfEstimatedAbundances.' mProphet_methods{i} '(length(hbp)-j+1)'])),'FontSize',10)
    end
    
    %   pause
    attachPlotToFile(h,'boxplots')
    %   pause
    close all
end

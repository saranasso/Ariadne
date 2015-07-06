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
    association_dil_file = association_file; % necessary to cope with different versions
end

molarities=sort(unique([association_dil_file{:,1}]),'descend');
molarities(molarities>1)=round(molarities(molarities>1));
mProphet_methods={'totalxic';'maxapex';'apexsum';'apexsum_outlier'};

methods=cat(1,{'Ariadne'},mProphet_methods);

for i=1:numel(methods)
    
    eval(['estimatedAbundances.' methods{i} '(estimatedAbundances.' methods{i} '<0) = nan;'])
    
end

disp('Ariadne')
[RE_A sd_A CV_A]=getAccPrec(estimatedAbundances.Ariadne,molarities');
m_RE_A=nanmedian(RE_A)
m_CV_A=nanmedian(CV_A)

if settings.mProphet
    disp('mProphet')
    for i=1:numel(mProphet_methods)
        disp(mProphet_methods{i})
        eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(estimatedAbundances.' mProphet_methods{i} ',transpose(molarities));']);
    end
    for i=1:numel(mProphet_methods)
        accuracyGain(i,:)=abs((RE(i,:)))-abs((RE_A'));
        precisionGain(i,:)=abs((CV(i,:)))-abs((CV_A'));
    end
    
    estimatedAbundancesTmp=estimatedAbundances;
    whichstatsOption = '{''numel''}';
    estimatedAbundancesTmp.Ariadne(isinf(estimatedAbundancesTmp.Ariadne))=nan;
    for numMethod=1:numel(methods)
        eval(['[x g] = straightenX(estimatedAbundancesTmp.' methods{numMethod} ',molarities);']);
        eval(['numberOfEstimatedAbundances.' methods{numMethod} ' = grpstats(x,g,' whichstatsOption ');']);
        
    end
    
    %% plot acc-prec
    h = figure;
    plot(log10(molarities),log2(abs([RE_A';RE])),'o-')
    hold on
    plot(log10(molarities),log2([CV_A';CV]),'*:')
    legend(methods)
    set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
    legend(methods, 'Location', 'southeast' )
    xlabel('Dilution groups (fmol)')
    ylabel('Log2(%)')
    title('Relative Error (solid line) & Coefficient of Variation (dotted line)')
    pause
    attachPlotToFile(h,'accPrec')
    close(h)
    
    % excluding totalxic
    h = figure;
    plot(log10(molarities),log2(abs([RE_A';RE([2,3,4],:)])),'o-')
    hold on
    plot(log10(molarities),log2([CV_A';CV([2,3,4],:)]),'*:')
    set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
    legend(methods([1,3,4,5]))
    xlabel('Dilution groups Log10(fmol)')
    ylabel('Log2(%)')
    title('Relative Error (solid line) & Coefficient of Variation (dotted line)')
    pause
    attachPlotToFile(h,'accPrec')
    close(h)
    
    
    %% efficiency
    for i=1:numel(methods)
        methods{i}
        eval(['efficiency.' methods{i} ' = (numberOfEstimatedAbundances.' methods{i} ')./numel(pepLibrarySRM)*100;']);
        eval(['ariadnesEfficiencyGain.' methods{i} '= ((numberOfEstimatedAbundances.Ariadne)-nanmedian(numberOfEstimatedAbundances.' methods{i} '))/nanmedian(numberOfEstimatedAbundances.' methods{i} ')*100;'])
    end
    
    %% plot quant efficiency
    for i=1:numel(methods)
        plotEfficiency(i,:) = eval(['numberOfEstimatedAbundances.' methods{i} ';']);
        % save('plotEff','plotEfficiency','molarities')
    end
  
    %% diary off
    diary off
    
    
    %% boxplots
    options={'extrememode','compress','colors'};
    cmap = lines(numel(molarities));
    YdataLim=[log2(molarities(end))-0.5*abs(log2(molarities(end))),log2(molarities(1))+0.5*abs(log2(molarities(1)))];
    if YdataLim(1)==0
        YdataLim(1)= -3;
        YdataLim(2)= 9.3;
    end
    
    YdataLim(2)= 11;
    YdataLim(1)= -7;
    
    for i=1:numel(mProphet_methods)
        h=figure;
        subplot 211
        boxplot(log2(estimatedAbundances.Ariadne),molarities,options{3},cmap)
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
            text(median(get(hbp(j),'XData')),(log2(molarities(1))+0.3*abs(log2(molarities(1)))),num2str(numberOfEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
        end
        
        subplot 212
        eval(['boxplot(log2(estimatedAbundances.' mProphet_methods{i} '),molarities,options{3},cmap)'])
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
            text(median(get(hbp(j),'XData')),(log2(molarities(1))+0.3*abs(log2(molarities(1)))),num2str(eval(['numberOfEstimatedAbundances.' mProphet_methods{i} '(length(hbp)-j+1)'])),'FontSize',10)
        end
        
        pause
        attachPlotToFile(h,'boxplots')
        %       pause
        close all
    end
else
    %% plot acc-prec
    h = figure;
    plot(log10(molarities),log2(abs([RE_A'])),'o-')
    hold on
    plot(log10(molarities),log2([CV_A']),'*:')
    legend('Ariadne')
    set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
    legend('Ariadne'           , 'Location', 'southeast' )
    xlabel('Dilution groups (fmol)')
    ylabel('Log2(%)')
    title('Relative Error (solid line) & Coefficient of Variation (dotted line)')
    pause
    attachPlotToFile(h,'accPrec')
    close(h)
    
    %% boxplots
    options={'extrememode','compress','colors'};
    cmap = lines(numel(molarities));
    YdataLim=[log2(molarities(end))-0.5*abs(log2(molarities(end))),log2(molarities(1))+0.5*abs(log2(molarities(1)))];
    if YdataLim(1)==0
        YdataLim(1)= -3;
        YdataLim(2)= 9.3;
    end 
    YdataLim(2)= 11;
    YdataLim(1)= -7;
    h=figure;
    boxplot(log2(estimatedAbundances.Ariadne),molarities,options{3},cmap)
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
        text(median(get(hbp(j),'XData')),(log2(molarities(1))+0.3*abs(log2(molarities(1)))),num2str(numberOfEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
    end
end

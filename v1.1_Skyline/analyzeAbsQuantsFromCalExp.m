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
function [estimatedAbundances pepLibrarySRM] = analyzeAbsQuantsFromCalExp(pepLibrarySRM, association_file, estimatedAbundances)


global settings
diary('Ariadne.log')

groups=sort(unique([association_file{:,1}]),'descend');
groups(groups>1)=round(groups(groups>1));

%% efficiency
targetsNumber=numel([pepLibrarySRM(:).dataFile]);

absquant.Ariadne=[pepLibrarySRM(:).estimatedAbundance];
absquant.Ariadne=absquant.Ariadne(absquant.Ariadne>0);
absquant.Ariadne=absquant.Ariadne(~isnan(absquant.Ariadne));
absquant_efficiency.Ariadne= numel(absquant.Ariadne)/targetsNumber*100;
ratios.Ariadne=[pepLibrarySRM(:).ratio];
ratios.Ariadne=ratios.Ariadne(ratios.Ariadne>0);
ratios.Ariadne=ratios.Ariadne(~isnan(ratios.Ariadne));
efficiency.Ariadne=numel(ratios.Ariadne)/targetsNumber*100;
for i=1:numel(mProphet_methods)
    eval(['ratios.' mProphet_methods{i} '=[];'])
end
for i=1:numel(mProphet_methods)
    for j=1:numel(pepLibrarySRM)
        try
            eval(['ratios.' mProphet_methods{i} '=cat(2, ratios.' mProphet_methods{i} ',[pepLibrarySRM(j).ratio_mProphet.' mProphet_methods{i} ']);'])
        end
    end
end
for i=1:numel(mProphet_methods)
    eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '(~isnan(ratios.' mProphet_methods{i} '));'])
    eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '((ratios.' mProphet_methods{i} ')>0);'])
    eval(['ratios.' mProphet_methods{i} '=ratios.' mProphet_methods{i} '((ratios.' mProphet_methods{i} ')~=1);'])
    eval(['efficiency.' mProphet_methods{i} '=numel(ratios.' mProphet_methods{i} ')/targetsNumber*100;'])
    eval(['ariadnesEfficiencyGain.' mProphet_methods{i} '=(numel(ratios.Ariadne)-numel(ratios.' mProphet_methods{i} '))/numel(ratios.' mProphet_methods{i} ')*100;'])
end


if settings.dilution_series
    disp('Ariadne')
    [RE_A sd_A CV_A]=getAccPrec(estimatedAbundances.Ariadne,groups');
    nanmean(RE_A)
    nanmean(CV_A)
    
    if settings.mProphet
        disp('mProphet')
        for i=1:numel(mProphet_methods)
            disp(mProphet_methods{i})
            eval(['[RE sd CV]=getAccPrec(estimatedAbundances.' mProphet_methods{i} ',transpose(groups));']);
            m_RE(i)=nanmean(RE)
            m_CV(i)=nanmean(CV)
        end
        [min_RE idx_best_RE_mProphet]=min(m_RE);
        [min_CV idx_best_CV_mProphet]=min(m_CV);
        eval(['[RE sd CV]=getAccPrec(estimatedAbundances.' mProphet_methods{idx_best_RE_mProphet} ',transpose(groups));']);
        accuracyGain=abs(nanmean(RE))-abs(nanmean(RE_A))
        precisionGain=abs(nanmean(CV))-abs(nanmean(CV_A))
        
        
        estimatedAbundancesTmp=estimatedAbundances;
        whichstatsOption = '{''numel''}';
        methods=cat(1,{'Ariadne'},mProphet_methods);
        for numMethod=1:numel(methods)
            estimatedAbundancesTmp.Ariadne(isinf(estimatedAbundancesTmp.Ariadne))=nan;
            eval(['[x g] = straightenX(estimatedAbundancesTmp.' methods{numMethod} ',groups);']);
            eval(['numberOfEstimatedAbundances.' methods{numMethod} ' = grpstats(x,g,' whichstatsOption ');']);
            
        end
        
        %% boxplots
        
        options={'extrememode','compress','colors'};
        cmap = lines(numel(groups));
        % options={[],[]}
        % YdataLim=[log10(groups(end))-0.5*abs(log10(groups(end))),log10(groups(1))+0.5*abs(log10(groups(1)))];
        YdataLim=[log10(groups(end))-1,log10(groups(1))+1];
        
        for i=1:numel(mProphet_methods)
            h=figure;
            subplot 211
            boxplot(log10(estimatedAbundances.Ariadne),groups,options{1},options{2},options{3},cmap)
            xlabel('True values (fmoles)')
            ylabel('Log10 values')
            title('Estimated  Absolute Quants: Ariadne')
            ylim('manual')
            ylim([YdataLim(1) YdataLim(2)])
            hold on
            xl=xlim;
            for j=1:numel(groups)
                plot([xl(1) xl(2)],[log10(groups(j)) log10(groups(j))],'--','Color',cmap((end-j+1),:))
            end
            hbp = findobj(gca,'Tag','Box');
            for j=1:length(hbp)
                text(median(get(hbp(j),'XData')),(log10(groups(1))+0.3*abs(log10(groups(1)))),num2str(numberOfEstimatedAbundances.Ariadne(length(hbp)-j+1)),'FontSize',10)
            end
            
            subplot 212
            eval(['boxplot(log10(estimatedAbundances.' mProphet_methods{i} '),groups,options{1},options{2},options{3},cmap)'])
            xlabel('True values (fmoles)')
            ylabel('Log10 values')
            title(['Estimated  Absolute Quants: mProphet method ' mProphet_methods{i}])
            ylim('manual')
            ylim([YdataLim(1) YdataLim(2)])
            hold on
            xl=xlim;
            for j=1:numel(groups)
                plot([xl(1) xl(2)],[log10(groups(j)) log10(groups(j))],'--','Color',cmap(end-j+1,:))
            end
            hbp = findobj(gca,'Tag','Box');
            for j=1:length(hbp)
                text(median(get(hbp(j),'XData')),(log10(groups(1))+0.3*abs(log10(groups(1)))),num2str(eval(['numberOfEstimatedAbundances.' mProphet_methods{i} '(length(hbp)-j+1)'])),'FontSize',10)
            end
            
            %   pause
            attachPlotToFile(h,'boxplots')
            %   pause
            close all
        end
    end
end
diary off

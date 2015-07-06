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
pause off
interferingBackground = 0;
if ~interferingBackground
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExp(pepLibrarySRM,association_dil_file,pepLibrarySRMCal)
else
    [estimatedAbundances pepLibrarySRM]=getAbsoluteQuantsFromCalExpWithSameBackground(pepLibrarySRM,association_file,pepLibrarySRMCal)
end

exportAriadneSkylineResults2Csv
exportResults2Txt %only Ariadne

try
    association_dil_file = association_file; % necessary to cope with different versions
end

molarities=sort(unique([association_dil_file{:,1}]),'descend');
molarities(molarities>1)=round(molarities(molarities>1));

mProphet_methods={'Skyline_mProphet'};



methods=cat(1,{'Ariadne'},mProphet_methods);

for i=1:numel(methods)
    
    eval(['estimatedAbundances.' methods{i} '(estimatedAbundances.' methods{i} '<=0) = nan;'])
    
end

disp('Ariadne')
[RE_A sd_A CV_A]=getAccPrec(estimatedAbundances.Ariadne,molarities');

if settings.mProphet
    disp('mProphet')
    for i=1:numel(mProphet_methods)
        disp(mProphet_methods{i})
        eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(estimatedAbundances.' mProphet_methods{i} ',transpose(molarities));']);
    end
    for i=1:numel(mProphet_methods)
        accuracyGain(i,:)=abs((RE(i,:)))./abs((RE_A'))
        precisionGain(i,:)=abs((CV(i,:)))./abs((CV_A'))
    end
    
    estimatedAbundancesTmp=estimatedAbundances;
    whichstatsOption = '{''numel''}';
    estimatedAbundancesTmp.Ariadne(isinf(estimatedAbundancesTmp.Ariadne))=nan;
    for numMethod=1:numel(methods)
        eval(['[x g] = straightenX(estimatedAbundancesTmp.' methods{numMethod} ',molarities);']);
        eval(['numberOfEstimatedAbundances.' methods{numMethod} ' = grpstats(x,g,' whichstatsOption ');']);        
    end
    
    
    %% efficiency
    for i=1:numel(methods)
        methods{i}
        eval(['efficiency.' methods{i} ' = (numberOfEstimatedAbundances.' methods{i} ')/numel(pepLibrarySRM)*100;']);
        %         eval(['ariadnesEfficiencyGain.' methods{i} '= ((numberOfEstimatedAbundances.Ariadne)-nanmedian(numberOfEstimatedAbundances.' methods{i} '))/nanmedian(numberOfEstimatedAbundances.' methods{i} ')*100;'])
        eval(['ariadnesEfficiencyGain.' methods{i} '= ((numberOfEstimatedAbundances.Ariadne))/(numberOfEstimatedAbundances.' methods{i} ')*100;'])
    end
end
%% plot quant efficiency
for i=1:numel(methods)
    plotEfficiency(i,:) = eval(['numberOfEstimatedAbundances.' methods{i} ';']);
    % save('plotEff','plotEfficiency','molarities')
end

%% diary off
diary off


%% boxplots
boxplots(estimatedAbundances,molarities,methods,mProphet_methods)


%% correlation/consensus ALL quants Ariadne and Skyline
corrValueAll = getConsensus(estimatedAbundances)


%% corrections by Skyline?
[i,j]=find(isnan(estimatedAbundances.Ariadne));
corrections.Skyline_mProphet = accumarray([i,j], estimatedAbundances.Skyline_mProphet(isnan(estimatedAbundances.Ariadne)))
corrections.Skyline_mProphet(corrections.Skyline_mProphet==0)=nan;
disp('mProphet')
for i=1:numel(mProphet_methods)
    disp(mProphet_methods{i})
    eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(corrections.' mProphet_methods{i} ',transpose(molarities));'])
end

%% corrections by Ariadne?
[i,j]=find(isnan(estimatedAbundances.Skyline_mProphet));
corrections.Ariadne = accumarray([i,j], estimatedAbundances.Ariadne(isnan(estimatedAbundances.Skyline_mProphet)))
corrections.Ariadne(corrections.Ariadne==0)=nan;
disp('Ariadne')
[RE_A sd_A CV_A]=getAccPrec(corrections.Ariadne,molarities')

%% corrections quality
for i=1:numel(mProphet_methods)
    accuracyGain(i,:)=abs((RE(i,:)))./abs((RE_A'))
    precisionGain(i,:)=abs((CV(i,:)))./abs((CV_A'))
end
[CV_A,CV']
[RE_A,RE']
boxplots(corrections,molarities,methods,mProphet_methods)

%% commonly quantified peps - gains
commonQuants = estimatedAbundancesTmp;
commonQuants.Ariadne(isnan(commonQuants.Skyline_mProphet)) = nan;
commonQuants.Skyline_mProphet(isnan(commonQuants.Ariadne)) = nan;
%% correlation/consensus COMMON quants Ariadne and Skyline
corrValueCommon = getConsensus(commonQuants)
disp('Ariadne')
[RE_A sd_A CV_A]=getAccPrec(commonQuants.Ariadne,molarities');
m_RE_A=nanmedian(RE_A)
m_CV_A=nanmedian(CV_A)
disp('mProphet')
for i=1:numel(mProphet_methods)
    disp(mProphet_methods{i})
    eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(commonQuants.' mProphet_methods{i} ',transpose(molarities));']);
end
for i=1:numel(mProphet_methods)
    accuracyGain(i,:)=abs((RE(i,:)))./abs((RE_A'))
    precisionGain(i,:)=abs((CV(i,:)))./abs((CV_A'))
end

boxplots(commonQuants,molarities,methods,mProphet_methods)

%% plot acc-prec on common peptides
h = figure;
plot(log10(molarities),(abs([RE_A';RE])),'o-')
hold on
plot(log10(molarities),([CV_A';CV]),'*--')
legend(methods)
set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
legend(methods, 'Location', 'southeast' )
xlabel('Dilution groups (fmol)')
ylabel('(%)')
title('Relative Error (solid line) & Quartile Coefficient of Dispersion (dashed line)')
pause
%     attachPlotToFile(h,'accPrec')
%     close(h)
plot(log10(molarities),(abs([RE])),'ko-','DisplayName','Skyline')
plot(log10(molarities),([CV]),'k*--','DisplayName','Skyline')
% to modify the legend
allDatah = flipud(get(gca,'children'));
% str = {'Ariadne','Skyline','totalxic','maxapex','apexsum','apexsum outlier'}
% size(allDatah)
str = {'Ariadne','Skyline'}
legend([allDatah(1); allDatah(size(allDatah,1)-1)],str)

% load('comparison.mat')
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

%% plot acc-prec
h = figure;
hold on
plot(log10(molarities),([CV_A';CV]),'*-')
set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
legend(methods, 'Location', 'southeast' )
xlabel('Dilution groups (fmol)')
ylabel('(%)')
title('Quartile Coefficient of Dispersion on common (dashed line) & additional (solid line) estimates')
break
%% commonly quantified peps - gains
commonQuants = estimatedAbundancesTmp;
commonQuants.Ariadne(isnan(commonQuants.Skyline_mProphet)) = nan;
commonQuants.Skyline_mProphet(isnan(commonQuants.Ariadne)) = nan;
%% correlation/consensus COMMON quants Ariadne and Skyline
corrValueCommon = getConsensus(commonQuants)
%complex = 0.8395
%noisy = 0.9887
disp('Ariadne')
[RE_A sd_A CV_A]=getAccPrec(commonQuants.Ariadne,molarities');
m_RE_A=nanmedian(RE_A)
m_CV_A=nanmedian(CV_A)
disp('mProphet')
for i=1:numel(mProphet_methods)
    disp(mProphet_methods{i})
    eval(['[RE(i,:) sd(i,:) CV(i,:)]=getAccPrec(commonQuants.' mProphet_methods{i} ',transpose(molarities));']);
end

    h = figure;
    plot(log10(molarities),(abs([RE_A';RE])),'o-')
    hold on
    plot(log10(molarities),([CV_A';CV]),'*:')
    legend(methods)
    set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))
    legend(methods, 'Location', 'southeast' )
    xlabel('Dilution groups (fmol)')
    ylabel('(%)')
    title('Relative Error (solid line) & Quartile Coefficient of Dispersion (dotted line)')

% %% completing graph
% 
% plot(log10(molarities),([CV_A';CV]),'+--')
% set(gca, 'XTick',sort(log10(molarities)),'XTickLabel',sort(molarities))











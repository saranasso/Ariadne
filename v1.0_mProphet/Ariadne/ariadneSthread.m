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
function [pepLibrary]=ariadneSthread(precursor, Q3s, rts,trans,pepLibrary,idx,sub_idx)
% this function performs the relative quantitation and consistently updates
% the pepLibrary
% idx - pep
% sub_idx - file

global settings

pepLibrary(idx).error=[];
if settings.mProphet
    try
        mRetTime=pepLibrary(idx).exp_RT_BestPeakGroup(sub_idx);
        if ~isempty(pepLibrary(idx).RT_time)
            expRetTime=pepLibrary(idx).RT_time;
        else
            expRetTime=max(min(rts.endo),min(rts.ref)); % just a value
        end
        mMaxApex.ref=pepLibrary(idx).abundance_mProphet(sub_idx).max_apex_ref;
        mMaxApex.endo=pepLibrary(idx).abundance_mProphet(sub_idx).max_apex_endo;
    catch
        if ~isempty(pepLibrary(idx).RT_time)
            mRetTime = pepLibrary(idx).RT_time;
            mMaxApex.ref = NaN;
            mMaxApex.endo = NaN;
        else
            mRetTime = NaN;
            mMaxApex.ref = NaN;
            mMaxApex.endo = NaN;
        end
    end
else
    expRetTime=pepLibrary(idx).RT_time;
end

% re-mapping (transposing) vectors
for i=1:size(trans.ref,2)
    endo{i}=trans.endo{i}';
    ref{i}=trans.ref{i}';
    rts_endo{i}=rts.endo';
    rts_ref{i}=rts.ref';
end

%% PEAK PICKING
if settings.mProphetFilter
    % using mProphet estimates and validation
    pepLibrary(idx).estRetTime(sub_idx)=mRetTime;
    pepLibrary(idx).estMaxApex(sub_idx)=mMaxApex;
    retTime = mRetTime;
    maxApex.ref = mMaxApex.ref*ones(size(trans.ref,2),1);
    maxApex.endo = mMaxApex.endo*ones(size(trans.ref,2),1); % ref and endo trans number MUST be the same since they are paired
    estimatedPeakWidth = settings.expPeakWidth;
elseif ~settings.mProphetFilter || isempty(mRetTime) || isempty(maxApex)
    % estimating which one is the best peak group for representing the
    % current peptide occurrence
    disp('peak picking')
    tic
    [retTime, maxApex, estimatedPeakWidth]=peakGroupSelection(rts_endo,endo,rts_ref,ref,pepLibrary(idx));
    toc
    if isempty(retTime) && isempty(maxApex) && isempty(estimatedPeakWidth)
        return
    end
    pepLibrary(idx).estRetTime(sub_idx)=retTime;
    pepLibrary(idx).estMaxApex(sub_idx)=maxApex;
end
% pepLibrary(idx).estError(sub_idx).endo=min(dist.endo);
% pepLibrary(idx).estError(sub_idx).ref=min(dist.ref);
% pepLibrary(idx).error=cat(2,pepLibrary(idx).error,[setdiff(dist.endo,min(dist.endo)),setdiff(dist.ref,min(dist.ref))]);

%% selecting region of interest (ROI) around the peak group
rti=retTime- estimatedPeakWidth;
rtf=retTime+ estimatedPeakWidth;
[endo ref rts_endo rts_ref] = getPeakGroupWindow(rti,rtf,rts,trans);

%% some needed checks on baseline signal
[endo, ref, rts_endo,rts_ref, rtsFilt, trans] = checkBaseline(endo,ref,rts_endo,rts_ref,trans,rts);
%
%% plot raw signal and estimates
plotRawTraces

%% SMOOTHING and plot smoothed signal in ROI
disp('baseline correction and smoothing')
[endo ref] = transProcessing (endo, ref, rts_endo, rts_ref);

plotProcessedTraces

%% force same retention times
[ref,endo,rts_endo,rts_ref] = setSameRetTimes(ref,endo,rts_endo,rts_ref);

% plotResampledTraces

%% computing cross correlation and "evaluating the lag"
for i=1:size(trans.ref,2)
    % [xcorrs(:,i) lags(:,i)]=xcorr(ref{i},endo{i},'unbiased'); 
    [xcorrs lags]=xcorr(ref{i}-mean(ref{i}),endo{i}-mean(endo{i}));
    [maxXcorr(i,1) idx_max]=max(xcorrs);
    lag=lags(idx_max);
end
% pepLibrary(idx).xCorr{sub_idx}=maxXcorr; 

warning on
%% wavelet identification of peaks in the trans'signal and correlation coefficients estimates
for i=1:size(trans.ref,2)
    %         figure
    %         subplot 211
    [Peaklist{i}.endo, PFWHH{i}.endo, PExt{i}.endo] =  mspeaks(rts_endo{i},endo{i},'HeightFilter',settings.minExpectPeakHeight,'OverSegmentationFilter', settings.expPeakWidth,'ShowPlot', 'false','Style','extline');
    %         subplot 212
    [Peaklist{i}.ref, PFWHH{i}.ref, PExt{i}.ref] =  mspeaks(rts_ref{i},ref{i},'HeightFilter',settings.minExpectPeakHeight,'OverSegmentationFilter', settings.expPeakWidth,'ShowPlot', 'false','Style','extline');
    %     pause
    [corr_coeff(:,:,i),pValue(:,:,i),corr_coeff_low(:,:,i),corr_coeff_up(:,:,i)]=corrcoef(ref{i},endo{i}); % at 0 lag! it penalizes deviation from co-elution
    
end

% store corr values for post-processing
pepLibrary(idx).corrCoef{sub_idx}=corr_coeff(1,2,:);

%% select best correlating trans pairs (if no pair is correlating then data rescue is executed)
idx_good_xCorr=find(corr_coeff(1,2,:)>settings.minXcorr);

if isempty(idx_good_xCorr) % if no trans pairs were highly correlated
    disp('data rescue')
    tic
    dataRescue
    toc
    if corr_coeff(1,2,1)>settings.minXcorr
        disp('corrected!')
        idx_good_xCorr=1;
    else
        return
    end
    if  pValue(1,2,idx_good_xCorr)>settings.xCorrSignificativity
        pepLibrary(idx).ratio(sub_idx)=NaN;
        return
    end
else % if some trans pairs were highly correlated
    for i=1:numel(idx_good_xCorr)
        if pValue(1,2,idx_good_xCorr(i))<settings.xCorrSignificativity
            significantOnes(i)=1;
        else
            significantOnes(i)=0;
        end
    end
    idx_good_xCorr=idx_good_xCorr(logical(significantOnes));
end

%% quantify only on best correlating trans pairs
% first, define peak borders through quantiles
% then, integrate the inlying ion counts

h=figure;
for j=1:numel(idx_good_xCorr)
    [idx_peak_ref meas]=getPeakIdx(Peaklist{idx_good_xCorr(j)}.ref,[retTime maxApex.ref(idx_good_xCorr(j))]);
    [idx_peak_endo meas]=getPeakIdx(Peaklist{idx_good_xCorr(j)}.endo,[retTime maxApex.endo(idx_good_xCorr(j))]);
    idx_peak.endo=idx_peak_endo;
    idx_peak.ref=idx_peak_ref;
    rts_tmp{idx_good_xCorr(j)}.endo=rts_endo{idx_good_xCorr(j)};
    rts_tmp{idx_good_xCorr(j)}.ref=rts_ref{idx_good_xCorr(j)};
    trans_tmp.endo=endo;
    trans_tmp.ref=ref;
    trans_abundance(j)=getPairedTransAbundance(rts_tmp{idx_good_xCorr(j)},trans_tmp,PExt{idx_good_xCorr(j)},idx_peak,idx_good_xCorr(j),estimatedPeakWidth);
end
title({['k is ' num2str(idx) ' and i is ' num2str(sub_idx)];['Peptide Sequence: ' pepLibrary(idx).sequence];['File name: ' strrep(pepLibrary(idx).dataFile{sub_idx},'_',' ')];['Corr values are: ' num2str(pepLibrary(idx).corrCoef{sub_idx})]})
attachPlotToFile(h,'Ariadne_plots')
close all
java.lang.System.gc();
transRatio=[trans_abundance(:).endo]./[trans_abundance(:).ref];

%% computing the pep overall abundance and weighted ratio
weights=[trans_abundance.endo]+[trans_abundance.ref];
for j=1:numel(idx_good_xCorr)
    wweights(j)=weights(j)*corr_coeff(1,2,idx_good_xCorr(j))*maxXcorr(idx_good_xCorr(j));
end
pepRatio=sum(wweights.*transRatio)/sum(wweights);
pepAbundance.endo=(prod([trans_abundance.endo])).^(1/numel(idx_good_xCorr));
pepAbundance.ref=(prod([trans_abundance.ref])).^(1/numel(idx_good_xCorr));
wpepAbundance.endo=(prod([trans_abundance.endo].*wweights)).^(1/numel(idx_good_xCorr));
wpepAbundance.ref=(prod([trans_abundance.ref].*wweights)).^(1/numel(idx_good_xCorr));

%% updating pepLibrary with quant estimates
pepLibrary(idx).ratio(sub_idx)=pepRatio;
pepLibrary(idx).abundance(sub_idx).endo=pepAbundance.endo;
pepLibrary(idx).abundance(sub_idx).ref=pepAbundance.ref;
pepLibrary(idx).wabundance(sub_idx).endo=wpepAbundance.endo;
pepLibrary(idx).wabundance(sub_idx).ref=wpepAbundance.ref;
disp('another ratio computed! ')


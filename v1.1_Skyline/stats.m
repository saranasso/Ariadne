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


estimatedAbundancesTmp=estimatedAbundances;
whichstatsOption = '{''numel'',''range''}';
plotOn=1; % to plot graphs
logScale = 1; % to evaluate stats on the log scale
logScaleOutRemoval = 0; % deprecated
removeOutliers = 1; % to remove outliers
pause off
for numMethod=1:numel(methods)
    
    eval(['[x(:,numMethod) g] = straightenX((estimatedAbundancesTmp.' methods{numMethod} '),molarities);']); % g is the same for all methods
    eval(['range.' methods{numMethod} ' = grpstats(x(:,numMethod), g,''range'');']);
    estimatedAbundancesTmp.Ariadne(isinf(estimatedAbundancesTmp.Ariadne))=nan;
    for i=1:numel(molarities)
        if plotOn
            histfit(eval(['estimatedAbundancesTmp.' methods{numMethod} '(:,i)']))
            h=gcf;
            title('Pre outliers removal')
            attachPlotToFile(h,['Normality assesment - qqplot for ' methods{numMethod}])
            pause(0.5)
            close(h)
        end
        if removeOutliers
            if logScale
                eval(['[estimatedAbundancesTmp.' methods{numMethod} '(:,i)] = rmOutliersOutOf1pt5IQR(log10(estimatedAbundancesTmp.' methods{numMethod} '(:,i)),0);']);
            else
                if logScaleOutRemoval
                    eval(['[estimatedAbundancesTmp.' methods{numMethod} '(:,i)] = rmOutliersOutOf1pt5IQR((estimatedAbundancesTmp.' methods{numMethod} '(:,i)),1);']);
                else
                    eval(['[estimatedAbundancesTmp.' methods{numMethod} '(:,i)] = rmOutliersOutOf1pt5IQR((estimatedAbundancesTmp.' methods{numMethod} '(:,i)),0);']);
                end
            end
            if plotOn
                histfit(eval(['estimatedAbundancesTmp.' methods{numMethod} '(:,i)']))
                pause(0.5)
                h=gcf;
                title('Post outliers removal')
                attachPlotToFile(h,['Normality assesment - qqplot for ' methods{numMethod}])
                close(h)
            end
        else logScale && ~removeOutliers
            eval(['estimatedAbundancesTmp.' methods{numMethod} '(:,i) = log10(estimatedAbundancesTmp.' methods{numMethod} '(:,i));']);
        end
        if plotOn
            h=figure
            %         bootstrap(100, @qqplot, eval(['estimatedAbundancesTmp.' methods{numMethod} '(:,i)']))
            qqplot(eval(['estimatedAbundancesTmp.' methods{numMethod} '(:,i)']))
            pause(0.5)
            h=gcf;
            attachPlotToFile(h,['Normality assesment - qqplot for ' methods{numMethod}])
            close(h)
        end
    end
    eval(['[x(:,numMethod) g] = straightenX((estimatedAbundancesTmp.' methods{numMethod} '),molarities);']); % g is the same for all methods
    
end
X_SPSS=x;
X_SPSS(isnan(X_SPSS)) = 10000; % matrix to be used in IBM SPSS to get the statistical analysis results as published


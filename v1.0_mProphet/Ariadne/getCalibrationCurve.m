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
function [regression]=getCalibrationCurve(abundance,concentration)
% it checks for linear relationship and fit the line

% check whether there is any linear relationship
[r p] = corr(abundance,concentration,'tail','gt','rows','pairwise');

if isnan(r)
    disp('r is nan!')
end

if p<0.05 & r^2>0.9 % if there is linear relationship
    
    N_dati=numel(concentration);
    
    axis([0,max(abundance),0,max(concentration)])
    axis manual
    axis(axis)
    hold on
    plot(abundance,concentration,'o')
    titolo=[{'Calibration Curve'}];
    title(titolo)
    x_SD=sqrt(abundance);
    y_SD=sqrt(concentration);
    
    [m,m_SE,q,q_SE,min_quad_distXY,residual,covar_post,J,SE]=fitlineXY(abundance,x_SD,concentration,y_SD,1,0);
    hold on
    x=(0:0.1:max(abundance));
    plot(x,m*x+q,'r')
    
    %% stats
    yhat=m*abundance+q;
    xhat=(concentration-q)/m;
    regression.m=m;
    regression.q=q;
    regression.m_SE=m_SE;
    regression.q_SE=q_SE;
    regression.R2=r^2;
    regression.pValue_R2=p;
    regression.RMSE=SE;
    regression.minEstimate=yhat(end);
    regression.maxEstimate=yhat(1);
    regression.localDynamicRange=log10(regression.maxEstimate/regression.minEstimate);


    legend({['Quantifications #=' num2str(N_dati) ' R2xy=' num2str(r^2,'%.2f') ' pValue=' num2str(p)],['Regression line: m=' num2str(m,'%.2f') ' q=' num2str(q,'%.2f')  ' SE(m)=' num2str(m_SE,'%.2f') ' SE(q)=' num2str(q_SE,'%.2f') ' RMSE=' num2str(SE,'%.2f')]})
    xlabel('estimate')
    ylabel('concentration')
    hold off
else
    regression=[];
    disp(['R2 is: ' num2str(r^2) ' and its p-value is: ' num2str(p) ' --> no linear relationship!'])
    return
end




	

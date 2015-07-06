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
function [m,m_SE,q,q_SE,min_quad_distXY,residual,covar_post,J,RMSE]=fitlineXY(x,x_SD,y,y_SD,m_start,q_start)
% it performs the actual fit

options=optimset('TolX',1e-10,'TolFun',1e-10);
x_apex=x./x_SD;
y_apex=y./y_SD;
start=[m_start,q_start];
parameters=[x_apex,x_SD,y_apex,y_SD];

[finish,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@(start) distanceXY(start,parameters),start,[],[],options);

min_quad_distXY=resnorm;
m=finish(1);
q=finish(2);
J=full(jacobian);   % lsqnonlin returns weighted Jacobian 
covpar=inv(J'*J);
dataNumber=length(residual);
paramNumber=length(start);
squared_sigma_post=resnorm/(dataNumber-paramNumber);
covar_post=squared_sigma_post.*covpar;
var_post=diag(covar_post);
SE_post=sqrt(var_post);
RMSE=sqrt(resnorm/(dataNumber-paramNumber));
m_SE=SE_post(1);
q_SE=SE_post(2);
end


	

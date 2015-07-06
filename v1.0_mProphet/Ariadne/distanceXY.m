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
function [res_XY]=distanceXY(variables, parameters)
% it computes residuals

m=variables(1);
q=variables(2);

x_apex=parameters(:,1);
x_SD=parameters(:,2);
y_apex=parameters(:,3);
y_SD=parameters(:,4);

% computing transformed coefficients
m_apex=m.*x_SD./y_SD;
q_apex=q./y_SD;

% applying residues formula
res_XY=sqrt((y_apex-(m_apex.*x_apex)-q_apex).^2./(1+m_apex.^2));	
	

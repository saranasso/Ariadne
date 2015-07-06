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

function input = getBooleanInput(prompt,dlg_title,num_lines,def)
% it keeps on asking for a boolean input value till it gets it

answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

input=str2double(answer{1});
while isempty(intersect(input,[0,1]))
    disp('Only 0 or 1 accepted!')
    input = getBooleanInput(prompt,dlg_title,num_lines,def);
end



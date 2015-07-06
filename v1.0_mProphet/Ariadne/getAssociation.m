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

function association_file = getAssociation(mzXMLFilesPath)
% it asks the user for the association "file to group"

global settings

for i=1:numel(mzXMLFilesPath)
    if settings.dilution_series
        prompt = {['Enter the real value to be estimated corresponding to data file ' mzXMLFilesPath{i}]};
        dlg_title = 'Input info';
        num_lines = 1;
        def = {'e.g., 0.5 for a file named 20120507_AQ1_05ftmol_1_rep2.mzXML'};
    else
        prompt = {['Enter an identifier for data file ' mzXMLFilesPath{i} ' (pep quants on files with the same identifier will be averaged)']};
        dlg_title = 'Input info';
        num_lines = 1;
        def = {'e.g.,  for a file named 20120507_AQ1_05ftmol_1_rep2.mzXML'};
    end
    try
        value = getNumericalInput(prompt,dlg_title,num_lines,def);
    catch
        disp('Only one numerical value accepted!')
    end
end
[sorted_dilutions idx_dil]=sort(value,2,'descend');
mzXMLFilesPath_trunc=strtok(mzXMLFilesPath(:),'.');
association_file=[num2cell(sorted_dilutions)', mzXMLFilesPath_trunc(idx_dil)];
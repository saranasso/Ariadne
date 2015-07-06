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
GPLdisclaimer=importdata('/Users/saran/work/SRM data/codice/GNU_GPL.txt');
dataDir = uigetdir(pwd, 'Select the directory where to add the license');
cd(dataDir)
fileM=struct2cell(dir('*.m'));
mFilesPath=fileM(1,:);
for i=1:numel(mFilesPath)
    filePath=[pwd '/' mFilesPath{i}];
    fid=fopen(filePath);
    [file]=fscanf(fid,'%c');
    data=cat(1,GPLdisclaimer,file);
    writeIt(filePath,data)
end	

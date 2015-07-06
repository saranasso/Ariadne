%----------------------------
% Standardize x and g.
% x can be a numeric vector or matrix.  It is guaranteed to be non-empty.
% g can be a numeric vector or matrix, categorical vector or matrix, char
%   array, cellstr, or a cell vector of vectors or cellstrs.  The vectors
%   should be the same length as x, if x is a column, or the same as the
%   number of columns in x if x is a matrix.
% Return xDat as a numeric column vector, with xlen rows.
% Return gDat as a non-empty cell vector of column vectors, each with
%   uniform length of xlen.  The column vectors are not necessarily
%   numeric.
% Return origRow, a column xlen long containing the row the observation
%   came from.
% Return xlen, the height of xDat.
% Return gexplicit, which is true if the g parameter is non-empty, due to
%   being specified by the user.

function [xDat,gDat,origRow,xlen,gexplicit,origInd,origNumXCols] =...
    straightenX(x,g)


if ~isnumeric(x) || ndims(x)>2
    error(message('stats:boxplot:XNumeric'));
end
[xRows,xCols] = size(x);

xDat = x(:);
xlen = numel(xDat);
origInd = (1:xlen)';
if xRows==1
    % If X is a vector, consider it to be a column vector.
    xRows = xlen;
    xCols = 1;
end
origNumXCols = xCols;

[gDat, gRows, gCols] = convertToCellarrayOfColumnVectors(g,'G');

if xCols==1
    % X is a vector.
    
    origRow = (1:xlen)';
    
    if gRows==0
        % Treat X as one big group, rather than xlen groups each with
        % one element.
        gexplicit = false;
        gDat = {ones(xlen,1)};
    else
        gexplicit = true;
        % If x was specified as a vector, gRows may be 1 or xlen.
        if ~(gRows==1 || gRows==xlen)
            error(message('stats:boxplot:XGLengthMismatch'));
        end
        % If gDat has only 1 row, copy the one group repeatedly
        % to expand it to xlen rows.
        if gRows==1
            for i=1:gCols
                gDat{i} = gDat{i}(ones(xlen,1));
            end
        end
        
    end
    
else
    % When X is a matrix, convert it to a vector.  Each column of X will be
    % treated as a distinct group.
    
    [origRow,xMatgroup] = ind2sub([xRows,xCols],(1:xlen)');
    
    if gRows==0
        % If G absent, assign G according to X matrix column.
        gexplicit = false;
        gDat = {xMatgroup};
    else       
        % If G is present, it is required to have one element (for one
        % grouping variable in G) or row (for multiple grouping variables)
        % per column of X or per element of X.
        %
        % The number of columns in X will often be the same as
        % numFlatGroups, when it is calculated later, but it may differ
        % depending on the contents of G and the fullfactors parameter...
        % duplicates values and NaN/'' in G will reduce the eventual value
        % of numFlatGroups, while fullfactors will increase it if G
        % contains multiple grouping variables and does not contain every
        % possible combination of factor levels.
        
        gexplicit = true;
        if gRows==xlen
            % Do nothing.
        elseif gRows==xCols
            % Expand vectors in gDat to give them xlen rows.
            for i=1:gCols
                gDat{i} = gDat{i}(xMatgroup);
            end
        else
            error(message('stats:boxplot:XGLengthMismatch'));
            
        end
    end
end




end


%----------------------------
% Standardize large args, and verify that lengths are self-consistent.
%
% Accepts argument in various forms:
% vector (numeric or categorical)
% matrix (numeric or categorical)
% cell vector of scalars or strings
% cell vector of vectors (vectors may be numeric, cellstrs, or categorical)
% character array
% empty
%
% Returns argument out in a standardized form:
% cell vector of equal-length column vectors
% out is empty if in is empty.

function [out,numrows,numcols] = ...
    convertToCellarrayOfColumnVectors(in,argname,outputType)

if nargin<3
    outputType = 'any';
end

if iscell(in) && ~isvector(in) && ~isempty(in)
    error(message('stats:boxplot:BadCellDataInput', argname));
end

% Empty [] -> return empty.
if isempty(in)
    out = {};
    numcols = 0;
    numrows = 0;
    return;
end

% Convert in into a cell vector of cells or a cell with one vector.
if ~iscell(in)
    in = {in};
else
    
    % Elements must be all cells, all scalars, or all chars.
    if ischar(in{1}) && size(in{1},1)<=1 && isvector(in)
        % Check that we have a cell vector of chars.
        for i=2:length(in)
            if ~ischar(in{i}) || size(in{1},1)>1
                error(message('stats:boxplot:AllCellsAllScalars', argname));
            end
        end
        in = cellstr(in); % Create a cell vector with one cellstr.
        in = {in(:)};
    elseif isscalar(in{1}) && ~iscell(in{1}) && isvector(in)
        % Check that we have a cell vector of scalars.
        for i=2:length(in)
            if ~isscalar(in{i}) || iscell(in{i})
                error(message('stats:boxplot:AllCellsAllScalars', argname));
            end
        end
        tmp = in{1};
        for i=2:length(in)
            tmp(i) = in{i};
        end
        in = {tmp(:)};% Create cell with one vector.
    elseif iscell(in{1}) || ndims(in{1})==2
        % Check that we have a cell vector that contains cells, vectors, or
        % matrices.
        for i=2:length(in)
            if ~iscell(in{i}) && ndims(in{i})~=2
                error(message('stats:boxplot:AllCellsAllScalars', argname));
            end
            % Do nothing, we have a cell vector of cells.
        end
    end
end

% Create out from in, making it a row cell vector of column vectors.
% Initialize out, this may grow if in has matrix elements.
out = cell(1,length(in));
j = 1;
for i=1:length(in)
    % Unpack one level if needed.
    if isscalar(in{i}) && iscell(in{i}) && ~iscellstr(in{i})
        in{i}=in{i}{1};
    end
    if isvector(in{i}) && ~ischar(in{i})% Includes vector numeric and cellstr.
        out{j} = in{i}(:);
        j = j+1;
    elseif ischar(in{i}) % Char array, convert to cellstr.
        out{j} = cellstr(in{i});
        j = j+1;
    elseif ndims(in{i})==2 && ~iscell(in{i}) % Matrix, split up columns.
        for k=1:size(in{i},2)
            out{j} = in{i}(:,k);
            j = j+1;
        end
    else
        error(message('stats:boxplot:ArgInvalid', argname));
    end
end


numcols = size(out,2);
numrows = size(out{1},1);
% Test that all vectors are the same length.
for i=2:numcols
    if size(out{i},1)~=numrows
        error(message('stats:boxplot:UnequalLength', argname));
    end
end

% Cast or complain about variable types.
switch outputType
    case 'numeric'
        for i=1:numcols
            if ~isnumeric(out{i})
                error(message('stats:boxplot:ArgNonNumeric', argname));
                
            end
        end
    case 'string'
        for i=1:numcols
            if isnumeric(out{i})
                out{i} = cellstr(num2str(out{i}));
            end
        end
    case 'any'
        % Omit checking or casting.
    otherwise
        error(message('stats:boxplot:badOutputtype'));
end

end	

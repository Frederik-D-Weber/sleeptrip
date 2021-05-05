function [Y,ndx,dbg] = natsortfiles(X,rgx,varargin)
% Natural-order / alphanumeric sort of filenames or foldernames.
%
% (c) 2014-2021 Stephen Cobeldick
%
% Sorts text by character code and by number value. File/folder names, file
% extensions, and path directories (if supplied) are sorted separately to
% ensure that shorter names sort before longer names. For names without
% file extensions (i.e. filenames without extensions, or foldernames) use
% the 'noext' option. Use the 'nopath' option to ignore any filepath.
%
%%% Example:
% P = 'C:\SomeDir\SubDir';
% S = dir(fullfile(P,'*.txt'));
% S = natsortfiles(S);
% for k = 1:numel(S)
%     F = fullfile(P,S(k).name)
% end
%
%%% Syntax:
%  Y = natsortfiles(X)
%  Y = natsortfiles(X,rgx)
%  Y = natsortfiles(X,rgx,<options>)
% [Y,ndx,dbg] = natsortfiles(X,...)
%
% To sort the elements of a string/cell array use NATSORT (File Exchange 34464)
% To sort the rows of a string/cell array use NATSORTROWS (File Exchange 47433)
%
%% File Dependency %%
%
% NATSORTFILES requires the function NATSORT (File Exchange 34464). The optional
% arguments <options> are passed directly to NATSORT. See NATSORT for case
% sensitivity, sort direction, number substring matching, and other options.
%
%% Explanation %%
%
% Using SORT on filenames will sort any of char(0:45), including the printing
% characters ' !"#$%&''()*+,-', before the file extension separator character '.'.
% Therefore this function splits the name and extension and sorts them separately.
%
% Similarly the path separator character within filepaths can cause longer
% directory names to sort before shorter ones, as char(0:46)<'/' and
% char(0:91)<'\'. This example on a PC demonstrates why this matters:
%
% >> X = {'A1\B', 'A+/B', 'A/B1', 'A=/B', 'A\B0'};
% >> sort(X)
% ans =   'A+/B'  'A/B1'  'A1\B'  'A=/B'  'A\B0'
% >> natsortfiles(X)
% ans =   'A\B0'  'A/B1'  'A1\B'  'A+/B'  'A=/B'
%
% NATSORTFILES splits filepaths at each path separator character and sorts
% every level of the directory hierarchy separately, ensuring that shorter
% directory names sort before longer, regardless of the characters in the names.
% On a PC separators are '/' and '\' characters, on Mac and Linux '/' only.
%
%% Examples %%
%
% >> A = {'a2.txt', 'a10.txt', 'a1.txt'}
% >> sort(A)
% ans = 'a1.txt'  'a10.txt'  'a2.txt'
% >> natsortfiles(A)
% ans = 'a1.txt'  'a2.txt'  'a10.txt'
%
% >> B = {'test_new.m'; 'test-old.m'; 'test.m'};
% >> sort(B) % Note '-' sorts before '.':
% ans =
%    'test-old.m'
%    'test.m'
%    'test_new.m'
% >> natsortfiles(B) % Shorter names before longer:
% ans =
%    'test.m'
%    'test-old.m'
%    'test_new.m'
%
% >> C = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
% >> sort(C) % Wrong number order:
% ans =
%    'test.m'
%    'test1.m'
%    'test10-old.m'
%    'test10.m'
%    'test2.m'
% >> natsortfiles(C) % Shorter names before longer:
% ans =
%    'test.m'
%    'test1.m'
%    'test2.m'
%    'test10.m'
%    'test10-old.m'
%
%%% Directory Names:
% >> D = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1archive.zip';'A1\test.m'};
% >> sort(D) % Wrong number order, and '-' sorts before '\':
% ans =
%    'A10\test.m'
%    'A1\test.m'
%    'A1archive.zip'
%    'A2-old\test.m'
%    'A2\test.m'
% >> natsortfiles(D) % Shorter names before longer:
% ans =
%    'A1archive.zip'
%    'A1\test.m'
%    'A2\test.m'
%    'A2-old\test.m'
%    'A10\test.m'
%
%% Input and Output Arguments %%
%
%%% Inputs (**=default):
% X   = Array of file or folder names to be sorted. Can be the struct
%       returned by DIR, a string array, or a cell array of char vectors.
% rgx = Regular expression to match number substrings, '\d+'**
%     = [] uses the default regular expression.
% <options> can be supplied in any order:
%     = 'noext' should be specified for names without extensions.
%     = 'nopath' sorts by name only, ignoring any preceding path.
%     = all remaining options are passed directly to NATSORT.
%
%%% Outputs:
% Y   = Array X sorted into alphanumeric order.
% ndx = NumericMatrix, same size as X. Indices such that Y = X(ndx).
% dbg = CellVectorOfCellArrays, size 1xMAX(2+NumberOfDirectoryLevels).
%       Each cell contains the debug cell array for directory names, filenames,
%       and file extensions. Helps debug the regular expression (see NATSORT).
%
% See also SORT NATSORT NATSORTROWS DIR FILEPARTS FULLFILE NEXTNAME CELLSTR REGEXP IREGEXP SSCANF

%% Input Wrangling %%
%
fun = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
%
if isstruct(X)
	assert(isfield(X,'name'),...
		'SC:natsortfiles:X:StructMissingNameField',...
		'If first input <X> is a struct then it must have field <name>.')
	fnm = {X.name};
	assert(all(fun(fnm)),...
		'SC:natsortfiles:X:StructNameInvalidType',...
		'First input <X> field <name> must contain only character vectors.')
	[pth,fnm,ext] = cellfun(@fileparts, fnm, 'uni',false);
	if isfield(X,'folder')
		pth = {X.folder};
		assert(all(fun(pth)),...
			'SC:natsortfiles:X:StructFolderInvalidType',...
			'First input <X> field <folder> must contain only character vectors.')
	end
elseif iscell(X)
	assert(all(fun(X(:))),...
		'SC:natsortfiles:X:CellContentInvalidType',...
		'First input <X> cell array must contain only character vectors.')
	[pth,fnm,ext] = cellfun(@fileparts, X(:), 'uni',false);
elseif ischar(X)
	[pth,fnm,ext] = cellfun(@fileparts, cellstr(X), 'uni',false);
else
	try
		isg = isstring(X); % ISSTRING introduced R2016b.
	catch
		isg = false;
	end
	assert(isg,...
		'SC:natsortfiles:X:InvalidType',...
		'First input <X> must be a structure, a cell array, or a string array.');
	[pth,fnm,ext] = cellfun(@fileparts, cellstr(X(:)), 'uni',false);
end
%
varargin = cellfun(@nsStr2Char, varargin, 'uni',false);
%
assert(all(fun(varargin)),...
	'SC:natsortfiles:option:InvalidType',...
	'All optional arguments must be character vectors or string scalars.')
%
ide = strcmpi(varargin,'noext');
assert(nnz(ide)<2,...
	'SC:natsortfiles:noext:Overspecified',...
	'File-extension handling is overspecified.')
varargin(ide) = [];
%
idp = strcmpi(varargin,'nopath');
assert(nnz(idp)<2,...
	'SC:natsortfiles:nopath:Overspecified',...
	'File-path handling is overspecified.')
varargin(idp) = [];
%s
if nargin>1
	varargin = [{rgx},varargin];
end
%
%% Path and Extension %%
%
% Path separator regular expression:
if ispc()
	rgx = '[^/\\]+';
else % Mac & Linux
	rgx = '[^/]+';
end
%
% No file extension:
if any(ide)
	fnm = strcat(fnm,ext);
	ext(:) = {''};
end
%
% No filepath:
if any(idp)
	num = 0;
else
	% Split path into {dir,subdir,subsubdir,...}:
	spl = regexp(pth,rgx,'match');
	len = cellfun('length',spl);
	num = max(len);
	vec = cell(numel(len),1);
end
%
%% Sort File Names/Paths %%
%
% Alphanumeric sort of the file extensions and filenames:
if isempty(num)
	ndx = [];
	ids = [];
	dbg = {};
elseif nargout<3 % faster:
	[~,ndx] = natsort(ext,varargin{:});
	[~,ids] = natsort(fnm(ndx),varargin{:});
else % for debugging:
	[~,ndx,dbg{num+2}] = natsort(ext,varargin{:});
	[~,ids,tmp] = natsort(fnm(ndx),varargin{:});
	[~,idd] = sort(ndx);
	dbg{num+1} = tmp(idd,:);
end
ndx = ndx(ids);
%
% Alphanumeric sort of the directory names:
for k = num:-1:1
	idx = len>=k;
	vec(:) = {''};
	vec(idx) = cellfun(@(c)c(k),spl(idx));
	if nargout<3 % faster:
		[~,ids] = natsort(vec(ndx),varargin{:});
	else % for debugging:
		[~,ids,tmp] = natsort(vec(ndx),varargin{:});
		[~,idd] = sort(ndx);
		dbg{k} = tmp(idd,:);
	end
	ndx = ndx(ids);
end
%
% Return the sorted input array and corresponding indices:
if ischar(X)
	ndx = ndx(:);
	Y = X(ndx,:);
else
	ndx = reshape(ndx,size(X));
	Y = X(ndx);
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles
function txt = nsStr2Char(txt)
% If scalar string extract the character vector, otherwise data unchanged.
try
	isg = isstring(txt); % ISSTRING introduced R2016b.
catch
	isg = false;
end
if isg && isscalar(txt)
	txt = txt{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsStr2Char
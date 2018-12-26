function mesh = rectMesh(domain, shape)
%RECTMESH	construct a rectangular mesh.
% mesh = RECTMESH(domain, shape) returns a 2D rectangle mesh defined by
% {(x,y) | domain(1)?x?domain(2), domain(3)?y?domain(4)}. The mesh is
% splitted into shape(1) cells in x direction and shape(2) cells in y 
% direction.
% 
% mesh = RECTMESH(domain, h) returns the same 2D rectangle mesh but with a
% cell size h(1) in x direction and h(2) in y direction.
%
% Author, Xiukun Hu, xiukun.hu@outlook.com
%
% See also Mesh

if numel(shape) == 1
	shape = [shape shape];
end

left = domain(1);
right = domain(2);
bottom = domain(3);
top = domain(4);

nx = shape(1);
ny = shape(2);

nnodex = nx+1;
nnodey = ny+1;

x = linspace(left, right, nnodex);
y = linspace(bottom, top, nnodey);

P = gen_pb(x, y);
T = gen_tb(nnodex, nnodey, 1);
mesh = MatFem.Mesh(P, T);

end


function Pb = gen_pb(x, y)
% Given coordinates
%	x = [x1 x2 x3 ...],
%	y = [y1 y2 y3 ...],
% returns
%	Pb = [x1 x1 x1 ... x2 x2 x2 ... x3 x3 x3 ... ...;
%		  y1 y2 y3 ... y1 y2 y3 ... y1 y2 y3 ... ...]

[X, Y] = meshgrid(x, y);
Pb = [X(:), Y(:)];
end

function Tb = gen_tb(nx, ny, basisType)
% returns matrix Tb as
%	Tb = [node 1 index array;
%		  node 2 index array;
%		  node 3 index array;
%		  node 4 index array if basis type > 1;
%		  ...]
% where ith column refers to all the nodes of element i.
%
% e.g., a finite element with basisType = 2:
%
%	5--10--15--20--25
%	| \  IV | \VIII |
%	4   9  14  19  24
%	| III \ | VII \ |
%	3---8--13--18--23
%	| \  II | \  VI |
%	2   7  12  17  22
%	|  I  \ |  V  \ |
%	1---6--11--16--21
%
% should give
%
%	Tb = [ 1  3  3  5 11 13 13 15;
%		  11 11 13 13 21 21 23 23;
%		   3 13  5 15 13 23 15 25;
%		   6  7  8  9 16 17 18 19;
%		   7 12  9 14 17 22 19 24;
%		   2  8  4 10 12 18 14 20]

nNodes = nx*ny;

idx = reshape(1:nNodes, ny, nx);

% Consider the first square
% NW---N---NE
% |	 \      |
% |   \     |
% W    C    E
% |	     \  |
% |       \ |
% SW---S---SE

SW = idx(1,1);
NW = idx(1+basisType, 1);
SE = idx(1, 1+basisType);
NE = idx(1+basisType, 1+basisType);

firstSqr = [SW SE NW
			NW SE NE];

switch basisType
	case 1
		% do nothing
	case 2
		W = idx(2,1);
		E = idx(2,3);
		S = idx(1,2);
		N = idx(3,2);
		C = idx(2,2);
		
		append = [S C W
				  C E N];
		firstSqr = [firstSqr append];
			
	otherwise
		error('MatFem:Unsupported', ['Function Space with basis type', ...
			'%d is not supported.'], basisType);
end

% Note: Each square is a repeat of the first square with a shift of index.

nySqr = (ny-1)/basisType;
nxSqr = (nx-1)/basisType;

% repeat and increment first square to get first column
firstCol = rep_inc(firstSqr, nySqr, basisType);
% repeat and increment first column to get Tb
Tb = rep_inc(firstCol, nxSqr, ny*basisType);
end

function r = rep_inc(A, nrepeat, inc)
% repeat matrix and increment each repeat with inc.
r = repmat(A, nrepeat, 1);
toAdd = 0:inc:(nrepeat-1)*inc;
toAdd = repelem(toAdd', size(A,1), size(A,2));
r = r+toAdd;
end
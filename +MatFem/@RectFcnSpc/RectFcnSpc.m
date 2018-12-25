classdef RectFcnSpc < handle
	% RectFcnSpc constructs a rectangular function space.
	%	Create a 2D rectangular function space with triangle elements.
	%	The space is a Continuous Galerkin's finite element space with
	%	order = 1 or 2, set by the argument basisType.
	
	properties (Access = private)
		qpts	% quadrature points
	end
	
	properties (SetAccess = private)
		mesh		% the mesh (rectangular) this function space builds on
		basisType	% basis type, 1 = linear, 2 = quadratic
		
		% # of FE nodes in x direction, different from mesh.nx if basis type
		% is not linear
		nx
		
		% # of FE nodes in y direction, different from mesh.ny if basis type
		% is not linear
		ny
		
		x	% array of x coordinates of FE nodes
		y	% array of y coordinates of FE nodes
		
		% node coordinate matrix, Pb(:, i) = [xi; yi] where (xi, yi) is the
		% coordinate of the ith node.
		Pb
		
		% element matrix, Tb(:, i) = [i1; i2; i3; ...] where i's are
		% indices of nodes in ith element.
		Tb
		
	end
	
	properties (Dependent)
		nLocalBases		% # of local bases in each element
		nBases			% # of total bases
		nElems			% # of total elements
	end
	
	methods 
		function obj = RectFcnSpc(rectMesh, basisType)
			% RectFcnSpc constructs a rectangular function space on mesh
			% rectMesh and with linear bases if basisType = 1, and
			% quadratic bases if basisType = 2.
			
			obj.mesh = rectMesh;
			obj.basisType = basisType;
			
			obj.nx = (rectMesh.nx-1)*basisType + 1;
			obj.ny = (rectMesh.ny-1)*basisType + 1;
			
			obj.x = linspace(rectMesh.left, rectMesh.right, obj.nx);
			obj.y = linspace(rectMesh.bottom, rectMesh.top, obj.ny);
			
			obj.Pb = gen_pb(obj.x, obj.y);
			obj.Tb = gen_tb(obj.nx, obj.ny, basisType);
			
			obj.qpts = calc_qpts(obj.mesh.hx, obj.mesh.hy);
		end
		
		function [X, Y] = meshgrid(spc)
			% meshgrid returns the x and y coordinates of all nodes, in a
			% matrix.
			
			[X, Y] = meshgrid(spc.x, spc.y);
		end
		
		function [vtc, idx] = elem(spc, i)
			idx = spc.Tb(:,i);
			vtc = spc.Pb(:,idx)';
		end
		
		function r = get.nLocalBases(spc)
			r = 3*spc.basisType;
		end
		
		function r = get.nElems(spc)
			r = size(spc.Tb, 2);
		end
		
		function r = get.nBases(spc)
			r = size(spc.Pb, 2);
		end
		
		function iE = ptInElem(spc, x, y)
			% Determine the index of the element point [x, y] is in
			m = spc.mesh;
			ixx = x/m.hx;
			iyy = y/m.hy;
			ix = floor(ixx);
			iy = floor(iyy);
			ixr = ixx - ix;
			iyr = iyy - iy;
			ithSqr = ix*(m.ny-1)+iy+1;
			iE = ithSqr*2 - (iyr < (1-ixr));
		end
		
		r = eval(spc, x, y, der)
				
		function r = intOverElem(spc, iE, ders, coefFcn, nqpts)
			vE = spc.elem(iE);
			vE = vE(1:3,:);
			q = spc.getQuadPts(iE, vE);
			r = MatFem.intOver(vE, spc.basisType, ders, coefFcn, q);
		end
		
		function r = intOverEdge(spc, edge, ders, coefFcn, nqpts)
			vE = spc.Pb(:, edge)';
			[X, Y, Wx, Wy] = MatFem.quadPts(nqpts, vE(1:2,:));
			r = MatFem.intOver(vE, spc.basisType, ders, coefFcn, {X, Y, Wx, Wy});
		end
		
		function r = evalFcn(spc, u, X, Y, iE, der)
			nls = spc.nLocalBases;
			vE = spc.elem(iE);
			
			r = 0;
			for ib = 1:nls
				iv = spc.Tb(ib, iE);
				r = r + u(iv)*MatFem.evalTriElem(X, Y, vE, ib, spc.basisType, der);
			end
		end
		
		function r = getQuadPts(spc, iE, vE)
			switch mod(iE, 2)
				case 0
					r = spc.qpts([3 4 5 6]);
				case 1
					r = spc.qpts([1 2 5 6]);
				otherwise
					error('iE must be a integer');
			end
			r{1} = r{1} + vE(1,1);
			r{2} = r{2} + vE(1,2);
		end
	end
	
end


function Pb = gen_pb(x, y)
% Given coordinates
%	x = [x1 x2 x3 ...],
%	y = [y1 y2 y3 ...],
% returns
%	Pb = [x1 x1 x1 ... x2 x2 x2 ... x3 x3 x3 ... ...;
%		  y1 y2 y3 ... y1 y2 y3 ... y1 y2 y3 ... ...]
	
	[X, Y] = meshgrid(x, y);
	Pb = [X(:), Y(:)]';
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

	firstSqr = [SW NW;
				SE SE;
				NW NE];

	switch basisType
		case 1
			% do nothing
		case 2
			W = idx(2,1);
			E = idx(2,3);
			S = idx(1,2);
			N = idx(3,2);
			C = idx(2,2);

			firstSqr = [firstSqr;
						S C;
						C E;
						W N];
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
	r = repmat(A, 1, nrepeat);
	toAdd = 0:inc:(nrepeat-1)*inc;
	toAdd = repelem(toAdd, size(A,1), size(A,2));
	r = r+toAdd;
end

function qpts = calc_qpts(hx, hy)
X1 = [0.023931132287081   0.106170269119576   0.188409405952072
   0.066554067839164   0.295266567779633   0.523979067720101
   0.102717654809626   0.455706020243648   0.808694385677670
	] *hx;

Y1 = [0.188409405952072   0.106170269119576   0.023931132287081
   0.523979067720101   0.295266567779633   0.066554067839165
   0.808694385677670   0.455706020243648   0.102717654809626
   ]*hy;

Wx = [0.069826979901454
   0.229241106359586
   0.200931913738959] * hx*hy;

Wy = [0.277777777777777
   0.444444444444444
   0.277777777777777];

X2 = hx - X1;

Y2 = - Y1;

qpts = {X1 Y1 X2 Y2 Wx Wy};
end
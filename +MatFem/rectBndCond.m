function bc = rectBndCond(spc, cond, range)
% SQRBNDCOND construct the boundary condition for a [square]
% mesh
%
%	spc: the [square] function space to define the boundary
%
%	cond: a 3x3 matrix that defines the boundary conditions,
%	each value corresponds to:
%
%		[left-bottom corner, bottom edge, right-bottom corner;
%		 left edge,			 none sense,  right edge;
%		 left-top corner,	 top edge,	  right-top corner]
%
%	and the values can be
%	1 or 'd': Dirichlet boundary
%	2 or 'n': Neumann boundary
%	3 or 'r': Robin boundary
%	or any mark, see applyNeuRob below
%
%	when cond is a scalar, cond is expand to ones(3)*cond

% only for rectangular space
if ischar(cond)
	cond = 1 * (cond=='d' | cond=='D') + ...
		2 * (cond=='n' | cond=='N') + ...
		3 * (cond=='r' | cond=='R');
end

if isequal(size(cond), [1 1])
	cond = ones(3)*cond;
end

if ~isequal(size(cond), [3 3])
	help RectFcnSpc.bndCond
	error('boundary condition size should be 3x3');
end

bi = getBndIdx(spc, range);
bn = genBndNodes(bi, cond);
be = genBndEdges(bi, cond, spc);
bc = MatFem.BndCond(bn, be, 2, spc);
end

function bi = getBndIdx(spc, range)
P = spc.Pb;
x = P(:,1); y = P(:,2);

Rl = range(1); Rr = range(2); Rb = range(3); Rt = range(4);

l = x==Rl; r = x==Rr; b = y==Rb; t = y==Rt;

lt = l & t; lb = l & b; rt = r & t; rb = r & b;

l = l - lt - lb; r = r - rt - rb; t = t - lt - rt; b = b - lb - rb;

bi = {lb; l; lt; b; t; rb; r; rt};
bi = cellfun(@find, bi, 'UniformOutput', 0);
end

function bn = genBndNodes(bi, cond)
cond = cond(:);
cond(5) = [];
condcell = num2cell(cond);

bn = cellfun(@(idx, cond) [ones(size(idx))*cond idx], bi, condcell, ...
	'UniformOutput', 0);
bn = cell2mat(bn);
end

function be = genBndEdges(bi, cond, spc)

Tb = spc.Tb;
switch spc.basisType
	case 1
		sides = [Tb(:,1) Tb(:,2);
			Tb(:,2) Tb(:,3);
			Tb(:,1) Tb(:,3)];
	case 2
		sides = [Tb(:,1) Tb(:,2) Tb(:,4);
			Tb(:,2) Tb(:,3) Tb(:,5);
			Tb(:,1) Tb(:,3) Tb(:,6)];
	otherwise
		error('function space with basis type %d not supported yet', ...
			  spc.basisType);
end
[lb, l, lt, b, t, rb, r, rt] = bi{:};
l = [lt;l;lb]; t = [lt;t;rt]; r = [rt;r;rb]; b = [lb;b;rb];

be = cellfun(@(x)sides(all(ismember(sides,x),2),:), {l;r;b;t}, 'UniformOutput', 0);
edgecond = num2cell([cond(2,1);cond(2,3);cond(1,2);cond(3,2)]);
be = cellfun(@(x,y) [ones(size(x,1),1)*y x], be, edgecond, 'UniformOutput', 0);
be = cell2mat(be);
end
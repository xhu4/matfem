function [V, bc] = setting1(h, basisType, cond, quadOrder)
if nargin < 4 || isempty(quadOrder)
	quadOrder = 9;
	if nargin < 3 || isempty(cond)
		cond = 'd';
	end
end

% mesh
domain = [-1 1 -1 1];
nx = (domain(2)-domain(1))/h+1;
ny = (domain(4)-domain(3))/h+1;
n = [nx ny];
m = MatFem.rectMesh(domain, n);

% function space
V = MatFem.mesh2spc(m, basisType, quadOrder);
bc = MatFem.rectBndCond(V, cond, domain);
end
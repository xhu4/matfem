tic;
h = [1/8 1/16 1/32 1/64];
basisType = 1;
T = MatFem.calcErrors();
for ih = h
	T = [T;evalElas(ih, basisType, 1)];
end
disp(T);
toc;

function errorTable = evalElas(h, basisType, plotresult)
if nargin < 3
	plotresult = 0;
end
import MatFem.*
quadOrder = 9;

% mesh
domain = [0 2 0 1];
nx = (domain(2)-domain(1))/h+1;
ny = (domain(4)-domain(3))/h+1;
n = [nx ny];
m = rectMesh(domain, n);

% function space & boundary condition
V = MatFem.mesh2spc(m, basisType, quadOrder);
bc = MatFem.rectBndCond(V, 'd', domain);

% problem parameters & functions
u = @(x, y) exp(x+y);
g = @(x, y) u(x,y)+onbnd(x,y);
f = @(x,y) -10*exp(x+y);
c = [1 2; 3 4];

% assemble the system
Axx = assemble(V, [1 0; 1 0], 1);
Ayy = assemble(V, [0 1; 0 1], 1);
Axy = assemble(V, [1 0; 0 1], 1);
Ayx = assemble(V, [0 1; 1 0], 1);
A = c(1,1)*Axx+c(1,2)*Ayx+c(2,1)*Axy+c(2,2)*Ayy;

b = assemble(V, [0 0], f);

[b, A] = bc.applyDir(1, g, b, A);

if plotresult
	figure(round(1/h))
	x = A\b;
	subplot(1,2,1);
	V.plotu(x);
	title('numeric');
	
	subplot(1,2,2);
	V.plotu(u);
	title('true');
end

errorTable = calcErrors(V, x, u, {u, u}, h);
end


function r = onbnd(x,y)
x = x(:); y = y(:);
assert(all(x==0|x==2|y==0|y==1));
r=0;
end
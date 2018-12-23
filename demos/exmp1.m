tic;
h = [1/8 1/16 1/32 1/64];
basisType = 2;
T = MatFem.calcErrors();
for ih = h
	T = [T;evalExmp1(ih, basisType, 0)];
end
disp(T);
toc;

function errorTable = evalExmp1(h, basisType, plotresult)
if nargin < 3
	plotresult = 0;
end
import MatFem.*
quadOrder = 9;

% mesh
domain = [-1 1 -1 1];
nx = (domain(2)-domain(1))/h+1;
ny = (domain(4)-domain(3))/h+1;
n = [nx ny];
m = rectMesh(domain, n);

% function space
V = MatFem.mesh2spc(m, basisType, quadOrder);

% problem
f = @(x,y) -y.*(1-y).*(1-x-x.^2/2).*exp(x+y) ...
		   -x.*(1-x/2).*(-3*y-y.^2).*exp(x+y);
u = @(x,y) x.*y.*(1-x/2).*(1-y).*exp(x+y);
g = @(x,y) u(x,y)*onbnd(x,y);

dudx = @(x,y) .5*(x.^2-2).*(y-1).*y.*exp(x+y);
dudy = @(x,y) .5*x.*(x-2).*(y.^2+y-1).*exp(x+y);

% boundary condition
bc = MatFem.rectBndCond(V, 'd', domain);

% assemble
A = assemble(V, [0 1; 0 1], 1);
A = A + assemble(V, [1 0; 1 0], 1);
b = assemble(V, [0 0], f);
[b, A] = bc.applyDir('d', g, b, A);

% solution
sol = A\b;

if plotresult
figure(1/h);
subplot(1,2,1);
V.plotu(sol, 'EdgeColor', [.7, .7, .7]);
title('numeric');

U = u(V.Pb(:,1), V.Pb(:,2));
subplot(1,2,2);
V.plotu(U, 'EdgeColor', [.7, .7, .7]);
title('true');
end

errorTable = calcErrors(V, sol, u, {dudx,dudy}, h);

end

function r = onbnd(x, y)
assert(all(x==1|x==-1|y==1|y==-1))
r = 1;
end
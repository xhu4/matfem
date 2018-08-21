tic;
hs = [1/8 1/16 1/32 1/64 1/128];
T = MatFem.calcErrors();
for h = hs
	T = [T; evalExmp2(h, 1, 0)];
end
disp(T);
toc;

function T = evalExmp2(h, basisType, plotresult)
if nargin < 3
	plotresult = 0;
end
import MatFem.*
quadOrder = 9;

domain = [-1 1 -1 1];
lx = domain(2)-domain(1);
ly = domain(4)-domain(3);
nx = round(lx/h)+1;
ny = round(ly/h)+1;

m = rectMesh(domain, [nx ny]);
V = mesh2spc(m, basisType, quadOrder);
boundMarker = [1 2 1
			   1 1 1
			   1 1 1];
bc =  rectBndCond(V, boundMarker, domain);

A = assemble(V, [1 0; 1 0], 1);
A = A + assemble(V, [0 1; 0 1], 1);
b = assemble(V, [0 0], @f);

b = b + bc.applyEdge(2, [0 0], @bn);
[b, A] = bc.applyDir(1, @bd, b, A);

sol = A\b;

T = calcErrors(V, sol, @u, {@u, @u}, h);

if plotresult
figure(1/h);
subplot(1,2,1);
V.plotu(sol, 'EdgeColor', 'interp');
title('numeric');

U = u(V.Pb(:,1), V.Pb(:,2));
subplot(1,2,2);
V.plotu(U, 'EdgeColor', 'interp');
title('true');
end

end

function r = u(x,y)
r = exp(x+y);
end

function r = f(x,y)
r = -2*exp(x+y);
end

function r = bd(x,y)
assert(all(x==1|x==-1|y==1));
r = exp(x+y);
end

function r = bn(x,y)
assert(all(y(:)==-1));
r = -exp(x+y);
end
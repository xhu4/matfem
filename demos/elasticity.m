tic;
h = [1/8 1/16 1/32 1/64];
basisType = 2;
T = MatFem.calcErrors();
for ih = h
	T = [T;evalElas(ih, basisType, 0)];
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
domain = [-1 1 -1 1];
nx = (domain(2)-domain(1))/h+1;
ny = (domain(4)-domain(3))/h+1;
n = [nx ny];
m = rectMesh(domain, n);

% function space
V = MatFem.mesh2spc(m, basisType, quadOrder);

% problem
% parameters
lambda = 2.0;
mu = 3.0;

% solutions:
u1 = @(x, y) exp(x+y);
u2 = @(x, y) x.*y;

% f & g
f1 = @(x,y) -((lambda+3*mu)*exp(x+y)+(lambda+mu));
f2 = @(x,y) -(lambda+mu)*exp(x+y);
g1 = @(x,y) exp(x+y);
g2 = @(x,y) x.*y;

% boundary condition
bc = MatFem.rectBndCond(V, 'd', domain);

% assemble
A1 = assemble(V, [1,0;1,0], @(x, y)lambda);
A2 = A1/lambda*mu;
A3 = assemble(V, [0,1;0,1], @(x, y)mu);

A4 = assemble(V, [1,0;0,1], @(x, y)lambda);
A5 = A4'/lambda*mu;

A6 = A5/mu*lambda;
A7 = A4/lambda*mu;

A8 = A3/mu*lambda;
A9 = A3;
A10 = A2;
% Assemble Load Vectors

b1 = assemble(V, [0 0], f1);
b2 = assemble(V, [0 0], f2);
% Build Linear System

Aup = [A1+2*A2+A3, A4+A5];
Adn = [A8+2*A9+A10, A6+A7];
% Apply Boundary Conditions

[b1, Aup, ~] = bc.applyDir(1, g1, b1, Aup);
[b2, Adn, ~] = bc.applyDir(1, g2, b2, Adn);
Adn = swap(Adn);		% TODO: Not efficient for sparse matrices
% Get Solution

A = [Aup;Adn];
b = [b1;b2];
x = A\b;

xl = numel(x);
x1 = x(1:xl/2);
x2 = x(xl/2+1:end);

if plotresult
figure(1/h);
subplot(2,2,1);
V.plotu(x1, 'EdgeColor', 'interp');
title('numeric u1');

subplot(2,2,2);
V.plotu(x2, 'EdgeColor', 'interp');
title('numeric u2');

U1 = u1(V.Pb(:,1), V.Pb(:,2));
subplot(2,2,3);
V.plotu(U1, 'EdgeColor', 'interp');
title('true u1');

U2 = u2(V.Pb(:,1), V.Pb(:,2));
subplot(2,2,4);
V.plotu(U2, 'EdgeColor', 'interp');
title('true u2');
end

errorTable = calcErrors(V, x1, u1, {u1, u1}, h);

end

function r = onbnd(x, y)
assert(all(x==1|x==-1|y==1|y==-1))
r = 1;
end

function A = swap(A)
    n = size(A, 2);
    assert(mod(n, 2)==0);
    h = n/2;
    A = [A(:, h+1:end), A(:,1:h)];
end

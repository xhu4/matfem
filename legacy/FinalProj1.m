%% Final Project Part I.
clear; close all
hx = 1/8;
hy = hx;
basisType = 2;

left = 0;
right = 2;
bottom = 0;
top = 1;

nx = (right-left)/hx*basisType+1;
ny = (top - bottom)/hy*basisType+1;

[Pb, Tb] = Gen2dTriMesh([left, right, bottom, top], hx, hy, basisType);
bc = GenBoundaryNodes(nx, ny);
be = GenBoundaryEdges([left, right, bottom, top], hx, hy, Tb);

Nb = size(Pb,2);
A = zeros(Nb);
b = zeros(Nb,1);
c = [1 2; 3 4];

A_nqpts = 5;
b_nqpts = 5;
v_nqpts = 3;



A = A + AssembleStiff2d(Pb, Tb, A_nqpts, [1,0;1,0], @(X,Y)c(1,1), basisType);
A = A + AssembleStiff2d(Pb, Tb, A_nqpts, [1,0;0,1], @(X,Y)c(1,2), basisType);
A = A + AssembleStiff2d(Pb, Tb, A_nqpts, [0,1;1,0], @(X,Y)c(2,1), basisType);
A = A + AssembleStiff2d(Pb, Tb, A_nqpts, [0,1;0,1], @(X,Y)c(2,2), basisType);

b = AssembleLoad2d(Pb, Tb, b_nqpts, [0 0], @(x,y)f(x,y), basisType);
[A, b] = ApplyBoundaryDirchlet2d(A, b, bc, Pb, @(x)g(x));
u = A\b;


U = reshape(u, ny, nx);

x = left:hx/basisType:right;
y = bottom:hy/basisType:top;

[X,Y] = meshgrid(x, y);
U_truth = solution(X, Y);

Ux = @(X,Y)solution(X,Y);
Uy = @(X,Y)solution(X,Y);

err_L2 = errorNorm(u, @(X,Y)solution(X,Y), Pb, Tb, basisType, 5, 'L2');
err_inf = errorNorm(u, @(X,Y)solution(X,Y), Pb, Tb, basisType, 5, 'inf');
err_H1 = errorNorm(u, {Ux, Uy}, Pb, Tb, basisType, 5, 'H1');

fprintf('h = %s, inf error, L2 error, H1 error\n', strtrim(rats(hx)))
fprintf('          %6.2g,    %6.2g,    %6.2g\n', err_inf, err_L2, err_H1)

figure(1);
subplot(1,2,1);
surf(X,Y,U);
title('Numerical solution');
xlabel('X');
ylabel('Y');
subplot(1,2,2);
surf(X,Y,U_truth);
title('True solution');


function U = solution(X, Y)
U = exp(X+Y);
end


function R = EvalC(~, ~)
R = 1;
end


function R = f(X,Y)

R = -10*exp(X+Y);

end


function r = g(X)
x = X(1,:);
y = X(2,:);

if ~all(x==0 | x==2 | y==0 | y==1)
	error('Point (%.2f, %.2f) is not on the boundary', x, y);
end

r = exp(x+y);

end


function r = q(X,Y)
assert(all(Y==-1));
r = zeros(size(X));
end


function re = ComputeR(X,Y)
assert(all(Y==-1));
re = ones(size(X));
end
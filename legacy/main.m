%% Final Project Part II.
clear; close all
hx = 1/8;
hy = hx;
basis_type = 2;

theta = 1;
dt = 8*hx^3;
T = 1;

left = 0;
right = 2;
bottom = 0;
top = 1;


nx = (right-left)/hx*basis_type+1;
ny = (top - bottom)/hy*basis_type+1;
nt = T/dt + 1;

[Pb, Tb] = Gen2dTriMesh(left, right, bottom, top, hx, hy, basis_type);
bc = GenBoundaryNodes(nx, ny);
be = GenBoundaryEdges(left, right, bottom, top, hx, hy, Tb);

Nb = size(Pb,2);
c = 2;

A_nqpts = 5;
b_nqpts = 5;
v_nqpts = 3;

x = left:hx/basis_type:right;
y = bottom:hy/basis_type:top;
[X,Y] = meshgrid(x,y);

Sol = zeros(ny*nx, nt);
Sol(:,1) = reshape(exp(X+Y),[],1);

M = AssembleStiff2d(Pb, Tb, A_nqpts, [0,0;0,0], @(X,Y)1, basis_type);

A = AssembleStiff2d(Pb, Tb, A_nqpts, [1,0;1,0], @(X,Y)c, basis_type) ...
	+ AssembleStiff2d(Pb, Tb, A_nqpts, [0,1;0,1], @(X,Y)c, basis_type);

b = AssembleLoad2d(Pb, Tb, b_nqpts, [0 0], @(x,y)f(x,y,0), basis_type);

[A, b, M] = ApplyBoundaryDirchlet2d(A, b, bc, Pb, @(x)g(x,0), M);

for i = 1:nt-1
	t = (i-1)*dt;
	tn = t+dt;
	
	bn = AssembleLoad2d(Pb, Tb, b_nqpts, [0 0], @(x,y)f(x,y,tn), basis_type);
	[~, bn, ~] = ApplyBoundaryDirchlet2d([], bn, bc, Pb, @(x)g(x,tn));
	
	A_t = M/dt+theta*A;
	b_t = theta*bn+(1-theta)*b+(M/dt-(1-theta)*A)*Sol(:,i);
	Sol(:,i+1) = A_t\b_t;
	
	b = bn;
end

U = reshape(Sol, ny, nx, nt);

U_truth = solution(X, Y, 0:dt:T);

error = norm(U(:)-U_truth(:),2);
display(error);

figure(1);

for i = 1:nt
subplot(1,2,1);
surf(X,Y,U(:,:,i));
title('Numerical solution');
xlabel('X');
ylabel('Y');
subplot(1,2,2);
surf(X,Y,U_truth(:,:,i));
title('True solution');
pause(dt);
drawnow;
end


function U = solution(X, Y, T)
U = zeros(size(X,1),size(X,2),numel(T));
for i = 1:numel(T)
	U(:,:,i) = exp(X+Y+T(i));
end
end


function R = EvalC(~, ~)
R = 1;
end


function R = f(X,Y,t)

R = -3*exp(X+Y+t);

end


function r = g(X,t)
x = X(1,:);
y = X(2,:);
r = exp(x+y+t);
end


function r = q(X,Y)
assert(all(Y==-1));
r = zeros(size(X));
end


function re = ComputeR(X,Y)
assert(all(Y==-1));
re = ones(size(X));
end
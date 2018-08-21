function T = stokesSteady

harray = [1/4 1/8 1/16];
T = MatFem.calcErrors();
basisType = 2;
for h = harray
	T = [T; evalStokes(h, basisType, 1)];
end
end

function errorTable = evalStokes(h, basisType, plotresult)
assert(basisType > 1);
import MatFem.*
[V, bc] = setting1(h, basisType);
[Vp, bcp] = setting1(h, basisType-1);

% define the problem
nu = 1;
p = @(x,y) x.*y;
u1 = @(x,y) 4*x.^2.*y;
u2 = @(x,y) -4*x.*y.^2;

g1 = @(x,y) onbnd(x,y)*u1(x,y);
g2 = @(x,y) onbnd(x,y)*u2(x,y);

px = @(x,y) y;
py = @(x,y) x;

u1x = @(x,y) 8*x.*y;
u1y = @(x,y) 4*x.^2;
u2x = @(x,y) -4*y.^2;
u2y = @(x,y) -8*x.*y;

u1xx = @(x,y) 8*y;
u1yy = @(x,y) 0;
u1xy = @(x,y) 8*x;
u2xx = @(x,y) 0;
u2yy = @(x,y) -8*x;
u2xy = @(x,y) -8*y;


f1 = @(x,y) - 2*nu*u1xx(x,y) + px(x,y) - nu*(u1yy(x,y)+u2xy(x,y));
f2 = @(x,y) - 2*nu*u2yy(x,y) + py(x,y) - nu*(u1xy(x,y)+u2xx(x,y));

if plotresult
	figure(1);
	hdl = zeros(2,3);
	for i = 1:6
		hdl(i) = subplot(3,2,i);
	end
	hdl = hdl';
	
end

% assemble stiffness
Auv = assemble(V, [0 0; 0 0], 1);
Axx = assemble(V, [1 0; 1 0], 1);
Ayy = assemble(V, [0 1; 0 1], 1);
Axy = assemble(V, [1 0; 0 1], 1);
Ayx = assemble(V, [0 1; 1 0], 1);

A11 = 2*nu*Axx + nu*Ayy;
A12 = nu*Ayx;
A21 = nu*Axy;
A22 = nu*Axx + 2*nu*Ayy;
A13 = -assemble([V, Vp], [1 0; 0 0], 1);
A23 = -assemble([V, Vp], [0 1; 0 0], 1);
A31 = -assemble([Vp, V], [0 0; 1 0], 1);
A32 = -assemble([Vp, V], [0 0; 0 1], 1);
O = sparse(Vp.nb,Vp.nb);

A = [A11 A12 A13; 
	 A21 A22 A23; 
	 A31 A32 O  ];

[~, A] = bc.applyDir(1, @(x,y)g1(x,y), [], A);
[~, A] = bc.applyDir(1, @(x,y)g2(x,y), [], A, [], V.nb);

b3 = zeros(Vp.nb, 1);
bndidx = bcp.bn(:,2);
cols = 1:size(A,2);
cols(bndidx(1)+2*V.nb) = [];
A(2*V.nb+bndidx(1), cols) = 0;
A(2*V.nb+bndidx(1), 2*V.nb+bndidx(1)) = 1;
b3(bndidx(1)) = p(-1,-1);

b1 = assemble(V, [0 0], f1);
b2 = assemble(V, [0 0], f2);

b1 = bc.applyDir(1, g1, b1);
b2 = bc.applyDir(1, g2, b2);
b = [b1; b2; b3];	
x = A\b;
x1 = x(1:V.nb);
x2 = x((V.nb+1):(2*V.nb));
pv = x((2*V.nb+1):end);

if plotresult
	plotus(hdl, x1, x2, pv);
end


errorTable = calcErrors(V, x1, u1, {u1x, u1y}, h);

	function plotus(hdl, x1, x2, pv)
		subplot(hdl(1,1));
		V.plotu(x1);
		subplot(hdl(2,1));
		V.plotu(x2);
		subplot(hdl(3,1));
		Vp.plotu(pv);
		subplot(hdl(1,2));
		V.plotu(@(x,y)u1(x,y));
		subplot(hdl(2,2));
		V.plotu(@(x,y)u2(x,y));
		subplot(hdl(3,2));
		Vp.plotu(@(x,y)p(x,y));
		drawnow
	end
end

function r = onbnd(x, y)
assert(all(x(:)==1|y(:)==1|x(:)==-1|y(:)==-1));
r=1;
end
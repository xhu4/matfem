function T = stokes

harray = [1/4 1/8 1/16];
T1 = MatFem.calcErrors();
T2 = T1; T3 = T1;
T = {T1, T2, T3};
basisType = 2;
for h = harray
	Tn = evalStokes(h, basisType, 0);
	T = cellfun(@(x,y) [x;y], T, Tn, 'UniformOutpu', false);
end
celldisp(T);
end

function T = evalStokes(h, basisType, plotresult)
assert(basisType > 1);
import MatFem.*
[V, bc] = setting1(h, basisType);
[Vp, bcp] = setting1(h, basisType-1);

dt = h^2;
t0 = 0;
tn = 1;

% define the problem
nu = 1;
[pp, ff] = dps1;
p = ff.p;
u1 = ff.uc1;
u2 = ff.uc2;

g1 = @(x,y,t) onbnd(x,y)*u1(x,y,t);
g2 = @(x,y,t) onbnd(x,y)*u2(x,y,t);

px = ff.dpdx;
py = ff.dpdy;

u1x = ff.duc1dx;
u1y = ff.duc1dy;
u1t = ff.duc1dt;
u2x = ff.duc2dx;
u2y = ff.duc2dy;
u2t = ff.duc2dt;

u1xx = ff.duc1dxx;
u1yy = ff.duc1dyy;
u1xy = ff.duc1dxy; 
u2xx = ff.duc1dxx;
u2yy = ff.duc2dyy;
u2xy = ff.duc2dxy;


f1 = ff.f1;
f2 = ff.f2;

x10 = V.project(@(x,y) u1(x, y, t0));
x20 = V.project(@(x,y) u2(x, y, t0));

if plotresult
	figure(1);
	hdl = zeros(2,3);
	for i = 1:6
		hdl(i) = subplot(3,2,i);
	end
	hdl = hdl';
	
	plotus(hdl, x10, x20, Vp.project(p), t0);
end

% assemble stiffness
Auv = assemble(V, [0 0; 0 0], 1);
Axx = assemble(V, [1 0; 1 0], 1);
Ayy = assemble(V, [0 1; 0 1], 1);
Axy = assemble(V, [1 0; 0 1], 1);
Ayx = assemble(V, [0 1; 1 0], 1);

A11 = Auv/dt + 2*nu*Axx + nu*Ayy;
A12 = nu*Axy;
A21 = nu*Ayx;
A22 = Auv/dt + nu*Axx + 2*nu*Ayy;
A13 = -assemble([V, Vp], [1 0; 0 0], 1);
A23 = -assemble([V, Vp], [0 1; 0 0], 1);
A31 = assemble([Vp, V], [0 0; 1 0], 1);
A32 = assemble([Vp, V], [0 0; 0 1], 1);
O = sparse(Vp.nb, Vp.nb);

A = [A11 A12 A13; 
	 A21 A22 A23; 
	 A31 A32 O  ];

[~, A] = bc.applyDir(1, @(x,y)g1(x,y,t0), [], A);
[~, A] = bc.applyDir(1, @(x,y)g2(x,y,t0), [], A, [], V.nb);

bndidx = bcp.bn(:,2);
cols = 1:size(A,2);
cols(bndidx(1)+2*V.nb) = [];
A(2*V.nb+bndidx(1), cols) = 0;
A(2*V.nb+bndidx(1), 2*V.nb+bndidx(1)) = 1;

t = t0:dt:tn;
x1old = x10;
x2old = x20;
for now = t(2:end)
	b = load(now, x1old, x2old);
	x = A\b;
	x1 = x(1:V.nb);
	x2 = x((V.nb+1):(2*V.nb));
	pv = x((2*V.nb+1):end);
	if plotresult
		plotus(hdl, x1, x2, pv, now);
	end
	x1old = x1;
	x2old = x2;
end

T1 = calcErrors(V, x1, @(x,y)u1(x,y,tn), ...
	{@(x,y)u1x(x,y,tn), @(x,y)u1y(x,y,tn)}, h);
T2 = calcErrors(V, x2, @(x,y)u2(x,y,tn), ...
	{@(x,y)u2x(x,y,tn), @(x,y)u2y(x,y,tn)}, h);
T3 = calcErrors(Vp, pv, @(x,y)p(x,y,tn), ...
	{@(x,y)px(x,y,tn), @(x,y)py(x,y,tn)}, h);
T = {T1, T2, T3};

	function b = load(t, u1old, u2old)
		import MatFem.*
		b1 = assemble(V, [0 0], @(x,y)f1(x,y,t)) ...
			+ Auv*u1old/dt;
		
		b2 = assemble(V, [0 0], @(x,y)f2(x,y,t)) ...
			+ Auv*u2old/dt;
		
		b3 = zeros(Vp.nb, 1);
		b3(1) = p(-1,-1,t);
		
		b1 = bc.applyDir(1, @(x,y)g1(x,y,t), b1);
		b2 = bc.applyDir(1, @(x,y)g2(x,y,t), b2);
		b = [b1; b2; b3];
	end

	function plotus(hdl, x1, x2, pv, t)
		subplot(hdl(1,1));
		V.plotu(x1);
		subplot(hdl(2,1));
		V.plotu(x2);
		subplot(hdl(3,1));
		Vp.plotu(pv);
		subplot(hdl(1,2));
		V.plotu(@(x,y)u1(x,y,t));
		subplot(hdl(2,2));
		V.plotu(@(x,y)u2(x,y,t));
		subplot(hdl(3,2));
		Vp.plotu(@(x,y)p(x,y));
		drawnow
	end
end

function r = onbnd(x, y)
assert(all(x(:)==1|y(:)==1|x(:)==-1|y(:)==-1));
r=1;
end
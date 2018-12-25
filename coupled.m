function Tables = coupled
clear; close all;
set(0, 'DefaultTextFontSize', 16)
set(0, 'DefaultAxesFontSize', 16)
harray = [1/16];
T_pm = MatFem.calcErrors();
T_pf = T_pm; T_uc = T_pm; T_p=T_pm;
basisType = 2;
% video = VideoWriter('coupled.mp4', 'MPEG-4');
% open(video)
gif = 'coupled.gif';
for h = harray
	tic;
	Ts = evalCoupled(h, basisType, 1, gif);
	[pm,pf,uc,p] = Ts{:};
	T_pm = [T_pm;pm];
	T_pf = [T_pf;pf];
	T_uc = [T_uc;uc];
	T_p = [T_p;p];
	toc;
end
Tables = {T_pm, T_pf, T_uc, T_p};
celldisp(Tables);
% close(video)
close
end

function errorTables = evalCoupled(h, basisType, plotresult, gif)
if nargin < 4
	gif = [];
if nargin < 3
	plotresult = false;
end
end

figure(1)
set(gcf, 'Position', [0 0 1280 500])

quadOrder = 9;

dom.D = [0 1 0 0.75];
dom.C = [0 1 -0.25 0];

spc = genspc1(dom, h, basisType, quadOrder);
bc = bndcond1(spc, dom);

t0 = 0; tt = 1; dt = h^2;
t = t0:dt:tt;

% define problem
[p, f] = dps1;
p.dt = dt;

t_onbnd_C = @(x,y) all(x(:)==0|x(:)==1|y(:)==-.25&y<=0);
t_onbnd_D = @(x,y) all(x(:)==0|x(:)==1|y(:)==.75&y>=0);

g.pm = @(x,y,t) f.pm(x,y,t)/t_onbnd_D(x,y);
g.pf = @(x,y,t) f.pf(x,y,t)/t_onbnd_D(x,y);
g.uc1 = @(x,y,t) f.uc1(x,y,t)/t_onbnd_C(x,y);
g.uc2 = @(x,y,t) f.uc2(x,y,t)/t_onbnd_C(x,y);

% assemble A's
[A, Auv] = assembleA1(spc, bc, p);

% boundary A
A = applybndA(bc, A, spc);

old = {spc.D.project(@(x,y)f.pm(x,y,t0));
	spc.D.project(@(x,y)f.pf(x,y,t0));
	spc.C.project(@(x,y)f.uc1(x,y,t0));
	spc.C.project(@(x,y)f.uc2(x,y,t0));
	0};

fprintf('%d DoFs\n', size(A,1));

if plotresult
	plotus(old, t0, spc, f, gif);
end

for now = t(2:end)
	% assemble b's
	b = assembleb1(spc, bc, p, f, g, Auv, old, now);
	
	x = A\b;
	x = mat2cell(x, [spc.D.nb,spc.D.nb,spc.C.nb,spc.C.nb,spc.Cm1.nb], 1);
	checkbnd(x, now, bc, f);
	
	xt = cellfun(@(x,y) x.project(@(a,b)y(a,b,now)), ...
		{spc.D; spc.D; spc.C; spc.C; spc.Cm1}, ...
		{f.pm; f.pf; f.uc1; f.uc2; f.p}, 'UniformOutput', false);
	xt = cell2mat(xt);
	error = A*xt-b;
	if plotresult
		plotus(x, now, spc, f, gif);
	end
	
	old = x;
end

import MatFem.calcErrors

e_pm = calcErrors(spc.D, x{1}, @(x,y)f.pm(x,y,tt), ...
	{@(x,y)f.dpmdx(x,y,tt), @(x,y)f.dpmdy(x,y,tt)}, h);
e_pf = calcErrors(spc.D, x{2}, @(x,y)f.pf(x,y,tt), ...
	{@(x,y)f.dpfdx(x,y,tt), @(x,y)f.dpfdy(x,y,tt)}, h);
e_uc = calcErrors(spc.C, x(3:4), {@(x,y)f.uc1(x,y,tt), @(x,y)f.uc2(x,y,tt)}, ...
	{@(x,y)f.duc1dx(x,y,tt), @(x,y)f.duc1dy(x,y,tt);
	@(x,y)f.duc2dx(x,y,tt), @(x,y)f.duc2dy(x,y,tt)}, h);
e_p = calcErrors(spc.Cm1, x{5}, @(x,y)f.p(x,y,tt), ...
	{@(x,y)f.dpdx(x,y,tt), @(x,y)f.dpdy(x,y,tt)}, h);
errorTables = {e_pm, e_pf, e_uc, e_p};

end

function [b1to2, b2to1] = idxmaps(Pb1, Pb2, ycond)
	b1 = find(Pb1(:,2)==ycond);
	b2 = find(Pb2(:,2)==ycond);
	Pbb1x = Pb1(b1,1);
	Pbb2x = Pb2(b2,1);
	[s1, i1] = sort(Pbb1x);
	[s2, i2] = sort(Pbb2x);
	assert(isequal(s1,s2));
	b1 = b1(i1);
	b2 = b2(i2);
	b1to2(b1) = b2;
	b2to1(b2) = b1;
end

function A = mapidx(A, map, dim, siz)
[i, j, v] = find(A);
switch dim
	case 1
		i = map(i);
		assert(~any(i==0));
	case 2
		j = map(j);
		assert(~any(j==0));
	otherwise
		error('dim should be 1 or 2.');
end
A = sparse(i,j,v,siz(1),siz(2));
end



function spc = genspc1(dom, h, basisType, quadOrder)
import MatFem.*

msh_D = rectMesh(dom.D, h);
msh_C = rectMesh(dom.C, h);

spc.D = mesh2spc(msh_D, basisType, quadOrder);
spc.C = mesh2spc(msh_C, basisType, quadOrder);
spc.Cm1 = mesh2spc(msh_C, basisType-1, quadOrder);
end

function bc = bndcond1(spc, dom)
cond_D = [1 2 1; 1 1 1; 1 1 1];
cond_C = [1 1 1; 1 1 1; 1 2 1];
bc.D = MatFem.rectBndCond(spc.D, cond_D, dom.D);
bc.C = MatFem.rectBndCond(spc.C, cond_C, dom.C);
end

function [A, Auv] = assembleA1(spc, bc, p)
import MatFem.assemble
Auv_D = assemble(spc.D, [0 0; 0 0]);
Axx_D = assemble(spc.D, [1 0; 1 0]);
Ayy_D = assemble(spc.D, [0 1; 0 1]);
Auv_C = assemble(spc.C, [0 0; 0 0]);
Axx_C = assemble(spc.C, [1 0; 1 0]);
Ayy_C = assemble(spc.C, [0 1; 0 1]);
Axy_C = assemble(spc.C, [1 0; 0 1]);
Ayx_C = assemble(spc.C, [0 1; 1 0]);

A11 = (p.phi_m*p.C_mt/p.dt+p.sigma*p.k_m/p.mu)*Auv_D+p.k_m/p.mu*(Axx_D+Ayy_D);
A12 = -p.sigma*p.k_m/p.mu*Auv_D;
A21 = A12;
A22 = (p.phi_f*p.C_ft/p.dt+p.sigma*p.k_m/p.mu)*Auv_D+p.k_f/p.mu*(Axx_D+Ayy_D);

A33 = p.eta*(Auv_C/p.dt+2*p.nu*(Axx_C+1/2*Ayy_C));
A44 = p.eta*(Auv_C/p.dt+2*p.nu*(Ayy_C+1/2*Axx_C));
A34 = p.eta*p.nu*Ayx_C;
A43 = p.eta*p.nu*Axy_C;
A35 = -p.eta*assemble([spc.C, spc.Cm1], [1 0; 0 0]);
A45 = -p.eta*assemble([spc.C, spc.Cm1], [0 1; 0 0]);
A53 = p.eta*assemble([spc.Cm1, spc.C], [0 0; 1 0]);
A54 = p.eta*assemble([spc.Cm1, spc.C], [0 0; 0 1]);

Auv_CD_C = bc.C.applyEdge(2, [0 0; 0 0]);
coef = p.alpha*p.nu*sqrt(p.dim)/sqrt(trace(p.PI));
A24 = -Auv_CD_C;
A42 = p.eta/p.rho*Auv_CD_C;
A33 = A33 + p.eta*coef*Auv_CD_C;
A32 = p.eta*coef*p.k_f/p.mu*bc.C.applyEdge(2, [0 0; 1 0]);

[idxc2d, ~] = idxmaps(spc.C.Pb, spc.D.Pb, 0);
A24 = mapidx(A24, idxc2d, 1, [spc.D.nb, spc.C.nb]);
A42 = mapidx(A42, idxc2d, 2, [spc.C.nb, spc.D.nb]);
A32 = mapidx(A32, idxc2d, 2, [spc.C.nb, spc.D.nb]);

Ocd = sparse(spc.C.nb, spc.D.nb);
Odc = Ocd';
Oqd = sparse(spc.Cm1.nb, spc.D.nb);
Odq = Oqd';
Oqq = sparse(spc.Cm1.nb, spc.Cm1.nb);

A = [A11 A12 Odc Odc Odq
	A21 A22 Odc A24 Odq
	Ocd A32 A33 A34 A35
	Ocd A42 A43 A44 A45
	Oqd Oqd A53 A54 Oqq];

Auv.D = Auv_D;
Auv.C = Auv_C;
end

function A = applybndA(bc, A, spc)
[~,A] = bc.D.applyDir(1,1,[],A);
[~,A] = bc.D.applyDir(1,1,[],A,[], spc.D.nb);
[~,A] = bc.C.applyDir(1,1,[],A,[], spc.D.nb*2);
[~,A] = bc.C.applyDir(1,1,[],A,[], spc.D.nb*2+spc.C.nb);
A(spc.D.nb*2+spc.C.nb*2+1,:) = 0;
A(spc.D.nb*2+spc.C.nb*2+1,spc.D.nb*2+spc.C.nb*2+1) = 1;
end

function b = assembleb1(spc, bc, p, f, g, Auv, old, time)
import MatFem.assemble
[pm_old, pf_old, uc1_old, uc2_old, ~] = old{:};
b1 = p.phi_m*p.C_mt/p.dt*(Auv.D*pm_old) + assemble(spc.D, [0 0], @(x,y)f.Qplus(x,y,time));
b2 = p.phi_f*p.C_ft/p.dt*(Auv.D*pf_old) + assemble(spc.D, [0 0], @(x,y)f.qp(x,y,time));
b3 = p.eta*(assemble(spc.C, [0 0], @(x,y)f.f1(x,y,time)) + Auv.C*uc1_old/p.dt);
b4 = p.eta*(assemble(spc.C, [0 0], @(x,y)f.f2(x,y,time)) + Auv.C*uc2_old/p.dt);

b1 = bc.D.applyDir(1,@(x,y)g.pm(x,y,time),b1);
b2 = bc.D.applyDir(1,@(x,y)g.pf(x,y,time),b2);
b3 = bc.C.applyDir(1,@(x,y)g.uc1(x,y,time),b3);
b4 = bc.C.applyDir(1,@(x,y)g.uc2(x,y,time),b4);

pt1 = spc.Cm1.Pb(1,:);
o = zeros(spc.Cm1.nb, 1);
o(1) = f.p(pt1(1), pt1(2), time);
b = [b1; b2; b3; b4; o];
end

function checkbnd(sol, time, bc, f)
[pm, pf, uc1, uc2, ~] = sol{:};
[iD, XD, YD] = bc.D.nodes(1);
[iC, XC, YC] = bc.C.nodes(1);
pmb = pm(iD);
pfb = pf(iD);
uc1b = uc1(iC);
uc2b = uc2(iC);
pmbt = f.pm(XD,YD,time);
pfbt = f.pf(XD,YD,time);
uc1bt = f.uc1(XC, YC, time);
uc2bt = f.uc2(XC, YC, time);
assert(isequal(pmb, pmbt));
assert(isequal(pfb, pfbt));
assert(isequal(uc1b, uc1bt));
assert(isequal(uc2b, uc2bt));
end


function plotus(sol, time, spc, f, video)
if nargin<5 || isempty(video)
	writevideo = false;
else
	writevideo = true;
end

figure(1)

qscale = 0.1;

[pm, pf, uc1, uc2, ~] = sol{:};
figure(1);
subplot(2,3,1);
spc.D.plotu(pm);
title('\fontsize{24}numerical p_m')
axis([0 1 0 0.75])
caxis([-0.4 0.2])
% axis equal;

subplot(2,3,4);
spc.D.plotu(@(x,y)f.pm(x,y,time));
title('\fontsize{24}true p_m')
axis([0 1 0 0.75])
caxis([-0.4 0.2])
% axis equal;

subplot(2,3,2);
spc.D.plotu(pf);
title('\fontsize{24}numerical p_f')
axis([0 1 0 0.75])
caxis([-2 2])
% axis equal;

subplot(2,3,5);
spc.D.plotu(@(x,y)f.pf(x,y,time));
title('\fontsize{24}true p_f')
axis([0 1 0 0.75])
caxis([-2 2])
% axis equal;

subplot(2,3,3);
spc.C.quiver(uc1*qscale, uc2*qscale);
axis([-0.3 1.3 -0.5 0.25])
title('\fontsize{24}numerical u_c')
% axis equal;

subplot(2,3,6);
spc.C.quiver(@(x,y)f.uc1(x,y,time)*qscale,@(x,y)f.uc2(x,y,time)*qscale);
axis([-0.3 1.3 -0.5 0.25])
title('\fontsize{24}true u_c')
% axis equal
% currentFigure = gcf;
% title(currentFigure.Children(end), ['t=',num2str(time)]);
% mtit(['time = ', num2str(time, '%10.5f')])
TimeBox = uicontrol('style', 'text');
set(TimeBox, 'String', ['time = ', num2str(time, '%10.5f')])
set(TimeBox, 'Position', [100,0,200,30])
set(TimeBox, 'FontSize', 24)
drawnow;

if writevideo
	frame = getframe(gcf);
	im = frame2im(frame);
	[imind, cm] = rgb2ind(im, 256);
	if time == 0
		imwrite(imind, cm, video, 'gif', 'Loopcount', inf);
	else
		imwrite(imind, cm, video, 'gif', 'WriteMode', 'append');
	end
% 	writeVideo(video, frame);
end
% pause
end

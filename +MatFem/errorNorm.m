function ret = errorNorm(fcnSpc, u, u_true, normType)
% ERRORNORM calculates the norm difference between u (vector) and u_true
% (function)
% 
%	See also INTOVERMESH


switch normType
	case 'inf'
		ret = errorNormInf(u, u_true, fcnSpc);
	
	case 'L2'
		ret = errorNormL2(u, u_true, fcnSpc);
		
	case 'H1'
		ret = errorNormH1(u, u_true, fcnSpc);
		
	otherwise
		error('FEM:InvalidArgument', 'unknown norm type %s', normType);
end

end


function ret = errorNormInf(u, u_true, fcnSpc)
ret = 0;
u = makecell(u);
u_true = makecell(u_true);
for i = 1:length(u)
	uu = u{i};
	uut = u_true{i};
	[X, Y] = fcnSpc.quad{:};
	ut = uut(X,Y);
	uv = fcnSpc.evalOnQuad(uu);
	ret = ret + (ut(:)-uv(:)).^2;
end
ret = sqrt(ret);
ret = norm(ret, 'inf');
end


function ret = errorNormL2(u, u_true, fcnSpc)
ret = 0;
u = makecell(u);
u_true = makecell(u_true);
for i = 1:length(u)
	uu = u{i};
	uut = u_true{i};
	[X, Y] = fcnSpc.quad{:};
	ut = uut(X,Y);
	uv = fcnSpc.evalOnQuad(uu);
	ret = ret + (ut-uv).^2;
end
ret = ret * fcnSpc.wt .*fcnSpc.area;
ret = sqrt(sum(ret));
end


function ret = errorNormH1(u, u_true, fcnSpc)
u = makecell(u);
u_true = makecell(u_true);
ret = 0;
for i = 1:length(u)
	uu = u{i};
	uut = u_true(i,:);
	dutdx = uut{1};
	dutdy = uut{2};
	[X,Y] = fcnSpc.quad{:};
	utx = dutdx(X,Y);
	uty = dutdy(X,Y);
	uvx = fcnSpc.evalOnQuad(uu, [1 0]);
	uvy = fcnSpc.evalOnQuad(uu, [0 1]);
	ret = ret + (utx-uvx).^2+(uty-uvy).^2;
end
ret = ret*fcnSpc.wt.*fcnSpc.area;
ret = sqrt(sum(ret));
end

function u = makecell(u)
if ~iscell(u)
	u = {u};
end
end
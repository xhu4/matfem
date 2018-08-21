function ret = errorNorm(u, u_true, Pb, Tb, basisType, nqpts, normType)
% ERRORNORM calculates the norm difference between u (vector) and u_true
% (function)
% 
%	See also INTOVERMESH

switch normType
	case 'inf'
		ret = errorNormInf(u, u_true, Pb, Tb, basisType, nqpts);
	
	case 'L2'
		ret = errorNormL2(u, u_true, Pb, Tb, basisType, nqpts);
		
	case 'H1'
		ret = errorNormH1(u, u_true, Pb, Tb, basisType, nqpts);
		
	otherwise
		error('FEM:InvalidArgument', 'unknown norm type %s', normType);
end

end


function ret = errorNormInf(u, u_true, Pb, Tb, basisType, nqpts)
Ne = size(Tb, 2);			% # of elements
ret = 0;
for iE = 1:Ne
	vE = VtcsOfElem(Pb, Tb, iE);
	[X, Y, ~, ~] = triquad(nqpts, vE);
	ut = u_true(X, Y);
	uh = EvalFcn(u, X, Y, Pb, Tb, iE, basisType, [0 0]);
	diff = reshape(ut-uh, [], 1);
	ret = max(ret, max(abs(diff)));
end
end


function ret = errorNormL2(u, u_true, Pb, Tb, basisType, nqpts)
EvalU = @(X, Y, iE) EvalFcn(u, X, Y, Pb, Tb, iE, basisType, [0 0]);
f = @(X, Y, iE) (u_true(X,Y)-EvalU(X, Y, iE)).^2;
ret = IntOverMesh(Pb, Tb, f, nqpts, 'iE');
ret = sqrt(ret);
end


function ret = errorNormH1(u, u_true, Pb, Tb, basisType, nqpts)
dudx = @(X, Y, iE) EvalFcn(u, X, Y, Pb, Tb, iE, basisType, [1 0]);
dudy = @(X, Y, iE) EvalFcn(u, X, Y, Pb, Tb, iE, basisType, [0 1]);
dutdx = u_true{1};
dutdy = u_true{2};
f = @(X, Y, iE) (dutdx(X,Y)-dudx(X,Y,iE)).^2+(dutdy(X,Y)-dudy(X,Y,iE)).^2;
ret = IntOverMesh(Pb, Tb, f, nqpts, 'iE');
ret = sqrt(ret);
end
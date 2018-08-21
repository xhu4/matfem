function ret = IntOverMesh(Pb, Tb, f, nqpts, ext)
% INTOVERMESH integrates function f over the mesh (Pb, Tb)
%
% See also INTOVERELEM, ERRORNORM

Ne = size(Tb, 2);			% # of elements

ret = 0;

for iE = 1:Ne
	vE = VtcsOfElem(Pb, Tb, iE);
	
	if nargin == 4
		g = f;
	elseif ext == 'iE'
		g = @(X, Y)f(X, Y, iE);
	end
	
	ret = ret + IntOverElem(vE, g, nqpts);
end
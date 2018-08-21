function A = genStiff(fcnSpc, ders, coefFcn, nqpts)
% genStiff assembles the stiffness matrix.
%	GENSTIFF builds the matrix A with
%	A(i,j) = \int phi_i? * phi_j? dxdy
%	where phi_i? is the derivative of phi_i determined by `derivatives`.
%
% See also GENLOAD

if nargin < 4
	nqpts = 5;
	if nargin < 3
		coefFcn = @(x,y)1;
	end
end

Nb = fcnSpc.nBases;				% # of basis functions,
Nlb = fcnSpc.nLocalBases;		% # of local bases and
Ne = fcnSpac.nElems;			% # of elements

A = zeros(Nb,Nb);


for i = 1:Ne
	[vE, iE] = fcnSpc.elem(i);
	[X, Y, Wx, Wy] = triquad(nqpts, vE(1:3));
	
	for b1 = 1:Nlb
		for b2 = 1:Nlb
			r = IntegrateGauss2d(vE, [b1, b2], ders, ...
								 coefFcn, X, Y, Wx, Wy, basis_type);
			A(iE(b2), iE(b1)) = A(iE(b2), iE(b1)) + r;
		end
	end
end
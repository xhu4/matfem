function A = AssembleStiff2d(Pb, Tb, nqpts, der_orders, c, basis_type)
% ASSEMBLESTIFF2D assembles the stiffness matrix.
%	ASSEMBLESTIFF2D builds the matrix A with
%	A(i,j) = \int phi_i? * phi_j? dxdy
%	where phi_i? is the derivative of phi_i determined by der_orders.
%
% See also ASSEMBLELOAD2D, INTEGRATEGAUSS2D

Nb = size(Pb, 2);				% # of basis functions
[Nlb, Ne] = size(Tb);			% # of local bases and # of elements

assert(Nlb == 3*basis_type);


% Initialize Matrix A
A = zeros(Nb,Nb);


% For each element
for iE = 1:Ne
	vE = Pb(:,Tb(1:3,iE))';
	[X, Y, Wx, Wy] = triquad(nqpts, vE);
	
	for b1 = 1:Nlb
		for b2 = 1:Nlb
			r = IntegrateGauss2d(vE, [b1, b2], der_orders, ...
								 c, X, Y, Wx, Wy, basis_type);
			A(Tb(b2, iE), Tb(b1, iE)) = A(Tb(b2, iE), Tb(b1, iE)) + r;
		end
	end
end
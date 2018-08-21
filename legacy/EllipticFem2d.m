function U = EllipticFem2d(Pb, Tb, bc, be, c, f, g, p, q, r, A_der_orders, ...
	b_der_orders, v_der_orders, R_der_orders, w_der_orders, A_nqpts, b_nqpts, v_nqpts, basis_type)

% [Pb, Tb] = Gen2dTriMesh(left, right, bottom, top, h_x, h_y, basis_type);

Nb = size(Pb,2);
A = zeros(Nb);
b = zeros(Nb,1);

for i = 1:size(A_der_orders,3)
	A = A + AssembleStiff2d(Pb, Tb, A_nqpts, A_der_orders(:,:,i), c, basis_type);
end

for i = 1:size(b_der_orders,1)
	b = AssembleLoad2d(Pb, Tb, b_nqpts, b_der_orders(i,:), f, basis_type);
end


% if has Robin condition
if ~isempty(q) || ~isempty(r)
	[R, w] = ApplyBoundaryRobin(be, Pb, Tb, v_nqpts, R_der_orders, ...
		w_der_orders, @(x,y)c(x,y).*r(x,y), @(x,y)c(x,y).*q(x,y), basis_type);
	A = A+R;
	b = b+w;
end


% if has Neumann condition
if ~isempty(p)
	p_tilde = @(x,y)c(x,y).*p(x,y);
	v = ApplyBoundaryNeumann(be, Pb, Tb, v_nqpts, v_der_orders, p_tilde, basis_type);
	b = b+v;
end

[A, b] = ApplyBoundaryDirchlet2d(A, b, bc, Pb, g);

U = A \ b;

end




function [A, b, M] = ApplyBoundaryDirchlet2d(A, b, bc, Pb, EvalG, M)

bn = bc(2,(bc(1,:)==1));

if ~isempty(b)
	b(bn) = EvalG(Pb(:,bn));
end

if ~isempty(A)
	A(bn,:)=0;
	for i = bn
		A(i,i)=1;
	end
else
	A=[];
end

if nargin==6 && ~isempty(M)
	M(bn,:)=0;
else
	M=[];
end

function r = assemble(spc, ders, f)
if nargin < 3
	f = 1;
end

if length(spc) == 1
	if isa(f, 'function_handle')
		F = f(spc.quad{:});
	elseif isnumeric(f)
		F = f;
	end
	
	wt = spc.wt(:);
	area = spc.area(:);
	
	n2int = size(ders, 1);
	switch n2int
		case 1
			Phi = spc.evalAll(ders);
			for i = 1:spc.nlb
				Phi{i} = ((Phi{i}.*F) * wt).*area;
			end
			r = accumarray(spc.Tb(:), cell2mat(Phi), [spc.nb, 1]);
		case 2
			i = zeros(spc.ne,spc.nlb,spc.nlb);
			j = i;
			val = i;
			
			Phi1 = spc.evalAll(ders(1,:));
			if isequal(ders(1,:), ders(2,:))
				Phi2 = Phi1;
			else
				Phi2 = spc.evalAll(ders(2,:));
			end
			for ii = 1:spc.nlb
				for jj = 1:spc.nlb
					i(:,ii,jj) = spc.Tb(:,ii);
					j(:,ii,jj) = spc.Tb(:,jj);
					val(:,ii,jj) = ((Phi1{ii}.*Phi2{jj}.*F) * wt) .* area;
				end
			end
			r = sparse(i(:), j(:), val(:), spc.nb, spc.nb);
	end
else
	assert(length(spc)==2 && size(ders,1)==2);
	spc1 = spc(1);
	spc2 = spc(2);
	assert(spc1.ne == spc2.ne);
	
	if isa(f, 'function_handle')
		F = f(spc.quad{:});
	elseif isnumeric(f)
		F = f;
	end
	
	i = zeros(spc1.ne,spc1.nlb,spc2.nlb);
	j = i;
	val = i;
	Phi1 = spc1.evalAll(ders(1,:));
	Phi2 = spc2.evalAll(ders(2,:));
	for ii = 1:spc1.nlb
		for jj = 1:spc2.nlb
			i(:,ii,jj) = spc1.Tb(:,ii);
			j(:,ii,jj) = spc2.Tb(:,jj);
			val(:,ii,jj) = ((Phi1{ii}.*Phi2{jj}.*F) * spc1.wt) .* spc1.area;
		end
	end
	r = sparse(i(:), j(:), val(:), spc1.nb, spc2.nb);
end
end


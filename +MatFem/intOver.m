% function r = intOver(v, basisType, ders, coefFcn, qpts)
% nv = size(v,1);
% nlb = nLocals(nv, basisType);
% 
% [X,Y,Wx,Wy] = qpts{:};
% 
% % evaluate c
% if isnumeric(coefFcn)
% 	c = coefFcn*ones(size(X));
% else
% 	c = coefFcn(X,Y);
% end
% 
% % can normalize vE: vE(:,1) = vE(:,1)-vE(1,1), ...
% 
% ndim = size(ders,1);
% phi = zeros(size(X,1), size(X,2), ndim, nlb);
% for b = 1:ndim
% 	for i = 1:nlb
% 		% evaluate phi_i / phi_j
% 		phi(:,:,b,i) = MatFem.evalTriElem(X,Y,v,i,basisType,ders(b,:));
% 	end
% end
% 
% switch ndim
% 	case 1
% 		r = zeros(nlb, 1);
% 		for i = 1:nlb
% 			r(i) = Wx'*(c.*phi(:,:,1,i))*Wy;
% 		end
% 	case 2
% 		r = zeros(nlb);
% 		for i = 1:nlb
% 			for j = 1:nlb
% 				r(i,j) = Wx'*(c.*phi(:,:,1,i).*phi(:,:,2,j))*Wy;
% 			end
% 		end
% 	otherwise
% 		error('cannot integrate more than 2 bases');
% end
% 
% end
% 
% function nlb = nLocals(nv, basisType)
% lookup = [2 3; 3 6];
% nlb = lookup(nv-1,basisType);
% end
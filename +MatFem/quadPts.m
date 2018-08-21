function[P, w, area] = quadPts(order, vertices, ndim)
% QUADPTS compute the quadrature points in a simplex
%
% vertices: each row of vertices is a input, e.g.
%	vertices = [v1 v2 v3 ...]
%	describes two 2D triangle elements, the first one contains three
%	vertices with coordinates (x11, y11), (x12, y12), (x13, y13).
%
% X: ?

nv = size(vertices, 2)/ndim;
ne = size(vertices, 1);
switch nv
	case 2
		[lambda, w] = quadpts1(order);
		area = vertices(:,1:ndim)-vertices(:,ndim+1:end);
		area = sqrt(sum(area.^2, 2));
	case 3
		[lambda, w] = quadpts2(order);
		switch ndim
			case 2
				area = polyarea(vertices(:,1:ndim:end), ...
								vertices(:,2:ndim:end), 2);
			case 3
				error('area in 3d not implemented yet');
				% TODO: implement 3D area calc
			otherwise
				assert(0, 'cannot have a %dD polygon', ndim)
		end
	case 4
		[lambda, w] = quadpts3(order);
		assert(ndim==3,'cannot have a tetrahedran in %dD',ndim);
		a = vertices(:,1:3)-vertices(:,10:12);
		b = vertices(:,4:6)-vertices(:,10:12);
		c = vertices(:,7:9)-vertices(:,10:12);
		area = dot(a, cross(b,c,2), 2)/6;
	otherwise
		error('cannot compute quadrature with %d vertices', nv);
end
v = reshape(vertices, [], nv);
X = v*lambda';
P = mat2cell(X, ne*ones(1,ndim));
w = w(:);
area = area(:);
end


function [lambda,weight] = quadpts1(order)
% QUADPTS1 quadrature points in 1-D.
%
% [lambda,weight] = QUADPTS1(order) return quadrature points with given
% order (up to 19) in the barycentric coordinates.
%
% The output lambda is a matrix of size nQ by 2, where nQ is the number of
% quadrature points. lambda(i,:) is two barycentric coordinate of the
% i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
% of all quadrature points. The coordinate of the p-th quadrature point
% can be computed as 
%
%    p = lambda(p,1)*a+ lambda(p,2)*b;
% where a and b are the end points of segment in 2D or 3D;
%
% References: 
% * Pavel Holoborodko
% http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/
% 
% See also quadpts, quadpts3, verifyquadpts
%
% Added by Huayi Wei
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

numPts = ceil((order+1)/2);

if numPts > 10
   numPts = 10; 
end

switch numPts
    case 1
        A = [0      2.0000000000000000000000000];
        
    case 2
        A = [0.5773502691896257645091488 	1.0000000000000000000000000
            -0.5773502691896257645091488 	1.0000000000000000000000000];
        
    case 3
        A = [0 	0.8888888888888888888888889
            0.7745966692414833770358531 	0.5555555555555555555555556
            -0.7745966692414833770358531 	0.5555555555555555555555556];
        
    case 4
        A = [0.3399810435848562648026658 	0.6521451548625461426269361
            0.8611363115940525752239465 	0.3478548451374538573730639
            -0.3399810435848562648026658 	0.6521451548625461426269361
            -0.8611363115940525752239465 	0.3478548451374538573730639];
        
    case 5
        A = [0 	                            0.5688888888888888888888889
            0.5384693101056830910363144 	0.4786286704993664680412915
            0.9061798459386639927976269 	0.2369268850561890875142640
            -0.5384693101056830910363144 	0.4786286704993664680412915
            -0.9061798459386639927976269 	0.2369268850561890875142640];
        
    case 6
        A = [0.2386191860831969086305017 	0.4679139345726910473898703
            0.6612093864662645136613996 	0.3607615730481386075698335
            0.9324695142031520278123016 	0.1713244923791703450402961
            -0.2386191860831969086305017 	0.4679139345726910473898703
            -0.6612093864662645136613996 	0.3607615730481386075698335
            -0.9324695142031520278123016 	0.1713244923791703450402961];
        
    case 7
        A = [0 	                            0.4179591836734693877551020
            0.4058451513773971669066064 	0.3818300505051189449503698
            0.7415311855993944398638648 	0.2797053914892766679014678
            0.9491079123427585245261897 	0.1294849661688696932706114
            -0.4058451513773971669066064 	0.3818300505051189449503698
            -0.7415311855993944398638648 	0.2797053914892766679014678
            -0.9491079123427585245261897 	0.1294849661688696932706114];
        
    case 8
        A = [0.1834346424956498049394761 	0.3626837833783619829651504
            0.5255324099163289858177390 	0.3137066458778872873379622
            0.7966664774136267395915539 	0.2223810344533744705443560
            0.9602898564975362316835609 	0.1012285362903762591525314
            -0.1834346424956498049394761 	0.3626837833783619829651504
            -0.5255324099163289858177390 	0.3137066458778872873379622
            -0.7966664774136267395915539 	0.2223810344533744705443560
            -0.9602898564975362316835609 	0.1012285362903762591525314];
        
    case 9
        A = [0 	                            0.3302393550012597631645251
            0.3242534234038089290385380 	0.3123470770400028400686304
            0.6133714327005903973087020 	0.2606106964029354623187429
            0.8360311073266357942994298 	0.1806481606948574040584720
            0.9681602395076260898355762 	0.0812743883615744119718922
            -0.3242534234038089290385380 	0.3123470770400028400686304
            -0.6133714327005903973087020 	0.2606106964029354623187429
            -0.8360311073266357942994298 	0.1806481606948574040584720
            -0.9681602395076260898355762 	0.0812743883615744119718922];
        
    case 10
        A = [0.1488743389816312108848260 	0.2955242247147528701738930
            0.4333953941292471907992659 	0.2692667193099963550912269
            0.6794095682990244062343274 	0.2190863625159820439955349
            0.8650633666889845107320967 	0.1494513491505805931457763
            0.9739065285171717200779640 	0.0666713443086881375935688
            -0.1488743389816312108848260 	0.2955242247147528701738930
            -0.4333953941292471907992659 	0.2692667193099963550912269
            -0.6794095682990244062343274 	0.2190863625159820439955349
            -0.8650633666889845107320967 	0.1494513491505805931457763
            -0.9739065285171717200779640 	0.0666713443086881375935688];
end
lambda1 = (A(:,1)+1)/2;
lambda2 = 1 - lambda1;
lambda = [lambda1, lambda2];
weight = A(:,2)/2;
end

function [lambda,weight] = quadpts2(order)
% QUADPTS quadrature points in 2-D.
%
% [lambda,weight] = quadpts2(order) return quadrature points with given
% order (up to 9) in the barycentric coordinates.
%
% The output lambda is a matrix of size nQ by 3, where nQ is the number of
% quadrature points. lambda(i,:) is three barycentric coordinate of the
% i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
% of all quadrature points. The x-y coordinate of the p-th quadrature point
% can be computed as
%
%     pxy = lambda(p,1)*node(elem(:,1),:) ...
%         + lambda(p,2)*node(elem(:,2),:) ...
%         + lambda(p,3)*node(elem(:,3),:);
%
% The weight of p-th quadrature point is given by weight(p). See
% verifyquadpts for the usage of qudrature rules to compute integrals over
% triangles.
%
% References:
%
% * David Dunavant. High degree efficient symmetrical Gaussian
%    quadrature rules for the triangle. International journal for numerical
%    methods in engineering. 21(6):1129--1148, 1985.
% * John Burkardt. DUNAVANT Quadrature Rules for the Triangle.
%    http://people.sc.fsu.edu/~burkardt/m_src/dunavant/dunavant.html
%
% See also quadpts1, quadpts3, verifyquadpts
%
% Order 6 - 9 is added by Huayi Wei, modify by Jie Zhou
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if order>9
	order = 9;
end
switch order
	case 1     % Order 1, nQuad 1
		lambda = [1/3, 1/3, 1/3];
		weight = 1;
	case 2     % Order 2, nQuad 3
		lambda = [2/3, 1/6, 1/6; ...
			1/6, 2/3, 1/6; ...
			1/6, 1/6, 2/3];
		weight = [1/3, 1/3, 1/3];
	case 3     % Order 3, nQuad 4
		lambda = [1/3, 1/3, 1/3; ...
			0.6, 0.2, 0.2; ...
			0.2, 0.6, 0.2; ...
			0.2, 0.2, 0.6];
		weight = [-27/48, 25/48, 25/48, 25/48];
	case 4     % Order 4, nQuad 6
		lambda = [0.108103018168070, 0.445948490915965, 0.445948490915965; ...
			0.445948490915965, 0.108103018168070, 0.445948490915965; ...
			0.445948490915965, 0.445948490915965, 0.108103018168070; ...
			0.816847572980459, 0.091576213509771, 0.091576213509771; ...
			0.091576213509771, 0.816847572980459, 0.091576213509771; ...
			0.091576213509771, 0.091576213509771, 0.816847572980459];
		weight = [0.223381589678011, 0.223381589678011, 0.223381589678011, ...
			0.109951743655322, 0.109951743655322, 0.109951743655322];
	case 5     % Order 5, nQuad 7
		alpha1 = 0.059715871789770;      beta1 = 0.470142064105115;
		alpha2 = 0.797426985353087;      beta2 = 0.101286507323456;
		lambda = [   1/3,    1/3,    1/3; ...
			alpha1,  beta1,  beta1; ...
			beta1, alpha1,  beta1; ...
			beta1,  beta1, alpha1; ...
			alpha2,  beta2,  beta2; ...
			beta2, alpha2,  beta2; ...
			beta2,  beta2, alpha2];
		weight = [0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506, ...
			0.125939180544827, 0.125939180544827, 0.125939180544827];
	case 6
		A =[0.249286745170910  0.249286745170910  0.116786275726379
			0.249286745170910  0.501426509658179  0.116786275726379
			0.501426509658179  0.249286745170910  0.116786275726379
			0.063089014491502  0.063089014491502  0.050844906370207
			0.063089014491502  0.873821971016996  0.050844906370207
			0.873821971016996  0.063089014491502  0.050844906370207
			0.310352451033784  0.636502499121399  0.082851075618374
			0.636502499121399  0.053145049844817  0.082851075618374
			0.053145049844817  0.310352451033784  0.082851075618374
			0.636502499121399  0.310352451033784  0.082851075618374
			0.310352451033784  0.053145049844817  0.082851075618374
			0.053145049844817  0.636502499121399  0.082851075618374];
		lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
		weight = A(:,3);
	case 7
		A =[0.333333333333333  0.333333333333333 -0.149570044467682
			0.260345966079040  0.260345966079040  0.175615257433208
			0.260345966079040  0.479308067841920  0.175615257433208
			0.479308067841920  0.260345966079040  0.175615257433208
			0.065130102902216  0.065130102902216  0.053347235608838
			0.065130102902216  0.869739794195568  0.053347235608838
			0.869739794195568  0.065130102902216  0.053347235608838
			0.312865496004874  0.638444188569810  0.077113760890257
			0.638444188569810  0.048690315425316  0.077113760890257
			0.048690315425316  0.312865496004874  0.077113760890257
			0.638444188569810  0.312865496004874  0.077113760890257
			0.312865496004874  0.048690315425316  0.077113760890257
			0.048690315425316  0.638444188569810  0.077113760890257];
		lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
		weight = A(:,3);
	case 8
		A =[0.333333333333333  0.333333333333333  0.144315607677787
			0.081414823414554  0.459292588292723  0.095091634267285
			0.459292588292723  0.081414823414554  0.095091634267285
			0.459292588292723  0.459292588292723  0.095091634267285
			0.658861384496480  0.170569307751760  0.103217370534718
			0.170569307751760  0.658861384496480  0.103217370534718
			0.170569307751760  0.170569307751760  0.103217370534718
			0.898905543365938  0.050547228317031  0.032458497623198
			0.050547228317031  0.898905543365938  0.032458497623198
			0.050547228317031  0.050547228317031  0.032458497623198
			0.008394777409958  0.263112829634638  0.027230314174435
			0.008394777409958  0.728492392955404  0.027230314174435
			0.263112829634638  0.008394777409958  0.027230314174435
			0.728492392955404  0.008394777409958  0.027230314174435
			0.263112829634638  0.728492392955404  0.027230314174435
			0.728492392955404  0.263112829634638  0.027230314174435];
		lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
		weight = A(:,3);
		
	case 9
		A =[0.333333333333333  0.333333333333333  0.097135796282799
			0.020634961602525  0.489682519198738  0.031334700227139
			0.489682519198738  0.020634961602525  0.031334700227139
			0.489682519198738  0.489682519198738  0.031334700227139
			0.125820817014127  0.437089591492937  0.07782754100474
			0.437089591492937  0.125820817014127  0.07782754100474
			0.437089591492937  0.437089591492937  0.07782754100474
			0.623592928761935  0.188203535619033  0.079647738927210
			0.188203535619033  0.623592928761935  0.079647738927210
			0.188203535619033  0.188203535619033  0.079647738927210
			0.910540973211095  0.044729513394453  0.025577675658698
			0.044729513394453  0.910540973211095  0.025577675658698
			0.044729513394453  0.044729513394453  0.025577675658698
			0.036838412054736  0.221962989160766  0.043283539377289
			0.036838412054736  0.741198598784498  0.043283539377289
			0.221962989160766  0.036838412054736  0.043283539377289
			0.741198598784498  0.036838412054736  0.043283539377289
			0.221962989160766  0.741198598784498  0.043283539377289
			0.741198598784498  0.221962989160766  0.043283539377289];
		lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
		weight = A(:,3);
end
end

function [lambda,weight] = quadpts3(order)
% QUADPTS3 quadrature points in 3-D.
%
% [lambda,weight] = quadpts(order) return quadrature points with given
% order (up to 5) in the barycentric coordinates.
%
% The output lambda is a matrix of size nQ by 3, where nQ is the number of
% quadrature points. lambda(i,:) is three barycentric coordinate of the
% i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
% of all quadrature points. The x-y coordinate of the p-th quadrature point
% can be computed as 
%
%     pxy = lambda(p,1)*node(elem(:,1),:) ...
%         + lambda(p,2)*node(elem(:,2),:) ... 
%         + lambda(p,3)*node(elem(:,3),:) ... 
%         + lambda(p,4)*node(elem(:,4),:);
%
% The weight of p-th quadrature point is given by weight(p). See
% verifyquadpts for the usage of qudrature rules to compute integrals over
% triangles.
% 
% References: 
%
% * Jinyun, Y. Symmetric Gaussian quadrature formulae for tetrahedronal
% regions. Comput. Methods Appl. Mech. Engrg.. 43(3):349--353, 1984.
%  
% See also quadpts1, quadpts, verifyquadpts
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if order>5
    order = 5;
end
switch order
    case 1    % Order 1, nQuad 1
        lambda = [1/4, 1/4, 1/4, 1/4];
        weight = 1;
    case 2    % Order 2, nQuad 4
        alpha = 0.5854101966249685; 
        beta =  0.138196601125015;
        lambda = [alpha beta beta beta; ....
                  beta alpha beta beta; ...
                  beta beta alpha beta; ...
                  beta beta beta alpha];
        weight = [1/4, 1/4, 1/4, 1/4];
    case 3    % Order 3, nQuad 5
        lambda = [1/4 1/4 1/4 1/4; ...
                  1/2 1/6 1/6 1/6; ...
                  1/6 1/2 1/6 1/6; ...
                  1/6 1/6 1/2 1/6; ...
                  1/6 1/6 1/6 1/2];
        weight = [-4/5, 9/20, 9/20, 9/20, 9/20];
    case 4    % Order 4, nQuad 16
        alpha1 = 0.7716429020672371; 
        beta1 =  0.7611903264425430e-1;
        w1 = 0.5037379410012282e-1;
        alpha = 0.4042339134672644;
        beta = 0.7183164526766925e-1;
        gamma = 0.11970052777978019;
        w2 = 0.6654206863329239e-1;
        lambda = [alpha1 beta1 beta1 beta1; ....
                  beta1 alpha1 beta1 beta1; ...
                  beta1 beta1 alpha1 beta1; ...
                  beta1 beta1 beta1 alpha1; ...
                  alpha alpha beta gamma; ...
                  alpha alpha gamma beta; ...
                  alpha beta alpha gamma; ...
                  alpha beta gamma alpha; ...
                  alpha gamma beta alpha; ...
                  alpha gamma alpha beta; ...
                  beta alpha alpha gamma; ...
                  beta alpha gamma alpha; ...
                  beta gamma alpha alpha; ...
                  gamma alpha alpha beta; ...
                  gamma alpha beta alpha; ...
                  gamma beta alpha alpha];                  
        weight = [w1, w1, w1, w1, ...
                  w2, w2, w2, w2, w2, w2, ...
                  w2, w2, w2, w2, w2, w2];
    case 5    % Order 5, nQuad 17
        alpha1 = 0.7316369079576180; 
        beta1 =  0.8945436401412733e-1;
        w1 = 0.6703858372604275e-1;
        alpha = 0.4214394310662522;
        beta = 0.2454003792903000e-1;
        gamma = 0.1325810999384657;
        w2 = 0.4528559236327399e-1;
        lambda = [1/4, 1/4, 1/4, 1/4; ...
                  alpha1 beta1 beta1 beta1; ....
                  beta1 alpha1 beta1 beta1; ...
                  beta1 beta1 alpha1 beta1; ...
                  beta1 beta1 beta1 alpha1; ...
                  alpha alpha beta gamma; ...
                  alpha alpha gamma beta; ...
                  alpha beta alpha gamma; ...
                  alpha beta gamma alpha; ...
                  alpha gamma beta alpha; ...
                  alpha gamma alpha beta; ...
                  beta alpha alpha gamma; ...
                  beta alpha gamma alpha; ...
                  beta gamma alpha alpha; ...
                  gamma alpha alpha beta; ...
                  gamma alpha beta alpha; ...
                  gamma beta alpha alpha];                  
        weight = [0.1884185567365411, ...
                  w1, w1, w1, w1, ...
                  w2, w2, w2, w2, w2, w2, ...
                  w2, w2, w2, w2, w2, w2];
end
end
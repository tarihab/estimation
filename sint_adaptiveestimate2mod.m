
function [dydt] = sint_adaptiveestimate2mod(t,y,K,glx,gly,npa,ind,xc,yc,sigma,sw)

	t  % print the current time
	dydt = zeros(size(y));
	% constants
	na = K(1);  % no.of agents
	k1 = K(2);  % controller gain
	np = K(3);  % no.of parameters
	sigma = K(4);
	g1 = K(5);
	k2 = K(6);  % consensus gain

	n = 2;  % ambient dimension

	% actual state of the agents
	p = zeros(n,na);
	%lambda = zeros(np,na);
        %ahat = zeros(np,na);

	for i=1:na
		nna(i) = ((npa(i)*npa(i) - npa(i))/2) + npa(i);
		ahat{i} = zeros(npa(i),na);  % estimate of i'th parameter set by all agents - 
		% ahat{i}(:,i) is the true estimate (by the i'th agent); others try to follow this 
		% estimate
		lambda{i} = zeros(npa(i),1);
		bi{i} = zeros(npa(i),1);
		ahatdot{i} = zeros(npa(i),1);
	end

	for i=1:na
		p(:,i) = y(((i-1)*n)+1:i*n);

		s = n*na;
		Lambda{i} = zeros(npa(i),npa(i));
		j1 = sum(nna(1:(i-1)));
		j2 = sum(nna(1:i));
		Lambdavec = y(s+j1+1:s+j2);
		Lambda{i} = vec2sm(Lambdavec,npa(i));

		s = s + sum(nna);
		j1 = sum(npa(1:(i-1)));
		j2 = sum(npa(1:i));
		lambda{i} = y(s+j1+1:s+j2);

		s = s + sum(npa);
		%j1 = sum(npa(1:(i-1)));
		%j2 = sum(npa(1:i));
		for j=1:na
			j1 = (np*(i-1)) + sum(npa(1:(j-1)));
			j2 = (np*(i-1)) + sum(npa(1:j));
			ahat{j}(:,i) = y(s+j1+1:s+j2);
		end
		%ahat{i} = y(s+j1+1:s+j2);
	end

	% sensory density function measurements of the agents
	phim = zeros(1,na);
	for i=1:na
		%phim(i) = phifcn(p(1,i),p(2,i),atrue,np);
		phim(i) = fieldvalue(p(1,i),p(2,i));
	end

%% control law
	%u = zeros(n,na);
	u = zeros(n,na);
	for i=1:na
		gl = [glx(i); gly(i)];
		u(:,i) = -k1*(p(:,i)-gl);
	end

	g2 = 1;
	%bi = zeros(np,na);
	%ahatdot = zeros(np,na);
	for i=1:na
		ahatdot{i} = zeros(npa(i),na);
	end

	for i=1:na
		bi{i} = - g1.*(Lambda{i}*ahat{i}(:,i) - lambda{i});	
		%ahatdot{i} = bi{i};
		%for j=1:na
			%if(i~=j)
			%	bi{i} = bi{i} - k2*(ahat{i}-ahat{i});
			%end
		%end
		P = zeros(npa(i),npa(i));
		for j=1:npa(i)
		 if((ahat{i}(j,i)<0) | ((ahat{i}(j,i)==0) & (bi{i}(j)<0)))
		  P(j,j) = 1;
		 end
		end
		%ahatdot(:,i) = bi(:,i);
		ahatdot{i}(:,i) = g2*(bi{i}-P*bi{i}); % projection

		for j=1:na
			if(j~=i)
				ahatdot{i}(:,j) = -(ahat{i}(:,j) - ahat{i}(:,i));
			end
		end
	end

	for i=1:na
		ahatdotvec{i} = [];
		for j=1:na
			ahatdotvec{i} = [ahatdotvec{i}; ahatdot{j}(:,i)];
		end
	end

	beta1 = 1;
	% derivative updates
	for i=1:na
		dydt(((i-1)*n)+1:i*n) = u(:,i);

		s = n*na;
		Lambdadot = zeros(npa(i),npa(i));
		Kt = Kvector(p(1,i),p(2,i),xc,yc,sigma);
		K = Kt(ind{i});
		Kbar = [];
		abar = [];
		for j=1:na
			if(j~=i)
				Kbar = [Kbar; Kt(ind{j})];
				abar = [abar; ahat{j}(:,i)];
			end
		end
		Kbar = Kbar';

		Lambdadot = sw(i).*(K*K');
		Lambdadotvec = mattovecmod(Lambdadot);
		j1 = sum(nna(1:(i-1)));
		j2 = sum(nna(1:i));
		dydt(s+j1+1:s+j2) = Lambdadotvec;

		s = s + sum(nna);
		lambdadot = sw(i).*(K.*(phim(i)-(Kbar*abar)));
		j1 = sum(npa(1:(i-1)));
		j2 = sum(npa(1:i));
		dydt(s+j1+1:s+j2) = lambdadot;

		s = s + sum(npa);
		%j1 = sum(npa(1:(i-1)));
		%j2 = sum(npa(1:i));
		%dydt(s+j1+1:s+j2) = ahatdot{i};
		dydt(s+((i-1)*np)+1:s+(i*np)) = ahatdotvec{i};
	end

end

function mat = vectomat(vec,l)
% converts vector into matrix by stacking each consecutive l elements as a column
% of the matrix; this requires the length of vec to be divisible by l

	n = length(vec)/l;
	mat = zeros(l,n);
	for(i=1:n)
		mat(:,i) = vec(((i-1)*l)+1:i*l);
	end

end

function vec = mattovec(mat)
% converts matrix into vector by stacking each column below the other

	[m,n] = size(mat);
	vec = [];
	for i=1:n
		vec = [vec; mat(:,i)];
	end

end

function vec = mattovecmod(mat)
% converts matrix into vector by stacking each column below the other

	[m,n] = size(mat);
	vec = [];
	for i=1:m
		vec = [vec; mat(i,i:end)'];
	end

end

% code for integration over an arbitrary polygon
% obtained from http://www.caam.rice.edu/~timwar/CAAM210/2D_Integrals.html
function intf = polygonint(f, xv, yv, a, b, c, tol)
% f is the handle to the function to be integrated
% xv and yv specify the vertices of the polygon
% a and b are extra parameters to the function f
% tol is the tolerance

  global gxv gyv gf ga gb gc
  gxv = xv; gyv = yv; gf = f;
  ga = a; gb = b; gc = c;

  % find bounding box for polygon
  xmin = min(xv); xmax = max(xv);
  ymin = min(yv); ymax = max(yv);
  %disp(xv);
  %disp(yv);

  %disp(xmin);
  %disp(xmax);
  %disp(ymin);
  %disp(ymax);
  % compute approximate area for polygon
  intf = dblquad(@(x,y) integrand(x,y), xmin, xmax, ymin, ymax, tol);
  %intf = dblquad(@(x,y) integrand(x,y,a,b), xmin, xmax, ymin, ymax, tol);
  %intf = integral2(@(x,y) integrand(x,y), xmin, xmax, ymin, ymax,'Method','iterated','AbsTol',0,'RelTol',tol);
  %intf = integral2(@(x,y) integrand(x,y), xmin, xmax, ymin, ymax,'AbsTol',0,'RelTol',tol);

end

%function chif = integrand(x, y, a, b)
function chif = integrand(x, y)

  global gxv gyv gf ga gb gc

  % handle the case when x is a vector and y is a scalar
  % (see help inpolygon)
  if(length(y)==1)
    y = ones(size(x))*y;
  end

  % scatter(x,y, 'r.')
  %disp(x)
  %disp(y)
  %chif = feval(gf, x,y).*inpolygon(x, y, gxv, gyv);
  chif = feval(gf,x,y,ga,gb,gc).*inpolygon(x, y, gxv, gyv);

end

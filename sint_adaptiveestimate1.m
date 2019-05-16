% double integrator dynamics with the adaptive coverage controller

function [dydt] = sint_adaptiveestimate1(t,y,K,glx,gly)

	t  % print the current time
	dydt = zeros(size(y));
	% constants
	na = K(1);  % no.of agents
	k1 = K(2);  % controller gain
	np = K(3);  % no.of parameters
	sigma = K(4);
	gamma = K(5);

	n = 2;  % ambient dimension

	% actual state of the agents
	p = zeros(n,na);
	lambda = zeros(np,na);
        ahat = zeros(np,na);
        nn = ((np*np - np)/2) + np;

	for i=1:na
		p(:,i) = y(((i-1)*n)+1:i*n);

		s = n*na;
		Lambda{i} = zeros(np,np);
		Lambdavec = y(s+((i-1)*nn)+1:s+(i*nn));
		Lambda{i} = vec2sm(Lambdavec,np);

		s = s + nn*na;
		lambda(:,i) = y(s+((i-1)*np)+1:s+(i*np));

		s = s + np*na;
		ahat(:,i) = y(s+((i-1)*np)+1:s+(i*np));
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
	gl = [glx; gly];
	for i=1:na
		u(:,i) = -k1*(p(:,i)-gl);
	end

	bi = zeros(np,na);
	ahatdot = zeros(np,na);
	for i=1:na
	end

	% derivative updates
	for i=1:na
		dydt(((i-1)*n)+1:i*n) = u(:,i);
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

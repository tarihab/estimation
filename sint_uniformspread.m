% double integrator dynamics with the adaptive coverage controller

function [dydt] = sint_uniformspread(t,y,K)

	t  % print the current time
	dydt = zeros(size(y));
	% constants
	na = K(1);
	%k2 = K(2);
	%b = K(3);

	k = K(2); % no. of agents

	% actual state of the agents
	p = zeros(n,na);

	n = 2;  % ambient dimension
	for i=1:na
		p(:,i) = y(((i-1)*n)+1:i*n);
	end

	%disp('Voronoi function called');
	%disp(vx);
	%disp(vy);
        for i=1:n
                for j=1:na
                        if(isnan(p(i,j)) | isinf(p(i,j)))
                                disp('Position NaN or Inf..');
                                disp(i);
                                disp(j);
                        end
                end
                if(max(p(i,:))>2.0 | min(p(i,:))<-1.0)
                %if(max(p(i,:))>1.0 | min(p(i,:))<0)
                        disp('Position bound exceeded..');
                        %disp(p(i,:));
                	b1 = p(i,:) > 2;
                	b2 = p(i,:) < -1;
                	p(i,b1) = 1.99;
                	p(i,b2) = -0.99;
                	%b1 = p(i,:) > 1;
                	%b2 = p(i,:) < 0;
                	%p(i,b1) = 0.999;
                	%p(i,b2) = 0.001;
                end
        end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%% DETERMINE THE PROPER VERTICES W.R.T SQUARE DOMAIN %%%%%%%%
	%%%% VoronoiLimit is an external function obtained from the %%%%%%
	%%%% internet. See the function m-file for more details. %%%%%%%%%
	%%%% The matlab built-in commands 'voronoi' and 'voronoin' %%%%%%%
	%%%% was not sufficient.                      		 %%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%s = [0 0; 0 1; 1 1; 1 0; 0 0];  % bounding surface - unit square
	%s = [-1 -1; -1 2; 2 2; 2 -1];  % bounding surface - unit square
	%[V,C] = VoronoiLimit(p(1,:)',p(2,:)',s);
	%[vornb,vorvx]=polybnd_voronoi(p',s);
	%[Vr,C] = MyVoronoiLimit(p(1,:)',p(2,:)',s);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% compute Lvi, Mvi and Cvi based on current parameter estimates
	Mvi = zeros(1,na);
	Lvi = zeros(n,na);
	Lvi1 = zeros(1,na);
	Lvi2 = zeros(1,na);
	Cvi = zeros(n,na);
	parfor i=1:na
		% get the voronoi partition of i'th agent - vertices
		%vorPar = Vr(C{i},:);
		%vorPar = unique(vorvx{i},'rows');
		% x-coordinates of voronoi partition vertices
		% last element repeated, hence excluding it
		%vx = vorPar(:,1);
		% y-coordinates of voronoi partition vertices
		% last element repeated, hence excluding it
		%vy = vorPar(:,2);
		[vx,vy] = compute_voronoi(i,xborder,yborder,p(1,:)',p(2,:)');
		% Compute Mvi
		Mvi(i) = polygonint(@phifcn,vx,vy,[],[],[],1e-4);
		% Compute Lvi.. component-wise integration..
		%Lvi(1,i) = polygonint(@qphifcn,vx,vy,ahat(:,i),1,1e-6);
		Lvi1(i) = polygonint(@qphifcn,vx,vy,1,[],[],1e-4);
		%Lvi(2,i) = polygonint(@qphifcn,vx,vy,ahat(:,i),2,1e-6);
		Lvi2(i) = polygonint(@qphifcn,vx,vy,2,[],[],1e-4);
		% Compute Cvi
		%Cvi(:,i) = Lvi(:,i)/Mvi(i);
		%Cvi(:,i)
	%end
	%Lvi = [Lvi1; Lvi2];
		Lvi(:,i) = [Lvi1(i); Lvi2(i)];
	%for i=1:na
		Cvi(:,i) = Lvi(:,i)/Mvi(i);
	end

	% control law
	%u = zeros(n,na);
	u = zeros(n,na);
	for i=1:na
		%u(:,i) = B{i}\(-k1*Mvi(i)*S{i}(1:2,:)'*(p(:,i)-Cvi(:,i)) - k2*v(:,i));
		%Bu(:,i) = -k1*Mvi(i)*S{i}(1:2,:)'*(p(:,i)-Cvi(:,i)) - k2*v(:,i);
		u(:,i) = -k1*Mvi(i)*(p(:,i)-Cvi(:,i));
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

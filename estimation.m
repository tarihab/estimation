%% Code for estimation of scalar field

% Domain is unit square
xborder = [0 1 1 0];
yborder = [0 0 1 1];

N = 5; % no.of agents

xcb = 0.05:0.1:1;  
ycb = 0.05:0.1:1;  
xc = [];  % x-coordinates of RBF centres
yc = [];  % y-coordinates of RBF centres
for i=1:length(ycb)
	if(rem(i,2)==0)
		xc = [xc fliplr(xcb)];
	else
		xc = [xc xcb];
	end
	
	yc = [yc ycb(i)*ones(1,length(ycb))];
end
sigma = 0.03;  % std. deviation of RBFs

np = length(xc);  % no.of parameters

k = 1; % controller gain

p0 = rand(2*N,1);  % initial positions at random

y0 = [p0];
tspan = [0 20];
RelTol = 1e-4;
AbsTol = 1e-4;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%{
K = [N; k];

[tout, yout] = ode45(@(t,y) sint_uniformspread(t,y,K,xborder,yborder),tspan,y0,options);

posout = yout(end,:)';
%}
posout = [0.6979; 0.1904; 0.7717; 0.8539; 0.2114; 0.2847; 0.2670; 0.7709; 0.5200; 0.5571];
pos = zeros(2,N);
for i=1:N
	pos(:,i) = posout(((i-1)*2)+1:i*2);
end
disp('Starting positions');
disp(pos);

%% Start estimation

%disp('Partitions');
for i=1:N		
	[vx,vy] = compute_voronoi(i,xborder,yborder,pos(1,:)',pos(2,:)');
	partitionx{i} = vx;
	partitiony{i} = vy;
	%disp(partitionx{i});
	%disp(partitiony{i});

	% locate which of the centres are in partition i
	ind{i} = [];
	for j=length(xc)
		if(inpolygon(xc(j),yc(j),vx,vy))
			ind{i} = [ind{i} j];
		end
	end
end

pause;

%% First algorithm

a0 = ones(N*np,1);
gamma = 1;

nn = ((np*np - np)/2) + np;
Lambda0 = zeros(nn*N,1);
lambda0 = zeros(N*np,1);

y0 = [posout; Lambda0; lambda0; a0];
tspan = [0 1];
k2 = 1;
K = [N; k; np; sigma; gamma; k2];

K = zeros(np,np);
for i=1:np
	for j=1:i
		K(i,j) = Kvalue(xc(i),yc(i),xc(j),yc(j),sigma);
		K(j,i) = K(i,j);
	end
end

flag = 1;
cntr_centre = ones(N,1);
c = 0
cx = zeros(N,1);
cy = zeros(N,1);
for i=1:N
	c = ind{i}(cntr_centre(i));
	cx(i) = xc(c);
	cy(i) = yc(c);
end
e = 0.1;
while(flag~=0)

	[tout, yout] = ode45(@(t,y) sint_adaptiveestimate1(t,y,K,cx,cy,xc,yc,sigma),tspan,y0,options);
	%[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,xborder,yborder,xc,yc),tspan,y0,options);
	y0 = yout(end,:)';
	f = 0;
	for i=1:N
		pi = y0(((i-1)*n)+1:i*n);
		chi = K\(Kvector(pi(1),pi(2),xc,yc,sigma));
		tmp = 0;
		for j=1:np
			if(j~=c)
				tmp = tmp + chi(j);
			end
		end
		if((chi(c)-tmp) > e)
			if(cntr_centre(i)<length(ind{i}))
				cntr_centre(i) = cntr_centre(i) + 1;
				c = ind{i}(cntr_centre(i));
				cx(i) = xc(c);
				cy(i) = yc(c);
			else
				f = f + 1;
			end
		end

		if(f==N)
			flag = 0;
		end
	end

	tspan = [tspan(2) tspan(2)+1];
end

%{
%% Second algorithm

a0 = ones(np,1);
gamma = 1;

for i=1:N
	npa(i) = length(ind{i});
	nna(i) = ((npa(i)*npa(i) - npa(i))/2) + npa(i);
end
nnas = sum(nna);
Lambda0 = zeros(nnas,1);
lambda0 = zeros(np,1);

y0 = [posout; Lambda0; lambda0; a0];
tspan = [0 1];
k2 = 1;
K = [N; k; np; sigma; gamma; k2];
cx = zeros(N,1);
cy = zeros(N,1);

flag = 1;
cntr_centre = 1;
for i=1:N
	c = ind{i}(cntr_centre);
	cx(i) = xc(c);
	cy(i) = yc(c);
end
e = 0.1;
while(flag~=0)

	[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,cx,cy,npa,ind,xc,yc,sigma),tspan,y0,options);
	%[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,xborder,yborder,xc,yc),tspan,y0,options);
	y0 = yout(end,:)';
	f = 0;
	for i=1:N
		pi = y0(((i-1)*n)+1:i*n);
		Kvec = Kvector(pi(1),pi(2),xc,yc,sigma);
		chi = K\(Kvec(ind{i}));
		tmp = 0;
		for j=1:npa(i)
			if((ind{i}(j))~=c)
				tmp = tmp + chi(j);
			end
		end
		if((chi(c)-tmp) > e)
			if(cntr_centre(i)<length(ind{i}))
				cntr_centre(i) = cntr_centre(i) + 1;
				c = ind{i}(cntr_centre(i));
				cx(i) = xc(c);
				cy(i) = yc(c);
			else
				f = f + 1;
			end
		end

		if(f==N)
			flag = 0;
		end
	end

	tspan = [tspan(2) tspan(2)+1];
end
%}

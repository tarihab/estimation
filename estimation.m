%% Code for estimation of scalar field

% Domain is unit square
xborder = [0 1 1 0];
yborder = [0 0 1 1];

N = 4; % no.of agents
n = 2; % ambient dimension = 2

xcb = 0.05:0.1:1;  
ycb = 0.05:0.1:1;  
%xcb = 0.04:0.07:1;  
%ycb = 0.04:0.07:1;  
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
sigmalist = [0.03; 0.04; 0.05];
% sigma = 0.05;  % std. deviation of RBFs
sigma = sigmalist(h);  % std. deviation of RBFs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% for exact parameterization case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xc = [0.2,0.35,0.6,0.85,0.7,0.75,0.15,0.35];
%yc = [0.25,0.26,0.18,0.3,0.75,0.9,0.75,0.6];
%atrue = [2,1,1.5,1.8,1.2,1.6,2.5,1.1];
%sigma = 0.1;
xc = [0.1,0.2,0.35,0.4,0.6,0.7,0.78,0.85,0.7,0.6,0.75,0.85,0.15,0.2,0.32,0.35];
yc = [0.18,0.25,0.26,0.32,0.18,0.25,0.3,0.45,0.75,0.83,0.9,0.60,0.65,0.72,0.75,0.6];
atrue = [2,1,2.6,0.9,1.5,1.6,1.1,1.8,1.2,0.8,2.0,1.6,2.0,1.5,2.5,1.1];
sigma = 0.1;
%xc = [0.1,0.2,0.4,0.6,0.7,0.85,0.7,0.6,0.85,0.15,0.2,0.35];
%yc = [0.18,0.25,0.32,0.18,0.3,0.45,0.75,0.9,0.60,0.65,0.72,0.6];
%atrue = [2,1,2.6,1.5,1.6,1.1,1.2,0.8,2.0,2.0,1.5,1.1];
%sigma = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np = length(xc);  % no.of parameters

k = 5; % controller gain

p0 = rand(2*N,1);  % initial positions at random

y0 = [p0];
tspan = [0 20];
RelTol = 1e-4;
AbsTol = 1e-4;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

%
%K = [N; k];
%
%[tout, yout] = ode45(@(t,y) sint_uniformspread(t,y,K,xborder,yborder),tspan,y0,options);
%
%posout = yout(end,:)';
%
%posout = [0.6979; 0.1904; 0.7717; 0.8539; 0.2114; 0.2847; 0.2670; 0.7709; 0.5200; 0.5571];
posout = [0.6979; 0.1904; 0.7717; 0.8539; 0.2114; 0.2847; 0.2670; 0.7709;];
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
	for j=1:length(xc)
		if(inpolygon(xc(j),yc(j),[vx vx(1)],[vy vy(1)]))
			ind{i} = [ind{i} j];
		end
	end
end

disp(ind);

%pause;

Km = zeros(np,np);
for i=1:np
	for j=1:i
		Km(i,j) = Kvalue(xc(i),yc(i),xc(j),yc(j),sigma);
		Km(j,i) = Km(i,j);
	end
end

%% First algorithm
%{

disp('First Algorithm');
pause;

a0 = 0.1.*ones(N*np,1);
gamma = 300;

nn = ((np*np - np)/2) + np;
Lambda0 = zeros(nn*N,1);
lambda0 = zeros(N*np,1);

y0 = [posout; Lambda0; lambda0; a0];
tspan = [0 1];
k2 = 1;
K = [N; k; np; sigma; gamma; k2];

flag = 1;
cntr_centre = ones(N,1);
c = zeros(N,1);
cx = zeros(N,1);
cy = zeros(N,1);
for i=1:N
	c(i) = ind{i}(cntr_centre(i));
	cx(i) = xc(c(i));
	cy(i) = yc(c(i));
end
e = 0.1;
stop = 0;
ctr = 1;
while(flag==1 || stop==0)

	if(flag==0)
		tspan = [tspan(1) tspan(1)+20];
		RelTol = 1e-4;
		AbsTol = 1e-4;
		options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
		stop = 1;
	end

	[tout{ctr}, yout{ctr}] = ode45(@(t,y) sint_adaptiveestimate1(t,y,K,cx,cy,xc,yc,sigma),tspan,y0,options);
	y0 = yout{ctr}(end,:)';
	f = 0;
	for i=1:N
		pi = y0(((i-1)*n)+1:i*n);
		chi = Km\(Kvector(pi(1),pi(2),xc,yc,sigma));
		tmp = 0;
		for j=1:np
			if(j~=c(i))
				tmp = tmp + abs(chi(j));
			end
		end
		if((abs(chi(c(i)))-tmp) > e)
			%disp(i);
			%disp(c(i));
			%pi
			%pause;
			if(cntr_centre(i)<length(ind{i}))
				cntr_centre(i) = cntr_centre(i) + 1;
				c(i) = ind{i}(cntr_centre(i));
				cx(i) = xc(c(i));
				cy(i) = yc(c(i));
			else
				f = f + 1;
				%disp('caught');
				%disp(f);
				%pause;
			end
		end

		if(f==N)
			flag = 0;
		end
	end

	tspan = [tspan(2) tspan(2)+0.2];

	ctr = ctr+1;

end

%}

% save('results1_sigma0-05_np100.mat');

% %{
%% Second algorithm

disp('Second Algorithm');
pause;

a0 = 0.1*ones(np,1);
gamma = 300;

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
cntr_centre = ones(N,1);
c = zeros(N,1);
cx = zeros(N,1);
cy = zeros(N,1);
sw = ones(N,1);  % filter switch

flag = 1;
for i=1:N
	c(i) = ind{i}(cntr_centre(i));
	cx(i) = xc(c(i));
	cy(i) = yc(c(i));
end
e = 0.1;
stop = 0;
ctr = 1;
while(flag==1 || stop==0)

	if(flag==0)
		tspan = [tspan(1) tspan(1)+20];
		stop = 1;
	end

	[tout{ctr}, yout{ctr}] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,cx,cy,npa,ind,xc,yc,sigma,sw),tspan,y0,options);
	%[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,xborder,yborder,xc,yc),tspan,y0,options);
	y0 = yout{ctr}(end,:)';
	f = 0;
	for i=1:N
		pi = y0(((i-1)*n)+1:i*n);
		Kvec = Kvector(pi(1),pi(2),xc,yc,sigma);
		Kmt = Km(ind{i},ind{i});
		chi = Kmt\(Kvec(ind{i}));
		tmp = 0;
		for j=1:npa(i)
			if((ind{i}(j))~=c(i))
				tmp = tmp + abs(chi(j));
			end
		end
		if((abs(chi(cntr_centre(i)))-tmp) > e)
			if(cntr_centre(i)<length(ind{i}))
				cntr_centre(i) = cntr_centre(i) + 1;
				c(i) = ind{i}(cntr_centre(i));
				cx(i) = xc(c(i));
				cy(i) = yc(c(i));
			else
				sw(i) = 0;
				f = f + 1;
			end
		end

		if(f==N)
			flag = 0;
		end
	end

	tspan = [tspan(2) tspan(2)+0.2];

	ctr = ctr+1;
end
% %}

%{
%% Modified second algorithm

disp('Modified second Algorithm');
pause;

a0 = 0.1*ones(N*np,1);
gamma = 300;

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
cntr_centre = ones(N,1);
c = zeros(N,1);
cx = zeros(N,1);
cy = zeros(N,1);
sw = ones(N,1);  % filter switch

flag = 1;
for i=1:N
	c(i) = ind{i}(cntr_centre(i));
	cx(i) = xc(c(i));
	cy(i) = yc(c(i));
end
e = 0.1;
stop = 0;
ctr = 1;
while(flag==1 || stop==0)

	if(flag==0)
		tspan = [tspan(1) tspan(1)+20];
		stop = 1;
	end


	[tout{ctr}, yout{ctr}] = ode45(@(t,y) sint_adaptiveestimate2mod(t,y,K,cx,cy,npa,ind,xc,yc,sigma,sw),tspan,y0,options);
	y0 = yout{ctr}(end,:)';
	f = 0;
	for i=1:N
		pi = y0(((i-1)*n)+1:i*n);
		Kvec = Kvector(pi(1),pi(2),xc,yc,sigma);
		Kmt = Km(ind{i},ind{i});
		chi = Kmt\(Kvec(ind{i}));
		tmp = 0;
		for j=1:npa(i)
			if((ind{i}(j))~=c(i))
				tmp = tmp + abs(chi(j));
			end
		end
		if((abs(chi(cntr_centre(i)))-tmp) > e)
			if(cntr_centre(i)<length(ind{i}))
				cntr_centre(i) = cntr_centre(i) + 1;
				c(i) = ind{i}(cntr_centre(i));
				cx(i) = xc(c(i));
				cy(i) = yc(c(i));
			else
				sw(i) = 0;
				f = f + 1;
			end
		end

		if(f==N)
			flag = 0;
		end
	end

	tspan = [tspan(2) tspan(2)+0.2];

	ctr = ctr+1;
end
%}

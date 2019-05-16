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

p = length(xc);  % no.of parameters

k = 1; % controller gain

p0 = rand(2*N,1);  % initial positions at random

y0 = [p0];
tspan = [0 20];
RelTol = 1e-4;
AbsTol = 1e-4;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

K = [N; k];

[tout, yout] = ode45(@(t,y) sint_uniformspread(t,y,K,xborder,yborder),tspan,y0,options);

posout = yout(end,:)';
pos = zeros(2,N);
for i=1:N
	pos(:,i) = posout(((i-1)*2)+1:i*2);
end
disp('Starting positions');
disp(pos);

%% Start estimation

for i=1:N		
	[vx,vy] = compute_voronoi(i,xborder,yborder,pos(1,:)',pos(2,:)');
	partitionx{i} = vx;
	partitiony{i} = vy;

	% locate which of the centres are in partition i
	ind{i} = [];
	for j=length(xc)
		if(inpolygon(xc(j),yc(j),vx,vy))
			ind{i} = [ind{i} j];
		end
	end
end

a0 = ones(N*p,1);
gamma = 1;

nn = ((p*p - p)/2) + p;
gamma0 = zeros(nn*N,1);
lambda0 = zeros(N*p,1);

y0 = [posout; gamma0; lambda0; a0];
tspan = [0 1];
K = [N; k; p; sigma; gamma];
cx = xc(1);
cy = yc(1);
[tout, yout] = ode45(@(t,y) sint_adaptiveestimate1(t,y,K,cx,cy),tspan,y0,options);
%[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,xborder,yborder,xc,yc),tspan,y0,options);

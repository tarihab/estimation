%% Code for estimation of scalar field

% Domain is unit square
xborder = [0 1 1 0];
yborder = [0 0 1 1];

N = 5; % no.of agents

p = 100;  % no.of parameters
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

k = 1; % controller gain

p0=[0.1010;
    0.1005;
    0.3920;
    0.1000;
    0.6020;
    0.1100;
    0.8820;
    0.1200;
    0.8930;
    0.5010];

y0 = [p0];
tspan = [0 30];
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
y0 = [yout(end,:)'; a0];
gamma = 1;
K = [N; k; p; sigma; gamma];

[tout, yout] = ode45(@(t,y) sint_adaptiveestimate1(t,y,K,xborder,yborder,xc,yc,ind),tspan,y0,options);
%[tout, yout] = ode45(@(t,y) sint_adaptiveestimate2(t,y,K,xborder,yborder,xc,yc),tspan,y0,options);

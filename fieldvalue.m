function phi = fieldvalue(x,y)

	phi =  3*(x).^2.*exp((-(x-0.7).^2 - (y-0.7).^2)/0.05) ...
  	+ 1.*exp((-(x-0.4).^2-(y-0.4).^2)/0.06) ...
   	+ 1/3*exp((-(x-0.2).^2 - (y-0.2).^2)/0.08);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% for exact parameterization case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xc = [0.2,0.35,0.6,0.85,0.7,0.75,0.15,0.35];
%yc = [0.25,0.26,0.18,0.3,0.75,0.9,0.75,0.6];
%atrue = [2,1,1.5,1.8,1.2,1.6,2.5,1.1];
%sigma = 0.1;
%xc = [0.1, 0.2, 0.35,0.4, 0.6, 0.7, 0.78,0.85,0.7, 0.6, 0.75,0.85,0.15,0.2, 0.32,0.35];
%yc = [0.18,0.25,0.26,0.32,0.18,0.25,0.3, 0.45,0.75,0.83,0.9, 0.60,0.65,0.72,0.75,0.6];
%atrue = [2,1,2.6,0.9,1.5,1.6,1.1,1.8,1.2,0.8,2.0,1.6,2.0,1.5,2.5,1.1];
%sigma = 0.1;
%xc = [0.1,0.2,0.4,0.6,0.7,0.85,0.7,0.6,0.85,0.15,0.2,0.35];
%yc = [0.18,0.25,0.32,0.18,0.3,0.45,0.75,0.9,0.60,0.65,0.72,0.6];
%atrue = [2,1,2.6,1.5,1.6,1.1,1.2,0.8,2.0,2.0,1.5,1.1];
%sigma = 0.1;

%phi = fieldestimate(x,y,xc,yc,sigma,atrue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

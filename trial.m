
%clear all
%close all

%[X,Y] = meshgrid(0:.005:0.2);
[X,Y] = meshgrid(0:.01:1);
%Z1 =  3*(1-X).^2.*exp(-(X.^2) - (Y+1).^2) ...
%   - 10*(X/5 - X.^3 - Y.^5).*exp(-X.^2-Y.^2) ...
%   - 1/3*exp(-(X+1).^2 - Y.^2);
%Z2 =  3*(X).^2.*exp((-(X-0.7).^2 - (Y-0.7).^2)/0.05) ...
%   + 1.*exp((-(X-0.4).^2-(Y-0.4).^2)/0.06) ...
%   + 1/3*exp((-(X-0.2).^2 - (Y-0.2).^2)/0.08);
%sigma = 0.03;
%Z3 = exp((-(X).^2 - (Y).^2)/(2*sigma*sigma));
% %Z3 = Z3.*(1/(sigma*sqrt(2*pi)));

[k,l] = size(X);


%Z4 = zeros(k,l);
% for first algorithm
%a = yout{end}(end,end-99:end)';

% compute and plot approximation
%for i=1:k
%	for j=1:l
%		Z4(i,j) = fieldestimate(X(i,j),Y(i,j),xc,yc,sigma,a);
%	end
%end
%mesh(X,Y,Z4);

%xc = [0.2,0.35,0.6,0.85,0.7,0.75,0.15,0.35];
%yc = [0.25,0.26,0.18,0.3,0.75,0.9,0.75,0.6];
%atrue = [2,1,1.5,1.8,1.2,1.6,2.5,1.1];
%sigma = 0.1;
%xc = [0.1,0.2,0.35,0.4,0.6,0.7,0.78,0.85,0.7,0.6,0.75,0.85,0.15,0.2,0.32,0.35];
%yc = [0.18,0.25,0.26,0.32,0.18,0.25,0.3,0.45,0.75,0.83,0.9,0.60,0.65,0.72,0.75,0.6];
%atrue = [2,1,2.6,0.9,1.5,1.6,1.1,1.8,1.2,0.8,2.0,1.6,2.0,1.5,2.5,1.1];
%sigma = 0.1;
%Z5 = zeros(k,l);
Z6 = zeros(k,l);
%a = yout{end}(end,end-15:end)';
a = parest(end,:)';
%a = parest{end}(end,:)';
for i=1:k
	for j=1:l
%		Z5(i,j) = fieldestimate(X(i,j),Y(i,j),xc,yc,sigma,atrue);
		Z6(i,j) = fieldestimate(X(i,j),Y(i,j),xc,yc,sigma,a);
	end
end
%
%mesh(X,Y,Z5);
figure;
mesh(X,Y,Z6);

%parest{4} = [yout{1}(:,end-99:end); yout{2}(:,end-99:end); yout{3}(:,end-99:end)];
%parest{3} = [yout{1}(:,end-15:end-8); yout{2}(:,end-15:end-8); yout{3}(:,end-15:end-8)];
%parest{2} = [yout{1}(:,end-23:end-16); yout{2}(:,end-23:end-16); yout{3}(:,end-23:end-16)];
%parest{1} = [yout{1}(:,end-31:end-24); yout{2}(:,end-31:end-24); yout{3}(:,end-31:end-24)];

% Algo 1
if(algo==1)
	ysize = length(yout);
	k1 = np-1;
	parest = cell(1,N);
	tset = [];
	for i=1:N
		parest{i} = [];
		k1 = (i*np)-1;
		k2 = (i-1)*np;
		for j=1:ysize
			parest{i} = [parest{i}; yout{j}(:,end-k1:end-k2)];
		end
	end
	for i=1:ysize
		tset = [tset; tout{i}];
	end
	avgperror = zeros(length(tset),np);
	for i=1:N
		avgperror = avgperror + abs(parest{i} - atrue);
	end
	avgperror = avgperror./N;
end

% Algo 2
%{
 	atmp = yout{end}(end,end-15:end)';
 	a = zeros(length(atmp),1);
 	s = 0;
 	for i=1:N
	 	for j=1:length(ind{i})
		 	idx = ind{i}(j);
		 	a(idx) = atmp(s+j);
	 	end
	 	s = s + length(ind{i});
 	end
%}

% Algo 2
if(algo==2)
	ysize = length(yout);
	parest = [];
	tset = [];
	k1 = np-1;
	%k2 = (i-1)*np;
	for j=1:ysize
		par_tmp = yout{j}(:,end-k1:end);
		[m_tmp,n_tmp] = size(par_tmp);
		par_tmp2 = zeros(m_tmp,n_tmp);
		s = 0;
		for i=1:N
			for k=1:length(ind{i})
				idx = ind{i}(k);
				par_tmp2(:,idx) = par_tmp(:,s+k);
			end
			s = s + length(ind{i});
		end
		parest = [parest; par_tmp2];
	end
	for i=1:ysize
		tset = [tset; tout{i}];
	end
	perror = zeros(length(tset),np);
	for i=1:N
		perror = abs(parest - atrue);
	end
end

% Algo 2 modified
if(algo==3)
	ysize = length(yout);
	parest = [];
	tset = [];
	k1 = np-1;
	%k2 = (i-1)*np;
	for j=1:ysize
		par_tmp = [];
		for i=1:N
			k1 = (i*np) - 1;
			k2 = (i-1)*np;
			k3 = sum(npa(1:i-1));
			k4 = sum(npa(1:i));
			par_tmp = [par_tmp, yout{j}(:,end-k1+k3:end-k1+k4-1)];
		end
		[m_tmp,n_tmp] = size(par_tmp);
		par_tmp2 = zeros(m_tmp,n_tmp);
		s = 0;
		for i=1:N
			for k=1:length(ind{i})
				idx = ind{i}(k);
				par_tmp2(:,idx) = par_tmp(:,s+k);
			end
			s = s + length(ind{i});
		end
		parest = [parest; par_tmp2];
	end
	for i=1:ysize
		tset = [tset; tout{i}];
	end
	perror = zeros(length(tset),np);
	for i=1:N
		perror = abs(parest - atrue);
	end
end

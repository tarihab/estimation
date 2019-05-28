
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
Z4 = zeros(k,l);

% for first algorithm
% a = yout{end}(end,end-195:end)';

% for second algorithm
 atmp = yout{end}(end,end-195:end)';
 a = zeros(length(atmp),1);
 s = 0;
 for i=1:N
	 for j=1:length(ind{i})
		 idx = ind{i}(j);
		 a(idx) = atmp(s+j);
	 end
	 s = s + length(ind{i});
 end

% compute and plot approximation
for i=1:k
	for j=1:l
		Z4(i,j) = fieldestimate(X(i,j),Y(i,j),xc,yc,sigma,a);
	end
end

mesh(X,Y,Z4);

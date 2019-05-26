function f = fieldestimate(x, y, xc, yc, sigma, a)
% evaluate the estimated value of the scalar field at point (x,y)

	p = length(xc);
	phi = zeros(1,p);
	
	f = 0;
	for i=1:p
		phi(i) = exp((-(x-xc(i))^2 - (y-yc(i))^2)/(2*sigma*sigma));
		%phi(i) = phi(i)*(1/(sigma*sqrt(2*pi)));

		phi(i) = a(i)*phi(i);
	end

	f = sum(phi);

end

function f = fieldestimate(x, y, xc, yc, sigma, a)
% evaluate the estimated value of the scalar field at point (x,y)

	p = length(xc);
	
	f = 0;
	for i=1:p
		phi = exp((-(x-xc)^2 - (y-yc)^2)/(2*sigma*sigma));
		phi = phi*(1/(sigma*sqrt(2*pi)));

		f = f + phi*a(i);
	end

end

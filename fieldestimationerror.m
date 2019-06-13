function e = fieldestimationerror(x, y, data)
% evaluate the estimated value of the scalar field at point (x,y)

	xc = data{1};
	yc = data{2};
	sigma = data{3};
	a = data{4};

	p = length(xc);
	phi = zeros(1,p);
	
	e = 0;
	for i=1:p
		phi(i) = exp((-(x-xc(i))^2 - (y-yc(i))^2)/(2*sigma*sigma));
		%phi(i) = phi(i)*(1/(sigma*sqrt(2*pi)));

		phi(i) = a(i)*phi(i);
	end

	e = abs(sum(phi) - fieldvalue(x,y));

end

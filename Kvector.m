function K = Kvector(x, y, xc, yc, sigma)

	np = length(xc);

	K = zeros(np,1);
	for i=1:np
		K(i) = exp((-(x-xc(i))^2 - (y-yc(i))^2)/(2*sigma*sigma));
		%K(i) = phi(i)*(1/(sigma*sqrt(2*pi)));
	end

end

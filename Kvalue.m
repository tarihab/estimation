function y = Kvalue(x1, y1, x2, y2, sigma)

	y = exp((-(x1-x2)^2 - (y1-y2)^2)/(2*sigma*sigma));
	y = y*(1/(sigma*sqrt(2*pi)));

end

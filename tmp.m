
%x = 10;
%k = func2(x,@func1);
%disp(k);

xborder = [0,2,2];
yborder = [0,0,1];
%xborder = [0,2,2,0];
%yborder = [0,0,1,1];

intv = polyintegrate_gl(xborder,yborder,@integrand_new,[],256);
disp(intv);

function [val] = integrand_const(x, y, data)
 
	val = 1.0;
end

function [val] = integrand_exp(x, y, data)

	val = exp(10*x)*tan(y);
end

function [val] = integrand_new(x, y, data)

	val = x*y*y;
end

function [y] = func1(x)
	y = x+2;
end

function [y] = func2(a,b)
	y = a + b(a+5);
end

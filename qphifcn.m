% computes the value of the density function at point q and parameter a

function qphi = qphifcn(qx,qy,a)

	[nr,nc] = size(qx);
	ind = a;

	phi = ones(nr,nc);
	qphi = zeros(nr,nc);
	for i=1:nr
	 for j=1:nc
	  q = [qx(i,j); qy(i,j)];
	  qphi(i,j) = q(ind).*phi(i,j);
	 end
	end

end

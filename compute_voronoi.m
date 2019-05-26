%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB code for computing the voronoi partition of the j-th site in a bounded polygon %%
%% Author: Rihab Abdul Razak, IITB-Monash Research Academy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xvert, yvert] = compute_voronoi(j, xborder, yborder, xsites, ysites)
% xborder - the x-coordinates of the boundary
% yborder - the y-coordinates of the boundary
% xsites - the x-coordinates of the sites
% ysites - the y-coordinates of the sites
% xvert - gives the x-coordinates of the voronoi cell of site j
% yvert - gives the y-coordinates of the voronoi cell of site j
% loop this function over j to get all the voronoi cells

	myx = xsites(j);
	myy = ysites(j);

	n = length(xsites);
	if(length(ysites) ~= n)
		disp('The length of xsites and ysites must be same..');
		xvert = [];
		yvert = [];
		return;
	end

	sortedxsites = xsites;
	sortedysites = ysites;

	% erase the j-th element
	sortedxsites(j) = [];
	sortedysites(j) = [];

	% sort elements based on distance to site j
	for i=1:n-1
		distancetoagents(i) = (myx-sortedxsites(i))^2 + (myy-sortedysites(i))^2;
	end
	[sortedxsites, sortedysites] = sortsites(j, n, distancetoagents, sortedxsites, sortedysites);

	borderx = xborder;
	bordery = yborder;

	% curr and bis are structures which represent a line segment between (x1,y1) and (x2,y2)
	curr = struct('x1',[],'y1',[],'x2',[],'y2',[]);
	%bis = struct('x1',[],'y1',[],'x2',[],'y2',[]);

	curr.x1 = myx;
	curr.y1 = myy;

	% loop over the sites other than j
	for i=1:n-1
	
		xcur = sortedxsites(i);
		ycur = sortedysites(i);

		curr.x2 = xcur;
		curr.y2 = ycur;
		[bis] = perpbisector(curr);
		[borderx, bordery] = lineintersection(borderx, bordery, bis, myx, myy);
	
	end

	xvert = borderx;
	yvert = bordery;

end

function [sx, sy] = sortsites(j, n, distancetoagents, sortedxsites, sortedysites)

	sx = sortedxsites;
	sy = sortedysites;

	% sortedxsites and sortedysites are of length (n-1)
	for i=1:(n-2)
		for k=(i+1):(n-1)
			if(distancetoagents(i)>distancetoagents(k))
				temp = distancetoagents(i);
				distancetoagents(i) = distancetoagents(k);
				distancetoagents(k) = temp;

				temp = sx(i);
				sx(i) = sx(k);
				sx(k) = temp;

				temp = sy(i);
				sy(i) = sy(k);
				sy(k) = temp;
			end
		end
	end

end


function [bis] = perpbisector(line1)

	TOL = 0.0000001;

	xnew1 = (line1.x1 + line1.x2)/2;
	ynew1 = (line1.y1 + line1.y2)/2;
	
	% slope of the perpendicular bisector
	if(((line1.x2 - line1.x1)<TOL) && ((line1.x2 - line1.x1)>-TOL))
		m = 0;
		xnew2 = xnew1 + 0.5;
		ynew2 = ynew1;
	else
		if(((line1.y2 - line1.y1)<TOL) && ((line1.y2 - line1.y1)>-TOL))
			% slope zero for line1 and consequently the perp bisector has infinity
			ynew2 = ynew1 + 0.5;
			xnew2 = xnew1;
		else
			m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
			m = -1/m1;
			if((m<=1.0) && (m>=-1.0))
				xnew2 = xnew1 + 0.5;
				ynew2 = (xnew1*(line1.x1-line1.x2) + ynew1*(line1.y1-line1.y2) - xnew2*(line1.x1-line1.x2))/(line1.y1-line1.y2);
			else
				ynew2 = ynew1 + 0.5;
				xnew2 = ((ynew2-ynew1) + m*xnew1)/m;
			end
		end
	end
	
	bis = struct('x1',xnew1,'y1',ynew1,'x2',xnew2,'y2',ynew2);

end

% intersection of a line with polygon defined by its vertices
function [borderx, bordery] = lineintersection(borderx, bordery, line, myx, myy)

	TOL = 0.0000001;
	condition = 0;
	% case where slope is infinity
	if(((line.x2 - line.x1)<TOL) && ((line.x2 - line.x1)>-TOL)) 
		condition = 1;
	end
	
	
	if(condition==1) 
		if((myx-line.x1)>0) 
			side = 1;
		else 
			side = -1;
		end
	else
	
		% slope of the line
		m = (line.y2 - line.y1)/(line.x2 - line.x1);

		% determine which side of the line the point is
		if(((myy-line.y1) - m*(myx-line.x1)) > 0) 
			side = 1;
		else 
			side = -1;
		end
	end

	% std::vector<double> borderxtmp(borderx);
	% std::vector<double> borderytmp(bordery);
	borderxtmp = borderx;
	borderytmp = bordery;

	% n = static_cast<int>(borderx.size());
	n = length(borderx);
	% determine the polygon vertices which lie on the same side as the point
	% std::vector<int> ind(n,0);
	ind = zeros(n,1);

	nochange = 1;
	for i=1:n
	
		currx = borderx(i);
		curry = bordery(i);
		switch side	
		 
			case 1 
				if(condition==1)  
					if((currx-line.x1)>=0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				else 
					if(((curry-line.y1) - m*(currx-line.x1)) >= 0) 
						ind(i) = 1;
					else  
						nochange = 0;
					end
				end

			case -1 
				if(condition==1) 
				
					if((currx-line.x1)<=0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				else 
					if(((curry-line.y1) - m*(currx-line.x1)) <= 0) 
						ind(i) = 1;
					else 
						nochange = 0;
					end
				end
		end
	end

	% if no changes to be done, return from the function
	if(nochange==1)
		return
	end

	%int i1,i2;
	%double ix, iy;
	%struct linedesc line2;
	line2 = struct('x1',[],'y1',[],'x2',[],'y2',[]);
	if(ind(1)==0)
	
		for i=1:n-1
		
			if((ind(i)==0) && (ind(i+1)==1)) 
			
				i1 = i;

				% compute intersection
				line2.x1 = borderx(i);
				line2.y1 = bordery(i);
				line2.x2 = borderx(i+1);
				line2.y2 = bordery(i+1);
				[ix,iy] = computeintersection(line, line2);
				
				% update the boundary variables
				borderxtmp(i) = ix;
				borderytmp(i) = iy;
			end
			if((ind(i)==1) && (ind(i+1)==0))
			
				i2 = i+1;

				% compute intersection
				line2.x1 = borderx(i);
				line2.y1 = bordery(i);
				line2.x2 = borderx(i+1);
				line2.y2 = bordery(i+1);
				[ix,iy] = computeintersection(line, line2);
				
				% update the boundary variables
				borderxtmp(i+1) = ix;
				borderytmp(i+1) = iy;
			end
		end
		if(ind(n)==1) 
		
			i2 = 1;
	
			% compute intersection
			line2.x1 = borderx(n);
			line2.y1 = bordery(n);
			line2.x2 = borderx(1);
			line2.y2 = bordery(1);
			[ix,iy] = computeintersection(line, line2);
				
			% update the boundary variables
			if(i1==i2)
				%borderxtmp.push_back(ix);
				%borderytmp.push_back(iy);
				borderxtmp(end+1) = ix;
				borderytmp(end+1) = iy;
			else
				borderxtmp(1) = ix;
				borderytmp(1) = iy;
			end

			if((i2+1)<(i1)) 
			%borderxtmp.erase(borderxtmp.begin() + (i2 + 1), borderxtmp.begin() + (i1));
			%borderytmp.erase(borderytmp.begin() + (i2 + 1), borderytmp.begin() + (i1));
			borderxtmp(i2+1:i1-1) = [];
			borderytmp(i2+1:i1-1) = [];
			end
		else
			if(1<(i1)) 
			%borderxtmp.erase(borderxtmp.begin(), borderxtmp.begin() + (i1));
			%borderytmp.erase(borderytmp.begin(), borderytmp.begin() + (i1));
			borderxtmp(1:i1-1) = [];
			borderytmp(1:i1-1) = [];
			end

			%if((i2+1)<=(n-1)) 
			if((i2+1)<=n) 
			%borderxtmp.erase(borderxtmp.begin() + (i2+1), borderxtmp.end());
			%borderytmp.erase(borderytmp.begin() + (i2+1), borderytmp.end());
			borderxtmp(i2+1:end) = [];
			borderytmp(i2+1:end) = [];
			end
		end
	else 
		if(ind(n)==1)
		
			for i=1:n-1
			
				if((ind(i)==1) && (ind(i+1)==0)) 
				
					i1 = i+1;

					% compute intersection
					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					borderxtmp(i+1) = ix;
					borderytmp(i+1) = iy;
				end
				if((ind(i)==0) && (ind(i+1)==1)) 

					i2 = i;

					% compute intersection
					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					if(i1==i2) 
						%borderxtmp.push_back(ix);
						%borderytmp.push_back(iy);
						%borderxtmp.insert(borderxtmp.begin()+(i+1),ix);
						%borderytmp.insert(borderytmp.begin()+(i+1),iy);
						borderxtmp = [borderxtmp(1:i)  ix  borderxtmp(i+1:end)];
						borderytmp = [borderytmp(1:i)  iy  borderytmp(i+1:end)];
					else 
						borderxtmp(i) = ix;
						borderytmp(i) = iy;
					end
				end
			end

			% remove elements from i1 to i2
			if((i1+1)<(i2)) 
				%borderxtmp.erase(borderxtmp.begin() + (i1+1), borderxtmp.begin() + (i2));
				%borderytmp.erase(borderytmp.begin() + (i1+1), borderytmp.begin() + (i2));
				borderxtmp(i1+1:i2-1) = [];
				borderytmp(i1+1:i2-1) = [];
			end

		else

			for i=1:n-1

				if((ind(i)==1) && (ind(i+1)==0)) 

					i1 = i+1;

					% compute intersection
					line2.x1 = borderx(i);
					line2.y1 = bordery(i);
					line2.x2 = borderx(i+1);
					line2.y2 = bordery(i+1);
					[ix,iy] = computeintersection(line, line2);
				
					% update the boundary variables
					borderxtmp(i+1) = ix;
					borderytmp(i+1) = iy;
				end
			end

			i2 = n;

			% compute intersection
			line2.x1 = borderx(1);
			line2.y1 = bordery(1);
			line2.x2 = borderx(n);
			line2.y2 = bordery(n);
			[ix,iy] = computeintersection(line, line2);
				
			% update the boundary variables
			if(i1==i2) 
				%borderxtmp.push_back(ix);
				%borderytmp.push_back(iy);
				borderxtmp(end+1) = ix;
				borderytmp(end+1) = iy;
			else 
				%borderxtmp[n-1] = ix;
				%borderytmp[n-1] = iy;
				borderxtmp(n) = ix;
				borderytmp(n) = iy;

				% remove elements from i1 to i2
				if((i1+1)<(i2))
					%borderxtmp.erase(borderxtmp.begin()+i1+1,borderxtmp.begin()+i2);
					%borderytmp.erase(borderytmp.begin()+i1+1,borderytmp.begin()+i2);
					borderxtmp(i1+1:i2-1) = [];
					borderytmp(i1+1:i2-1) = [];
				end
			end
		

		end
	end
	borderx = borderxtmp;
	bordery = borderytmp;

end

%% find intersection of two lines assuming they intersect ie. they are not parallel
function [ix, iy] = computeintersection(line1, line2)

	TOL = 0.0000001;

	condition1 = 0;
	condition2 = 0;

	% check if lines have slope infinity
	if(((line1.x2 - line1.x1)<TOL) && ((line1.x2 - line1.x1)>-TOL))
		condition = 1;
	end
	if(((line2.x2 - line2.x1)<TOL) && ((line2.x2 - line2.x1)>-TOL))
		condition2 = 1;
	end

	if(condition1==1 || condition2==1)
	
		% both conditions will not be 1 at the same time since it is assumed that the lines are not parallel
		if(condition1==1)
			m2 = (line2.y2 - line2.y1)/(line2.x2 - line2.x1);
			ix = line1.x1;
			iy = line2.y1 + m2*(ix-line2.x1);
		else % condition2 = 1
			m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
			ix = line2.x1;
			iy = line1.y1 + m1*(ix-line1.x1);
		end
	else
		% slopes of the two lines
		m1 = (line1.y2 - line1.y1)/(line1.x2 - line1.x1);
		m2 = (line2.y2 - line2.y1)/(line2.x2 - line2.x1);

		% computing the intersection
		ix = (m2*line2.x1 - m1*line1.x1 + line1.y1 - line2.y1)/(m2-m1);
		iy = line1.y1 + m1*(ix-line1.x1);
	end

end

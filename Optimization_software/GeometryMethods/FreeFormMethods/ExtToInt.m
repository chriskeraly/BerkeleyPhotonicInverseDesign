

function intParam = ExtToInt(obj,extParam,extDir,extPW)
	% obj = object (type Geometry, LevelSet, Freeform) with defined external coordinates
	% extParam = vector of values to convert
	% extDir = vector of directions of each (i.e. [x,y,x,z] would be [1,2,1,3])
	% extPW = vector of whether point or width to change ( [1,0,1] would be [P,W,P] )
	intParam = ...
		( extParam / obj.dx - (obj.x0/obj.dx - 1)*extPW ) .* (extDir==1) + ...
		( extParam / obj.dx - (obj.y0/obj.dx - 1)*extPW ) .* (extDir==2) + ...
		( extParam / obj.dx - (obj.z0/obj.dx - 1)*extPW ) .* (extDir==3);
	
	if(obj.testFlag)
		% fprintf('Converted external parameter %g to internal parameter %g \n',[extParam; intParam]);
	end
end
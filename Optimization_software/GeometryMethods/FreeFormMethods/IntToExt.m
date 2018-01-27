

function extParam = IntToExt(obj,intParam,intDir,intPW)
	% obj = object (type Geometry, LevelSet, Freeform) with defined external coordinates
	% intParam = vector of values to convert
	% intDir = vector of directions of each (i.e. [x,y,x,z] would be [1,2,1,3])
	% intPW = vector of whether point or width to change ( [1,0,1] would be [P,W,P] )
	extParam = ...
		( intParam * obj.dx + (obj.x0 - obj.dx)*intPW ) .* (intDir==1) + ...
		( intParam * obj.dx + (obj.y0 - obj.dx)*intPW ) .* (intDir==2) + ...
		( intParam * obj.dx + (obj.z0 - obj.dx)*intPW ) .* (intDir==3);
	
	if(obj.testFlag)
		% fprintf('Converted internal parameter %g to external parameter %g \n',[intParam; extParam]);
	end
end
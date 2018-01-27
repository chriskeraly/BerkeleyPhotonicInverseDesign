% DO NOT CHANGE: Main material setup file

if(strcmp(shapeType,'Polygon'))
   matArr = {opt.geo.newShapeEps, opt.geo.epsClad};
   for i=1:opt.geo.numShapes
       matArr = [matArr {opt.geo.shapes{i}.eps}];
   end
   save('SetupMaterials.mat','matArr','freqVec');
   [status,result] = system([lumerical ' -run LumericalMethods/SetupMaterials.lsf']);
   load('SetupMaterials.mat');
   opt.grad.epsVec(freqInd) = epsArr{1};
   opt.grad.epsOutVec(freqInd) = epsArr{2};
   for i=1:opt.geo.numShapes
       opt.geo.shapes{i}.epsVec = epsArr{2+i};
   end
else
   [~, eps, epsOut, ~, ~, ~, ~, ~] = opt.geo.getGeometry;
   matArr = {eps, epsOut};
   save('SetupMaterials.mat','matArr','freqVec');
   [status,result] = system([lumerical ' -run LumericalMethods/SetupMaterials.lsf']);
   load('SetupMaterials.mat');
   opt.grad.epsVec(freqInd) = epsArr{1};
   opt.grad.epsOutVec(freqInd) = epsArr{2};
end

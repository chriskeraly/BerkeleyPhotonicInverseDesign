function [epsGridImport,x_min,x_span,y_min,y_span,z_center,z_span] = autoImportGeometry(eps, freq, lumerical, baseFile, velocityMon, queueName, dx)

fprintf(['Import Geometry (' baseFile ')...\n']); %multiWaitbar(['Import Geometry (' baseFile ')'],'Busy');

freq=freq(1);

done=0;
save('autoImportGeometry.mat','done','eps','freq','baseFile','velocityMon','queueName');
runLumericalScript(lumerical, 'autoImportGeometry.mat', 'LumericalMethods/autoImportGeometry.lsf');

load('autoImportGeometry.mat','eps','epsGrid','x','y','vel_x_min','vel_x_span','vel_y_min','vel_y_span','vel_z','vel_z_span');

[x_grid, y_grid] = meshgrid(x,y);

x_min = vel_x_min;
x_span = vel_x_span;
y_min = vel_y_min;
y_span = vel_y_span;
z_center = vel_z;
z_span = vel_z_span;
numX = floor(1e-4*round(1e4*x_span / dx));
numY = floor(1e-4*round(1e4*y_span / dx));
numZ = floor(1e-4*round(1e4*z_span / dx));

[xGridImport, yGridImport] = meshgrid(x_min+dx*(0:numX-1)+dx/2, y_min+dx*(0:numY-1)+dx/2);
epsGridImport = interp2(x_grid, y_grid, epsGrid, xGridImport, yGridImport);
epsGridImport = 1*(epsGridImport==eps);

%multiWaitbar(['Import Geometry (' baseFile ')'],'Close');

end

%% calcTheta(grid)
% grid = binary bitmap of shapes
function theta = calcTheta(grid)
[Ny, Nx] = size(grid);

%flipud because y increases downward
filterSize = 11; % must be ODD
filterMid = (filterSize+1)/2;

[filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
filterX = filterX - filterMid;
filterY = filterY - filterMid;
filter = exp(1i*angle(filterX+1i*filterY));
filter(filterMid,filterMid)=0;
%filter = filter .* (abs(filterX+1i*filterY)<filterMid);

theta = mod(2*pi+angle(1e-4*round(1e4*conv2(filter,1*grid))),2*pi);
theta = theta(filterMid:filterMid+Ny-1, filterMid:filterMid+Nx-1).*(getBoundaryIn(grid)+getBoundaryOut(grid));
end
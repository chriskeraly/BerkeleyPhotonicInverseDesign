%% getBoundaryIn(grid)
% grid = binary bitmap of shapes
function gridIn = getBoundaryIn(grid)
[Ny, Nx] = size(grid);

epsilon_pad=zeros(Ny+2, Nx+2);

epsilon_pad(2:Ny+1,2:Nx+1)=(grid==1);

if(Ny==1)
    epsilon_pad2=epsilon_pad(2:Ny+1,2:Nx+1)...
        &epsilon_pad(2:Ny+1,1:Nx)...
        &epsilon_pad(2:Ny+1,3:Nx+2);
else
    epsilon_pad2=epsilon_pad(2:Ny+1,2:Nx+1)&epsilon_pad(1:Ny,1:Nx)...
        &epsilon_pad(1:Ny,3:Nx+2)&epsilon_pad(3:Ny+2,3:Nx+2)...
        &epsilon_pad(3:Ny+2,1:Nx)&epsilon_pad(1:Ny,2:Nx+1)...
        &epsilon_pad(2:Ny+1,1:Nx)&epsilon_pad(3:Ny+2,2:Nx+1)...
        &epsilon_pad(2:Ny+1,3:Nx+2);
end

gridIn = epsilon_pad(2:Ny+1,2:Nx+1) - epsilon_pad2;
end
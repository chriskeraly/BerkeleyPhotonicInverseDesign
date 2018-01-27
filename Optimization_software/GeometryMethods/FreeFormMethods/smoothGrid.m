%% 1*dx SMOOTHING FUNCTION
% gets rid of tiny 1 mesh cell features
function gridNew = smoothGrid(grid, flag)
[Ny, Nx] = size(grid);

grid_pad=ones(size(grid)+2);

for i=1:10
    grid_pad(2:Ny+1,2:Nx+1)=grid;
    
    grid_pad2=grid_pad(1:Ny,2:Nx+1)...
        +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
        +grid_pad(2:Ny+1,3:Nx+2);
    
    if(flag>0)
        gridNew = 1*(grid | (grid_pad2>2));
    else
        gridNew = 1*(grid & (grid_pad2>1));
    end
    
    grid_pad=ones(size(grid)+2);
    grid_pad(2:Ny+1,2:Nx+1)=gridNew;
    
    grid_pad2=grid_pad(1:Ny,2:Nx+1)...
        +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
        +grid_pad(2:Ny+1,3:Nx+2);
    
    if(flag<0)
        gridNew = 1*(gridNew | (grid_pad2>2));
    else
        gridNew = 1*(gridNew & (grid_pad2>1));
    end
    
    grid=gridNew;
end

end
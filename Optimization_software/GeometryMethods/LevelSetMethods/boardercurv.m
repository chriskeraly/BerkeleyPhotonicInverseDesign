function sidemask=boardercurv(phi,newShapeRad,dx)


%% Boarder management

shape=phi<0;

newShapeRad=floor(newShapeRad/dx);


[Ny,Nx]=size(shape);
sidemask=zeros(Ny,Nx);

currentstate=0;
transitionpoints=[];

for i=1:Ny
    if shape(i,1)~=currentstate;
        transitionpoints(end+1)=i;
        currentstate=shape(i,1);
    end
end

for k=1:length(transitionpoints);
    
    for i=1:Ny;
        for j=1:Nx;
            if ((i-transitionpoints(k)-newShapeRad)^2+j^2)<(newShapeRad^2)
                sidemask(i,j)=-2*(shape(i,1)-0.5);
            end
            if (i-transitionpoints(k)+newShapeRad)^2+j^2<newShapeRad^2
                sidemask(i,j)=-2*(shape(i,1)-0.5);
            end
        end
    end
end

currentstate=0;
transitionpoints=[];
    
for i=1:Nx
    if shape(1,i)~=currentstate;
        transitionpoints(end+1)=i;
        currentstate=shape(1,i);
    end
end

for k=1:length(transitionpoints);
    
    for i=1:Nx;
        for j=1:Ny;
            if ((i-transitionpoints(k)-newShapeRad)^2+j^2)<(newShapeRad^2)
                sidemask(j,i)=-2*(shape(1,i)-0.5);
            end
            if (i-transitionpoints(k)+newShapeRad)^2+j^2<newShapeRad^2
                sidemask(j,i)=-2*(shape(1,i)-0.5);
            end
        end
    end
end

currentstate=0;
transitionpoints=[];


for i=1:Ny
    if shape(i,Nx)~=currentstate;
        transitionpoints(end+1)=i;
        currentstate=shape(i,Nx);
    end
end

for k=1:length(transitionpoints);
    
    for i=1:Ny;
        for j=1:Nx;
            if ((i-transitionpoints(k)-newShapeRad)^2+(j-Nx)^2)<(newShapeRad^2)
                sidemask(i,j)=-2*(shape(i,Nx)-0.5);
            end
            if ((i-transitionpoints(k)+newShapeRad)^2+(j-Nx)^2)<newShapeRad^2
                sidemask(i,j)=-2*(shape(i,Nx)-0.5);
            end
        end
    end
end

currentstate=0;
transitionpoints=[];
    
for i=1:Nx
    if shape(Ny,i)~=currentstate;
        transitionpoints(end+1)=i;
        currentstate=shape(Ny,i);
    end
end



for k=1:length(transitionpoints);
    
    for i=1:Nx;
        for j=1:Ny;
            if ((i-transitionpoints(k)-newShapeRad)^2+(j-Ny)^2)<(newShapeRad^2)
                sidemask(j,i)=-2*(shape(Ny,i)-0.5);
            end
            if ((i-transitionpoints(k)+newShapeRad)^2+(j-Ny)^2)<(newShapeRad^2)
                sidemask(j,i)=-2*(shape(Ny,i)-0.5);
            end
        end
    end
end

sidemask=-sidemask;

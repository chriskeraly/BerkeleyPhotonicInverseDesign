dx=1;
xvec=[-50:dx:50];
yvec=xvec';

[xgrid,ygrid]=meshgrid(xvec,yvec);

phi0=sqrt(xgrid.^2+ygrid.^2)-30;
%phi0=ygrid.^2-25^2;
figure(1);imagesc(phi0);
figure(2);imagesc(phi0<0);

vel=ones(size(xgrid));
vel=(sin(-xgrid/150*8*pi)+sin(-ygrid/150*8*pi)+0.3)/2.3;%
%vel=sign(xgrid);
figure(4);imagesc(vel);
alpha=1;
phi=phi0;
r0=13;
vel
b=1.011*r0;
sidemask=boardercurv(phi,r0,dx);
vel=vel.*(sidemask==0)+sidemask*max(max(abs(vel)));
vel=vel/max(max(abs(vel)));

dt=alpha/(2*max(max(abs(vel)))/dx+4*b/dx^2);






for i=1:3500
    
    
    vPhi=vDotPhi(phi,vel,dx);
    curv=(abs(curvterm(phi,dx))>1/r0)*1/r0.*sign(curvterm(phi,dx));
    gradphi=sqrt(FODiff(phi,'x','central',dx).^2+FODiff(phi,'y','central',dx).^2);
    
    
    
    
    phi=phi-dt*(vPhi-b*gradphi.*curv);
    
    if floor(i/100)==i/100
        phi=signedDist(phi,dx);
        lastreinit=phi;
    end
    
    figure(3);imagesc(phi); pause(0.01);
    %figure(4);plot(phi(:,95));pause(0.1);
end

    
    
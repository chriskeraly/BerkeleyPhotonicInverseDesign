function phi=imposecurve(phi,newShapeRad)

figure(666)
phi=reinit(phi);
deltaT=0.05;
sidemask=boardercurv(phi,newShapeRad,1);
for i=1:10000
    curvTerm=curvaturechris(phi);
    curvTerm=curvTolFunc(curvTerm,newShapeRad);
    
    gradPhi=sqrt(diffx(phi).^2+diffy(phi).^2);
    phi = phi - deltaT*(-gradPhi.*(curvTerm+sidemask));
    
    if floor(i/100)==i/100; phi=reinit(phi);figure(200);imagesc(phi>0);pause(0.01);end
end
imagesc(phi);
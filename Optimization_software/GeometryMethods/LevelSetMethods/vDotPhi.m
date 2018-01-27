function [ vdotphi ] = vDotPhi( phi, V ,dx)

V=V+(V==0)*10^-99;

gradPhix=(min(V.*FODiff(phi,'x','up',dx),0)+max(V.*FODiff(phi,'x','down',dx),0))./V;
gradPhiy=(min(V.*FODiff(phi,'y','up',dx),0)+max(V.*FODiff(phi,'y','down',dx),0))./V;

vdotphi=V.*sqrt(gradPhix.^2+gradPhiy.^2);

%vdotphi=min(V.*FODiff(phi,'x','up',dx),0).*sign(FODiff(phi,'x','up',dx))+max(V.*FODiff(phi,'x','down',dx),0).*sign(FODiff(phi,'x','down',dx))+min(V.*FODiff(phi,'y','up',dx),0).*sign(FODiff(phi,'y','up',dx))+max(V.*FODiff(phi,'y','down',dx),0).*sign(FODiff(phi,'y','down',dx));


end


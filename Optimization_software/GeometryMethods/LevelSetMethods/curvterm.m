function [ curv ] = curvterm( phi,dx )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

eps=10^-99;
gradientmag= sqrt(FODiff(phi,'x','central',dx).^2+FODiff(phi,'y','central',dx).^2+eps);

curv=(FODiff(phi,'x','central',dx).^2.*FODiff(FODiff(phi,'y','down',dx),'y','up',dx)-2*FODiff(phi,'x','central',dx).*FODiff(phi,'y','central',dx).*FODiff(FODiff(phi,'x','central',dx),'y','central',dx)+FODiff(phi,'y','central',dx).^2.*FODiff(FODiff(phi,'x','down',dx),'x','up',dx))./(gradientmag.^3);

end


function G=energy(B,Ban,K2toK1,orientTheta,orientPhi,theta,phi)
    G = -B .* (sin(orientTheta).*sin(theta).*cos(orientPhi-phi) + cos(orientTheta).*cos(theta)) + .5*Ban.*(sin(theta).^2 + K2toK1.*sin(theta).^4);

end
R = 20;
% for ring
L = 2 * R * pi;
% cavity
%L = 10 * 2;
neff = 2.61189; %3.45;

lambda0 = 1.55;

m = round(neff*L/lambda0);

lambda1 = (neff*L)/m;
lambda2 = (neff*L)/(m+1);
dlambda = lambda1 - lambda2;


fprintf(1, 'lambda1=%f, lambda2=%f, dlambda=%f\n', lambda1, lambda2, dlambda);


delta = lambda0*lambda0 / (neff*L - lambda0);

fprintf(1, 'lambda0=%f, delta=%f\n', lambda0, delta);

FSR =lambda0/(L*neff)
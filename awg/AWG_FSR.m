pi = 3.141592654;
u0 = 4 * pi * 10^(-7);
e0 = 8.8854188 * 10^(-12);
wlc=1.550;
h=0.22;
neff_st=2.331302;
nu=1.444;
nc=3.472;
nl=nu;
a=h/2;
lu=3;
k0=2*pi/wlc;
ku=k0*sqrt(neff_st^2-nu^2);
kc=k0*sqrt(nc^2-neff_st^2);
kl=k0*sqrt(neff_st^2-nl^2);
we_te=a+1/(2*ku)+1/(2*kl);
theta0=wlc/nc/we_te*sqrt(2*pi);
theta_out_max=sqrt(lu*theta0^2/8.7);
theta_out_max


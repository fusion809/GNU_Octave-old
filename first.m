# Brenton Horne's script 20161204
clear all

# constants
N         = 1000;
NN        = 100000;
n         = 0:N;
a         = 0;
b         = 200;
k         = 1;
M         = round(N/6);
nsub      = n(2:N);
format long g;

# Basic matrices
x         = -cos(pi*n'/N);
xx        = linspace(-1,1,NN+1)';
yy        = (b-a)/2*xx+(a+b)/2;
xsub      = x(2:N);
y         = (b-a)/2*x+(a+b)/2;
ysub      = (b-a)/2*xsub+(a+b)/2;
T         = cos(acos(x)*n);
TT        = cos(acos(xx)*n);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dTsub     = Usub*diag(n);
dT        = [-((-1).^n).*n.^2; dTsub; n.^2];
d2Tsub    = diag(1./(1-xsub.^2))*(diag(xsub)*Usub-Tsub*diag(n))*diag(n);
d2T       = [((-1).^n).*(n.^2).*(n.^2-1)/3; d2Tsub; (n.^2).*(n.^2-1)/3];
D2        = d2T/T;
E2        = D2(2:N,2:N);

# Computation
# H is the matrix version of the Sturm-Liouville operator.
H         = - 4/((b-a)^2)*E2 + k*diag(ysub);
[Y, LAM]  = eig(H);
# Convert LAM (diagonal matrix) to vertical vector.
Lam       = diag(LAM);
[Lam, IX] = sort(Lam, 'ascend');
Y         = Y(:,IX);

# Add boundary values
Y         = [zeros(1,N-1); Y; zeros(1,N-1)];

# Interpolate
aa        = T\Y;
YY        = TT*aa;

# Error analysis
# the exact solution calls for the eigenvalues to negative the zeros of the Airy Ai function
err       = airy(0,-Lam);
rms       = sqrt(sum(err(1:M).^2)/(M));

# Plotting
figure(1);
plot(nsub,log10(err),"linewidth",2)
title("Log10 of Airy Ai of - the eigenvalues; lower (more negative) the better the solution")
figure(2);
plot(yy,YY(:,1),'-r',"linewidth",2)
title("Plot of the first eigenfunction on [0,200]")
figure(3);
plot(yy,YY(:,150),'-g',"linewidth",2)
title("Plot of the 150th eigenfunction on [0,200]")
figure(4);
plot(yy,YY(:,300),'-m',"linewidth",2)
title("Plot of the 300th eigenfunction on [0,200]")

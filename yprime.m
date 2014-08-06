function dydt = yprime(t,y)

dydt(1,1) = y(2);
dydt(2,1) = -0.2*y(2) - sin(y(1));

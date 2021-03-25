%C2 main
N=201;
epsilonlist = [1,0.1 ,0.001,0 ] ;
n= length(epsilonlist) ;
c=2;
a=0;
B=2*pi;
T=2;
TOL= 0.001 ;
errorlist = [] ;

for i = 1:length(epsilonlist)
epsilon = epsilonlist(i);
[u,u_star,x,t]= main_test2(2,N,a,B,c,TOL,epsilon,T); %anropas med 2 för C2boundary
figure(1)
subplot(n,1,i)
hold on;
plot(x,u(:,end)) %approximated solution

%plot(x, u_star(:,end)) %we don't want this for C2

xlabel('x')
ylabel('u approx')
title(num2str(epsilon))%different epsilonvalues
end
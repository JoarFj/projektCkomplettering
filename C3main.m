
N_list=[41,81,161,321,641];
n= length(N_list ) ;
c = 2;
a = 0;
B=2* pi ;
T=2;
TOL= 0.0001;
errorlist=[] ;
for i = 1:length(N_list )
N = N_list( i ) ;
h = (B-a)/(N-1) ;
epsilon = h/2 ;

[u,u_star,x,t]= main_test2(2,N,a,B,c,TOL,epsilon,T) ;

figure(1)
subplot(n,1,i )
hold on ;
plot(x,u(:,end))
%plot( x , u_star( : , end ))
title([ 'N = ' ,num2str(N), '    Epsilon=' ,num2str(epsilon)])
xlabel('x')
ylabel('u approx')
end

%C1 main
N_list=[41,81,161,321,641 ];
n=length(N_list ) ;
c=2;
a=-1;
B= 1;
epsilon=0.1 ;
T= 0.4;
TOL= 0.001;
errorlist=[];
 
for i = 1: length(N_list)
N = N_list(i) ;
[u,u_star,x,t]= main_test2(1,N,a,B,c,TOL,epsilon, T) ; %anropas med 1 för C1 boundary+function
error=norm(u(:,end) - u_star(:,end),2) ;
errorlist = [errorlist,error] ;
figure(1)
subplot(n,1,i)
hold on ;
plot(x,u(:,end))

plot(x,u_star(:,end))
title(num2str(N)) %different N
legend('u approx' ,'u exact' ) ;
end

figure(2)
plot( N_list , errorlist , ' x ' )
title( 'Errorplot as a function of Gridsize' )
xlabel( 'Grid size N' )
ylabel( 'l2-norm' )


%difference between two points (order of convergence):
q1 = log(errorlist(1)/errorlist(2)) /log(N_list(2)/N_list(1)) ;
q2 = log(errorlist(3)/errorlist(4)) /log(N_list(4)/N_list(3)) ;

disp( [ 'Convergence order: ' , num2str(q2) , '.' ] ) ;
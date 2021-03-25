function[UU, u_star , x , t ]= main_test2(C_info,N, a ,B,c ,TOL, epsilon , T)

x=linspace(a,B,N) ; %a,B decides spacedomain x
h=(B-a ) / ( N-1 ) ;
t=0;
k = 0.00001; %timestep

n=round(T/k ) ;

if C_info == 1
u0 = u_init(c,epsilon,x) ;% u0 for C1
end
if C_info==2
u0 = u_initC2(x); %u0 for C2-C3
end


u_exavect = zeros(N,N) ; %analytic solution vector
A = zeros(N) ;
A = A + diag( ones(N-1 , 1 ) , -1 ) * -1/ 2;
A = A + diag( ones(N-1 , 1 ) , 1 ) * 1 / 2 ;
A( 1 , 1 ) = -1/ 2;
A(N,N) = 1 / 2 ;
A = sparse(A) ;

M = eye(N) *2*h / 3 ;
M = M + diag( ones(N-1 , 1 ) , -1 ) *h / 6 ;
M = M + diag( ones(N-1 , 1 ) , 1 ) *h / 6 ;
M( 1 , 1 ) = h / 3 ;
M(N,N) =h / 3 ;
M = sparse(M) ;


S = eye(N) *2/ h ;
S = S + diag( ones(N-1 , 1 ) , -1 ) * -1/ h ;
S = S + diag( ones(N-1 , 1 ) , 1 ) * -1/ h ;
S( 1 , 1 ) = 1 / h ;
S(N,N) = 1/ h ;
S = sparse( S ) ;

U_list = zeros( length( x ) , length( t ) ) ;
for i = 1: length( x )
uu= u_exact( c , epsilon , x (i) , t ) ;
U_list( i , : ) =uu ;
end

u_vector = zeros(N, n ) ;
u_vector( : , 1 ) = u0' ;
 
for jj = 1 : n
U_temp = u_vector( : , jj ) ;

b = -A*U_temp .^2/2 - epsilon*S*U_temp ;

w1 = cg(M,b,TOL) ;
w = -A*( U_temp+k*w1 / 2 ).^ 2 / 2 - epsilon*S *( U_temp+k*w1 / 2 ) ;

w2 = cg(M,w,TOL) ;
w = -A*( U_temp+k*w2 / 2 ).^ 2 / 2 - epsilon*S *( U_temp+k*w2 / 2 ) ;

w3 = cg(M,w,TOL) ;
w = -A*( U_temp+k*w3 ) .^ 2 / 2 - epsilon *S *( U_temp+k*w3) ;

w4 = cg(M,w,TOL) ;
U_temp = U_temp+ k / 6 * ( w1 + 2*w2 + 2*w3 + w4 ) ;

u_vector(:,jj+1) = U_temp ;
t = t + k ;

if C_info ==1 %C1 boundary
u_vector( 1 , jj + 1 ) = u_exact ( c , epsilon , a , t ) ;
u_vector( end , jj + 1 ) = u_exact ( c ,epsilon ,B, t ) ;
end
 
if C_info == 2 %C2-C3 boundary
u_vector( 1 , jj + 1 ) = 0;
u_vector( end , jj + 1 ) = 0 ;
end

u_exavect(: , jj + 1 ) = u_exact( c ,epsilon , x , t ) ;

end
UU= u_vector ;
u_star=u_exavect ;

end

%u_exact used for boundaries in C1
function U = u_exact(c,eps,x,t)

U=c-tanh((x+1/2-c*t)/( 2* eps)) ;
 
 end
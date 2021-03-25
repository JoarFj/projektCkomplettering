%initial u for C1
function u = u_init(c,eps,x)

u = c-tanh((x+1/2)/(2*eps));

end
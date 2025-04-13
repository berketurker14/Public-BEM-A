function zeta = tellesTransformation(zeta0,gp)                
    lengp = length(gp);
    q=3;
    zeta = zeros(lengp,1);   
    for n=1:lengp
        t0=(pow((1.0+zeta0),1.0/q)-pow((1.0-zeta0),1.0/q))/(pow((1.0+zeta0),(1.0/q))+pow((1.0-zeta0),(1.0/q)));
        zeta(n,1)=zeta0+pow(2.0,-q)*pow((pow((1.0+zeta0),(1.0/q))+pow((1.0-zeta0),(1.0/q))),q)*pow((gp(n)-t0),q);
    end
end
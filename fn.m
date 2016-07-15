function k=fn(v,n,k1)
 %caculating kn
 %v is voltage n is memberane varity
 an=0.1*(1-0.1*v)/(exp(1-0.1*v)-1.0);
 bn=0.125*exp(-v/80.0);
 k=k1*(an*(1-n)-bn*n);
function k=fh(v,h,k1)
 %caculating kn
 %v is voltage n is memberane varity
 ah=0.07*exp(-v/20.0);
 bh=1.0/(exp(3.0-0.1*v)+1.0);
 k=k1*(ah*(1-h)-bh*h);
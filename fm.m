function k=fm(v,m,k1)
 %caculating km
 %v is voltage m is memberane varity
 am=(2.5-0.1*v)/(exp(2.5-0.1*v)-1);
 bm= 4*exp(-v/18.0);
 k=(am*(1-m)-bm*m)*k1;
 
 
function [kv,I_na,I_k,I_l,Ii]=fv1(v1,v,v3,ve1,ve,ve3,m,n,h,R,C,dC)
 %caculating  [kv,I_na,I_k,I_l,Ii]
 %v1=V_i-1,v=V_i,v3=V_i+1, m n h p is memberane varity,
 %ve1=Ve_i-1,ve=Ve_i,ve3=Ve_i+1
 %R is resistance
 v_l = 10.589; %in mV
 V_rest=-70; %in mv
 gl = 0.3; %in /kohm*cm^2 **also gmax_l
 %c = 1; %in microfarad/cm^2
 gna=120.0;%in /kohm*cm^2 **also gmax_l
 gk=36.0;%in /kohm*cm^2 **also gmax_l
 v_na=115.0 ;%in mV
 v_k=-12.0;%in mV
 I_na=gna*(m^3)*h*(v-v_na);
 I_k=gk*(n^4)*(v-v_k);
 I_l=gl*(v-v_l);
 Ii= I_na+I_k+I_l;
 Im=(v1-2*v+v3+ve1-2*ve+ve3)/R;
 kv=(Im-Ii-(v-70)*dC)/C;
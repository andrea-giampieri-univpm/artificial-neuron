syms Ic1 Iy1 Ix0 Iy0 C E Ic0 I1 Vc1 Vc0 Ik Ana0 Ana1 Ak0 Ak1 Ak1 T l y Ina Ia0 Ia1 Is0 Is1 e Ena Vna Inap1 Inap0 Iext

 eq1= (Vc1-Vc0)*C/T == -(Inap0+Ik+Ia-Iext);
 eq2= (Is1-Is0)/T == l*Ia0*(Vc0-y*Ia0); %rc indutt potassio
 eq3= (Ia1-Ia0)/T == l*Is0*(Vc0-y*Ia0); %rc indutt sodio
 eq4= e*(Inap1-Inap0)/T==Vc0-Ena-Vna; %topologia
 sol = solve([eq1,eq2,eq3,eq4],[Vc1,Ia1,Is1,Inap1]);

sol

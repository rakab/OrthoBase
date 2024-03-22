#-
Auto I a;
Auto I b;
F f1,f2;
CF P8;
CF SUNT,SUNF,SUND;

Tensor Tr(cyclic),Tp,f(antisymmetric),ff(rcyclic), d(symmetric);
Function TM,TT;
Symbols a,nf,Nc,Na,cF,cA,[cF-cA/6];
Dimension NA;
AutoDeclare  I mu=NA,mmu=NA,mnu=NA,nu=NA;
Dimension Nc;
AutoDeclare I i=Nc,mi=Nc,ii=Nc;
AutoDeclare I j=Nc,mj=Nc,jj=Nc;

*L t = distrib_(1,3,f1,f1,a1,a2,a3,a4,a5,a6);
*L t = distrib_(1,2,f1,f2,a1,a2,a3,a4);
*L t = f1(a1,b1)*f2(a2,b2);
*L t = f1(a1,b1)*f1(a2,b2)*f1(a3,b3);
*L t1 = e_(a1,a2,a3,a4)*e_(b1,b2,b3,b4);
L t2 = perm_(1,f1,a1,a2)*f2(b1,b2);
*L t3 = t2-t1;
*L t = perm_(f1,j1,j2)*f2(mj1,mj2)*perm_(f1,j3)*f2(mj3)*perm_(1,f1,mj1,mj3)*f2(jj1,jj3)*perm_(1,f1,mj2)*f2(jj2);
L P27(mu1,nu1,mu2,nu2)=P8(mu1,mmu1)*SUNT(mmu1,i1,j1)*SUNT(mu2,i2,j2)*P8(nu1,mnu1)*SUNT(mnu1,jj1,ii1)*SUNT(nu2,jj2,ii2)*perm_(f1,i1,i2)*f2(mi1,mi2)*perm_(1,f1,mi1)*f2(ii1)*perm_(1,f1,mi2)*f2(ii2)*perm_(f1,j1,j2)*f2(mj1,mj2)*perm_(1,f1,mj1)*f2(jj1)*perm_(1,f1,mj2)*f2(jj2);
*d_(mu1,mu2)*d(nu1,nu2);

L t1= 1/2*(perm_(f1,mu1,mu2)*f2(nu1,nu2)+perm_(f1,mu1,mu2)*f2(mmu1,mmu2)*Tr(mmu1,mnu2,mmu2,mnu1)*perm_(f1,mnu1,mnu2)*f2(nu1,nu2));
*d_(mu1,mu2)*d_(nu1,nu2);
L t2= 1/4*(perm_(f1,mu1,mu2)*f2(nu1,nu2)+perm_(f1,mu1,mu2)*f2(mmu1,mmu2)*Tr(mmu1,nu2,mmu2,nu1)/a^2)*d_(mu1,mu2)*d_(nu1,nu2);

L P27t = 1/4*(d_(mu1,mmu1)*d_(mu2,mmu2)+d_(mu1,mmu2)*d_(mu2,mmu1))*(d_(mmu1,nu1)*d_(mmu2,nu2)+1/a^2*Tr(mmu1,nu2,mmu2,nu1));

contract;
*L a = d_(a1,a3)*d_(a2,a4)+d_(a1,a4)*d_(a2,a3);
*Multiply 1/2;
.sort
*id f1(a1?,b1?)*f1(a2?,b2?) = d_(a1,b1)*d_(a2,b2) + d_(a1,b2)*d_(a2,b1);
repeat id f1(a1?,?a)*f2(b1?,?b) = d_(a1,b1)*f1(?a)*f2(?b);
id f1 = 1;
id f2 = 1;
id P8(mu?,nu?) = d_(mu,nu);

#call SUn
id a = 1/2;
id a^-1 = 2;
*id a =1;
id nf=1;
id NA = Nc^2-1;

.sort
Print +s;
.end

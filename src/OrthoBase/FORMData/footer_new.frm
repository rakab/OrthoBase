repeat id f1(a1?,?a)*f2(b1?,?b) = d_(a1,b1)*f1(?a)*f2(?b);
id f1 = 1;
id f2 = 1;
id [P(8)](mu?,nu?) = d_(mu,nu);

#call SUn
id a = 1/2;
id a^-1 = 2;
id nf=1;
id nf^-1=1;
id Na = Nc^2-1;

.sort
Print;
.end

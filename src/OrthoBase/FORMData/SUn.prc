#procedure SUn
*
*	Procedure to compute color traces for the SU(NF) groups
*	We follow the article by Cvitanovic (Phys.Rev.D14(1976)1536)
*
*	We use [T(i),T(j)] = i_*f(i,j,k)*T(k)    (f is the C of Cvitanovic)
*
*	We use the indices i in the space of the fundamental representation
*	The indices j are in the space of the adjoint.
*	The dimension should be the dimension of the fundamental representation.
*
*	We need: (C)Function Tr(cyclic). Indicates the traces.
*	CFunction T, Tp, f(antisymmetric);
*	Symbols a,nf,NF,NA;
*	Indices i1=NF,i2=NF,i3=NF,i4=NF;
*	Indices j1=NA,j2=NA,j3=NA;
*	Dimension NF;
*
*	Usually the value of a is taken to be 1/2;
*	nf is the number of flavors in the fundamental representation.
*	NF is the dimension of the fundamental representation.
*	NA is the dimension of the adjoint representation.
*
*	Routine by J.Vermaseren, 7-jan-1997
*
repeat;
	id,once,Tr(?a) = SUNT(?a,i1,i1)*nf;
	sum i1;
	repeat;
		id,once,SUNT(j1?,j2?,?a,i1?,i2?) = SUNT(j1,i1,i3)*SUNT(j2,?a,i3,i2);
		sum i3;
	endrepeat;
endrepeat;
#do i = 1,1
if ( count(SUNF,1) || match(SUNT(j1?,i1?,i2?)*SUNT(j1?,i3?,i4?)) )
		redefine i "0";
id,once,SUNF(j1?,j2?,j3?) = 1/a/i_*SUNT(j1,i1,i2)*SUNT(j2,i2,i3)*SUNT(j3,i3,i1)
			-1/a/i_*SUNT(j3,i1,i2)*SUNT(j2,i2,i3)*SUNT(j1,i3,i1);
sum i1,i2,i3;
id	SUNT(j1?,i1?,i2?)*SUNT(j1?,i3?,i4?) = Tp(i1,i2,i3,i4);
#do j = 1,1
if ( count(Tp,1) ) redefine j "0";
.sort
id,once,Tp(i1?,i2?,i3?,i4?) =
			a*(d_(i1,i4)*d_(i2,i3)-d_(i1,i2)*d_(i3,i4)/Nc);
*renumber;
#enddo
#enddo
#do i = 1,1
if ( count(SUND,1) || match(SUNT(j1?,i1?,i2?)*SUNT(j1?,i3?,i4?)) )
		redefine i "0";
id,once,SUND(j1?,j2?,j3?) = 1/a*SUNT(j1,i1,i2)*SUNT(j2,i2,i3)*SUNT(j3,i3,i1)
			+1/a*SUNT(j3,i1,i2)*SUNT(j2,i2,i3)*SUNT(j1,i3,i1);
sum i1,i2,i3;
id	SUNT(j1?,i1?,i2?)*SUNT(j1?,i3?,i4?) = Tp(i1,i2,i3,i4);
#do j = 1,1
if ( count(Tp,1) ) redefine j "0";
.sort
id,once,Tp(i1?,i2?,i3?,i4?) =
			a*(d_(i1,i4)*d_(i2,i3)-d_(i1,i2)*d_(i3,i4)/Nc);
*renumber;
#enddo
#enddo

repeat;
	id	SUNT(j1?,?a,i1?,i2?)*SUNT(j2?,?b,i2?,i3?) = SUNT(j1,?a,j2,?b,i1,i3);
endrepeat;
id	SUNT(?a,i1?,i1?) = Tr(?a)/nf;
id	Tr(j1?) = 0;
id	Tr(j1?,j2?) = a*d_(j1,j2)*nf;
.sort
id	cF = a*(Nc^2-1)/Nc;
id	cA = 2*a*Nc;
id	[cF-cA/6] = a*(2*Nc/3-1/Nc);
.sort
*
#endprocedure

*// vim: set filetype=form:

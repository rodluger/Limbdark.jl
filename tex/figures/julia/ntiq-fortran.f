CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ntiq.f								      C
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C Functions related to the calculations of obscured stellar flux taking       C
C into account quadratic limb darkening. These functions also provide the     C
C parametric partial derivatives of the obscuring function. A bit redundant   C
C coding, however, straight and hopefully clean.                              C
C Note that this particular implementation is an all-in-one stuff, i.e. in    C
C this module, the implementations of the elliptical integrals and a 'main'   C
C program are also provided.						      C
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C (c) 2007, 2008; Pal, A. (apal@szofi.elte.hu)				      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function carlson_elliptic_rc(x,y)
	real*8	carlson_elliptic_rc
	real*8	x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG
	real*8	COMP1,COMP2,THIRD,C1,C2,C3,C4
	parameter (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37)
	parameter (TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0)
	parameter (THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
	real*8	alamb,ave,s,w,xt,yt

	if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+
     *		abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))pause 
     *		'invalid arguments in carlson_elliptic_rc'
	if(y.gt.0.d0)then
		xt=x
		yt=y
		w=1.
	else
		xt=x-y
		yt=-y
		w=sqrt(x)/sqrt(xt)
	endif

1	continue
		alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
		xt=.25d0*(xt+alamb)
		yt=.25d0*(yt+alamb)
		ave=THIRD*(xt+yt+yt)
		s=(yt-ave)/ave
	if(abs(s).gt.ERRTOL)goto 1

	carlson_elliptic_rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)

	return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      function carlson_elliptic_rj(x,y,z,p)
	real*8	carlson_elliptic_rj
	real*8	p,x,y,z,ERRTOL
	real*8	TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
	parameter (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0)
	parameter (C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3)
	parameter (C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)

	real*8	carlson_elliptic_rc,carlson_elliptic_rf
	real*8	a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee
	real*8	fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt

	if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,
     *		z,abs(p)).gt.BIG)pause 'invalid arguments in carlson_elliptic_rj'
	sum=0.d0
	fac=1.d0

	if(p.gt.0.d0)then
		xt=x
		yt=y
		zt=z
	   	pt=p
	else
		xt=min(x,y,z)
		zt=max(x,y,z)
		yt=x+y+z-xt-zt
		a=1.d0/(yt-p)
		b=a*(zt-yt)*(yt-xt)
		pt=yt+b
		rho=xt*zt/yt
		tau=p*pt/yt
		rcx=carlson_elliptic_rc(rho,tau)
	endif

1	continue
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
	        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
	       	alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
	        beta=pt*(pt+alamb)**2
	        sum=sum+fac*carlson_elliptic_rc(alpha,beta)
	        fac=.25d0*fac
	        xt=.25d0*(xt+alamb)
	        yt=.25d0*(yt+alamb)
	        zt=.25d0*(zt+alamb)
	        pt=.25d0*(pt+alamb)
	        ave=.2d0*(xt+yt+zt+pt+pt)
	        delx=(ave-xt)/ave
	        dely=(ave-yt)/ave
	        delz=(ave-zt)/ave
	        delp=(ave-pt)/ave
	if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1

	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.d0*ec
	ee=eb+2.d0*delp*(ea-ec)

	carlson_elliptic_rj=3.d0*sum+
     +		fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *		(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
	if (p.le.0.d0) carlson_elliptic_rj=a*(b*carlson_elliptic_rj+
     +		3.d0*(rcx-carlson_elliptic_rf(xt,yt,zt)))
	return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      function carlson_elliptic_rf(x,y,z)
	real*8	carlson_elliptic_rf
	real*8	x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
	parameter (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0)
	parameter (C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
	real*8	alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

	if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *		z).gt.BIG)pause 'invalid arguments in carlson_elliptic_rf'

	xt=x
	yt=y
	zt=z

1	continue
		sqrtx=sqrt(xt)
	        sqrty=sqrt(yt)
	        sqrtz=sqrt(zt)
	        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
	        xt=.25d0*(xt+alamb)
	        yt=.25d0*(yt+alamb)
	        zt=.25d0*(zt+alamb)
	        ave=THIRD*(xt+yt+zt)
	        delx=(ave-xt)/ave
	        dely=(ave-yt)/ave
	        delz=(ave-zt)/ave
	if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1

	e2=delx*dely-delz**2
	e3=delx*dely*delz

	carlson_elliptic_rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)

	return

      end

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      function carlson_elliptic_rd(x,y,z)
	real*8	carlson_elliptic_rd
	real*8	x,y,z
	parameter (C1=3.0d0/14.0d0,C2=1.0d0/6.0d0,C3=9.0d0/22.0d0)
	parameter (C4=3.0d0/26.0d0,C5=0.25d0*C3,C6=1.5d0*C4)
	real*8	alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac
	real*8	sqrtx,sqrty,sqrtz,sum,xt,yt,zt,ans

	xt=x
	yt=y
	zt=z
	sum=0.0d0
	fac=1.0d0

1	continue
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		sum=sum+fac/(sqrtz*(zt+alamb))
		fac=0.25d0*fac
		xt=0.25d0*(xt+alamb)
		yt=0.25d0*(yt+alamb)
		zt=0.25d0*(zt+alamb)
		ave=0.2d0*(xt+yt+3.0d0*zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
	if(max(abs(delx),abs(dely),abs(delz)).gt.0.0015d0)goto 1

	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0d0*eb	
	ee=ed+ec+ec
	ans=3.0d0*sum+fac*(1.0d0+ed*(-C1+C5*ed-C6*delz*ee)+
     +		delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))

	carlson_elliptic_rd=ans
	return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function ntiq(p,z0,g1,g2)
	implicit none

	real*8	ntiq

	real*8	p,z0,g1,g2,z
	real*8	w,w0,w1,w2,f0,f1,fk,fe,fp,f2,n,k,df

	real*8	a,b,ci,cik,cie,cip,cg,cgk,cge,cgp,k0,k1,t
	real*8	q,ek,ee,ep

	real*8	pi,SMALL

	real*8	carlson_elliptic_rf,carlson_elliptic_rj,carlson_elliptic_rd

	w=6.0d0-2.0d0*g1-g2
	w0=(6.0d0-6.0d0*g1-12.0d0*g2)/w
	w1=(   6.0d0*g1+12.0d0*g2   )/w
	w2=(        6.0d0*g2        )/w

	z=abs(z0)

	pi=2.0d0*acos(0.0d0)

	SMALL=1e-13

C Case 'A'
	if(z.le.0.0d0.and.p.le.1.0d0) then
		f0=p*p
		f1=2.0d0/3.0d0*(1.0d0-(1.0d0-f0)*sqrt(1.0d0-f0))
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		f2=0.5d0*f0*f0
		k=0.0d0
		n=0.0d0
C Case 'A_G'
	elseif(z.le.p-1.0d0) then
		f0=1.0d0
		f1=2.0d0/3.0d0
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		f2=0.5d0
		k=0.0d0
		n=0.0d0
C Case 'B'
	elseif(z.lt.p.and.z.lt.1.0d0-p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2
		f0=p*p
		f1=2.0d0/3.0d0
		ci=2.0d0/(9.0d0*pi*sqrt(1.0d0-a))
		cik=(1.0d0-5.0d0*z*z+f0+a*b)
		cie=(z*z+7.0d0*f0-4.0d0)*(1.0d0-a)
		cip=-3.0d0*(p+z)/(p-z)
		fk=ci*cik
		fe=ci*cie
		fp=ci*cip
		f2=0.5d0*f0*(f0+2.0d0*z*z)
		k=sqrt(4.0d0*z*p/(1.0d0-a))
		n=(a-b)/a
C Case 'B_T'
	elseif(z.lt.p.and.abs(z-1.0d0+p).le.SMALL) then
		f0=p*p
		f1=2.0d0/(3.0d0*pi)*acos(1.0d0-2.0d0*p)-4.0d0/(9.0d0*pi)*(3.0d0+2.0d0*p-8.0d0*f0)*sqrt(p*(1.0d0-p))
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		f2=0.5d0*f0*(f0+2.0d0*z*z)
		k=0.0d0
		n=0.0d0
C Case 'B_G'
	elseif(z.lt.p) then
		a=(p-z)**2
		b=(p+z)**2
		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
		k1=acos((1.0d0+z*z-p*p)/(2.0d0*z))
		f0=(p*p*k0+k1-sqrt(z*z-0.25d0*(1+z*z-p*p)**2))/pi
		f1=2.0d0/3.0d0
		cg=1.0d0/(9.0d0*pi*sqrt(p*z))
		cgk=((1.0d0-b)*(2.0d0*b+a-3.0d0)-3.0d0*(p+z)*(p-z)*(b-2.0d0))
		cge=4.0d0*p*z*(z*z+7*p*p-4.0d0)
		cgp=-3.0d0*(p+z)/(p-z)
		fk=cg*cgk
		fe=cg*cge
		fp=cg*cgp
		f2=(k1+p*p*(p*p+2.0d0*z*z)*k0-0.25d0*(1.0d0+5.0d0*p*p+z*z)*sqrt((1.0d0-a)*(b-1.0d0)))/(2.0d0*pi)
		k=sqrt((1.0d0-a)/(4.0d0*z*p))
		n=(a-1.0d0)/a
C Case 'C'
	elseif(abs(z-p).le.SMALL.and.z.lt.1.0d0-p-SMALL) then
		f0=p*p
		f1=1.0d0/3.0d0
		t=2.0d0/(9.0d0*pi)
		fk=t*(1.0d0-4.0d0*f0)
		fe=4.0d0*t*(2.0d0*f0-1.0d0)
		fp=0.0d0
		f2=1.5d0*f0*f0
		k=2.0d0*p
		n=0.0d0
C Case 'C_T'
	elseif(abs(z-p).le.SMALL.and.abs(z-1.0d0+p).le.SMALL) then
		f0=0.25d0
		f1=1.0d0/3.0d0-4.0d0/(9.0d0*pi)
		f2=3.0d0/32.0d0
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		n=0.0d0
		k=0.0d0
C Case 'C_G'
	elseif(abs(z-p).le.SMALL) then
		a=(p-z)**2
		b=(p+z)**2
		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
		k1=acos((1.0d0+z*z-p*p)/(2.0d0*z))
		f0=(p*p*k0+k1-sqrt(z*z-0.25d0*(1.0d0+z*z-p*p)**2))/pi
		f1=1.0d0/3.0d0
		fk=-(1.0d0-4.0d0*p*p)*(3.0d0-8.0d0*p*p)/(9.0d0*pi*p)
		fe=16.0d0*p*(2.0d0*p*p-1.0d0)/(9.0d0*pi)
		fp=0.0d0
		f2=(k1+p*p*(p*p+2.0d0*z*z)*k0-0.25d0*(1.0d0+5.0d0*p*p+z*z)*sqrt((1.0d0-a)*(b-1.0d0)))/(2.0d0*pi)
		k=1.0d0/(2.0d0*p)
		n=0.0d0
C Case 'D'
	elseif(z.lt.1.0d0-p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2
		f0=p*p
		f1=0.0d0
		ci=2.0d0/(9.0d0*pi*sqrt(1.0d0-a))
		cik=(1.0d0-5.0d0*z*z+f0+a*b)
		cie=(z*z+7.0d0*f0-4.0d0)*(1.0d0-a)
		cip=-3.0d0*(p+z)/(p-z)
		fk=ci*cik
		fe=ci*cie
		fp=ci*cip
		f2=0.5d0*f0*(f0+2.0d0*z*z)
		k=sqrt(4.0d0*z*p/(1.0d0-a))
		n=(a-b)/a
C Case 'E'
	elseif(abs(z-1.0d0+p).le.SMALL) then
		f0=p*p
		f1=2.0d0/(3.0d0*pi)*acos(1.0d0-2.0d0*p)-4.0d0/(9.0d0*pi)*(3.0d0+2.0d0*p-8.0d0*f0)*sqrt(p*(1.0d0-p))
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		f2=0.5d0*f0*(f0+2.0d0*z*z)
		k=0.0d0
		n=0.0d0
C Case 'F'
	elseif(z.lt.1.0d0+p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2
		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
		k1=acos((1.0d0+z*z-p*p)/(2.0d0*z))
		f0=(p*p*k0+k1-sqrt(z*z-0.25d0*(1.0d0+z*z-p*p)**2))/pi
		f1=0.0d0
		cg=1.0d0/(9.0d0*pi*sqrt(p*z))
		cgk=((1.0d0-b)*(2.0d0*b+a-3.0d0)-3.0d0*(p+z)*(p-z)*(b-2.0d0))
		cge=4.0d0*p*z*(z*z+7.0d0*p*p-4.0d0)
		cgp=-3.0d0*(p+z)/(p-z)
		fk=cg*cgk
		fe=cg*cge
		fp=cg*cgp
		f2=(k1+p*p*(p*p+2.0d0*z*z)*k0-0.25d0*(1.0d0+5.0d0*p*p+z*z)*sqrt((1.0d0-a)*(b-1.0d0)))/(2.0d0*pi)
		k=sqrt((1.0d0-a)/(4.0d0*z*p))
		n=(a-1.0d0)/a
	else
		f0=0.0d0
		f1=0.0d0
		fk=0.0d0
		fe=0.0d0
		fp=0.0d0
		f2=0.0d0
		k=0.0d0
		n=0.0d0
	endif

	df=w0*f0+w1*f1+w2*f2
	if(fk.ne.0.0d0.or.fe.ne.0.0d0) then
		q=(1.0d0-k)*(1.0d0+k)
		ek=carlson_elliptic_rf(0.0d0,q,1.0d0)
		df=df+w1*fk*ek
		ee=ek-k*k*carlson_elliptic_rd(0.0d0,q,1.0d0)/3.0d0
		df=df+w1*fe*ee
		if(fp.ne.0.0d0) then
			ep=ek+n*carlson_elliptic_rj(0.0d0,q,1.0d0,1.0d0-n)/3.0d0
			df=df+w1*fp*ep
		endif
	endif

	ntiq=df

	return
      end

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      function ntiq_diff(p,z0,g1,g2,diff)
	implicit none

	real*8	ntiq_diff

	real*8	p,z0,g1,g2,z,diff(4)
	real*8	w,w0,w1,w2,w0g1,w0g2,w1g1,w1g2,w2g1,w2g2
	real*8	f0,f1,fk,fe,fp,f2,n,k,df
	real*8	f0pp,f0zz,f1pp,f1zz,fkpp,fkzz,fepp,fezz,fppp,fpzz,f2pp,f2zz	
	real*8	npp,nzz,kpp,kzz
	integer	calc_eke,calc_epi

	real*8	a,b,ci,cik,cie,cip,cg,cgk,cge,cgp,k0,k1,t,t1,t2
	real*8	q,ek,ee,ep,signz

	real*8	pi,SMALL

	real*8	carlson_elliptic_rf,carlson_elliptic_rj,carlson_elliptic_rd

	w=6.0d0-2.0d0*g1-g2
	w0=(6.0d0-6.0d0*g1-12.0d0*g2)/w
	w0g1=(2*w0-6)/w
	w0g2=(w0-12)/w
	w1=(   6.0d0*g1+12.0d0*g2   )/w
	w1g1=(6+2*w1)/w
	w1g2=(12+w1)/w
	w2=(        6.0d0*g2        )/w
	w2g1=2*w2/w
	w2g2=(6+w2)/w

	if(z0.lt.0.0d0)then
		signz=-1.0d0
		z=-z0
	else
		signz=1.0d0
		z=z0
	endif

	pi=2.0d0*acos(0.0d0)

	SMALL=1e-13

C Case 'A'
	if(z.le.0.0d0.and.p.le.1.0d0) then
		f0=p*p
		f0pp=2*p
		f0zz=0.0d0

		f1=2.0d0/3.0d0*(1.0d0-(1.0d0-f0)*sqrt(1.0d0-f0))
		f1pp=2*p*sqrt(1.0d0-f0)
		f1zz=0.0d0

		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=0.5d0*f0*f0
		f2pp=2*f0*p
		f2zz=0.0d0

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=0
		calc_epi=0
C Case 'A_G'
	elseif(z.le.p-1.0d0) then
		f0=1.0d0
		f0pp=0.0d0
		f0zz=0.0d0

		f1=2.0d0/3.0d0
		f1pp=0.0d0
		f1zz=0.0d0

		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=0.5d0
		f2pp=0.0d0
		f2zz=0.0d0

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=0
		calc_epi=0
C Case 'B'
	elseif(z.lt.p.and.z.lt.1.0d0-p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2

		f0=p*p
		f0pp=2*p
		f0zz=0.0d0

		f1=2.0d0/3.0d0
		f1pp=0.0d0
		f1zz=0.0d0

		ci=2.0d0/(9.0d0*pi*sqrt(1.0d0-a))
		cik=(1.0d0-5.0d0*z*z+f0+a*b)
		cie=(z*z+7.0d0*f0-4.0d0)*(1.0d0-a)
		cip=-3.0d0*(p+z)/(p-z)
		fk=ci*cik
		fkpp=+ci*2*p*(1+2*(p*p-z*z))+ci*cik*(p-z)/(1-a)
		fkzz=-ci*2*z*(5+2*(p*p-z*z))-ci*cik*(p-z)/(1-a)
		fe=ci*cie
		fepp=+ci*(14*p*(1-a)-2*(p-z)*(z*z+7*f0-4))+ci/(1-a)*(p-z)*cie
		fezz=+ci*( 2*z*(1-a)+2*(p-z)*(z*z+7*f0-4))-ci/(1-a)*(p-z)*cie
		fp=ci*cip
		fppp=+ci*3.0*(+(2*z)/a-(p+z)/(1-a))
		fpzz=+ci*3.0*(-(2*p)/a+(p+z)/(1-a))

		f2=0.5d0*f0*(f0+2.0d0*z*z)
		f2pp=p*(f0+2*z*z)+f0*p
	        f2zz=2*f0*z

		k=sqrt(4.0d0*z*p/(1.0d0-a))
	        kpp=2*z*(1+p*p-z*z)/((1-a)*(1-a)*k)
	        kzz=2*p*(1-p*p+z*z)/((1-a)*(1-a)*k)
		n=(a-b)/a
	        npp=+4*z*(p+z)/((p-z)*a)
	        nzz=-4*p*(p+z)/((p-z)*a)

	        calc_eke=1
	        calc_epi=1
C Case 'B_T'
	elseif(z.lt.p.and.abs(z-1.0d0+p).le.SMALL) then
		f0=p*p
	        f0pp=2*p
	        f0zz=0.0

		f1=2.0d0/(3.0d0*pi)*acos(1.0d0-2.0d0*p)-4.0d0/(9.0d0*pi)*(3.0d0+2.0d0*p-8.0d0*f0)*sqrt(p*(1.0d0-p))
	        f1pp=(8.0d0/pi)*p*sqrt(p*(1-p))
	        f1zz=-f1pp/3.0d0

		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=0.5d0*f0*(f0+2.0d0*z*z)
	        f2pp=p*(f0+2*z*z)+f0*p
	        f2zz=2*f0*z

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=0
		calc_epi=0
C Case 'B_G'
	elseif(z.lt.p) then
		a=(p-z)**2
		b=(p+z)**2
		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
		k1=acos((1.0d0+z*z-p*p)/(2.0d0*z))

	        t=sqrt(z*z-0.25*(1+z*z-p*p)**2)
	        f0=(p*p*k0+k1-t)/pi
	        f0pp=(2*p*k0)/pi
	        f0zz=(-2*p*sin(k0))/pi

		f1=2.0d0/3.0d0
	        f1pp=0.0
	        f1zz=0.0

		cg=1.0d0/(9.0d0*pi*sqrt(p*z))
		cgk=((1.0d0-b)*(2.0d0*b+a-3.0d0)-3.0d0*(p+z)*(p-z)*(b-2.0d0))
		cge=4.0d0*p*z*(z*z+7*p*p-4.0d0)
		cgp=-3.0d0*(p+z)/(p-z)

		fk=cg*cgk
	        fkpp=-fk/(2*p)-cg*2*(p*p*(12*p+21*z)+z*(z*z-4)+2*p*(5*z*z-6))
	        fkzz=-fk/(2*z)-cg*2*p*(-4+7*p*p+10*p*z+3*z*z)
		fe=cg*cge
	        fepp=-fe/(2*p)+cg*4*z*(-4+21*p*p+  z*z)
	        fezz=-fe/(2*z)+cg*4*p*(-4+ 7*p*p+3*z*z)
		fp=cg*cgp
	        fppp=-fp/(2*p)+cg*6*z/a
	        fpzz=-fp/(2*z)-cg*6*p/a

	        t=sqrt((1-a)*(b-1))
	        f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*t)/(2.0*pi)
	        f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/pi
	        f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/pi

	        k=sqrt((1-a)/(4*p*z))
	        kpp=-(1+p*p-z*z)/(8*k*p*p*z)
	        kzz=-(1-p*p+z*z)/(8*k*p*z*z)
	        n=(a-1)/a
	        npp=2/(a*(p-z))
	        nzz=-npp

	        calc_eke=1
	        calc_epi=1
C Case 'C'
	elseif(abs(z-p).le.SMALL.and.z.lt.1.0d0-p-SMALL) then
		f0=p*p
		f0pp=2*p
		f0zz=0.0d0

		f1=1.0d0/3.0d0
		f1pp=0.0d0
		f1zz=0.0d0

		t=2.0d0/(9.0d0*pi)
		fk=t*(1.0d0-4.0d0*f0)
	        fkpp=+2.0/(9.0*pi)* 2*p-1.0/(3*pi*p)
	        fkzz=-2.0/(9.0*pi)*10*p+1.0/(3*pi*p)
		fe=4.0d0*t*(2.0d0*f0-1.0d0)
	        fepp=+t*(14*p)+1.0/(3*pi*p)
	        fezz=+t*( 2*p)-1.0/(3*pi*p)

		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=1.5d0*f0*f0
	        f2pp=4*f0*p
	        f2zz=2*f0*p

		k=2.0d0*p
	        kpp=1.0d0
	        kzz=1.0d0
		n=0.0d0
	        npp=0.0d0
	        nzz=0.0d0

	        calc_eke=1
	        calc_epi=0
C Case 'C_T'
	elseif(abs(z-p).le.SMALL.and.abs(z-1.0d0+p).le.SMALL) then
		f0=0.25d0
	        f0pp=1.0d0
	        f0zz=0.0d0

		f1=1.0d0/3.0d0-4.0d0/(9.0d0*pi)
	        f1pp=+2.0/pi
	        f1zz=-2.0/(3*pi)

		f2=3.0d0/32.0d0
		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0

		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

	        calc_eke=0
	        calc_epi=0
C Case 'C_G'
	elseif(abs(z-p).le.SMALL) then
		a=(p-z)**2
		b=(p+z)**2

		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
	        t=(1+z*z-p*p)/(2*z)
		k1=acos(t)
		t=sqrt(z*z-0.25*(1+z*z-p*p)*(1+z*z-p*p))

	        f0=(p*p*k0+k1-t)/pi
	        f0pp=(2*p*k0)/pi
	        f0zz=(-2*p*sin(k0))/pi

		f1=1.0d0/3.0d0
		f1pp=0.0d0
		f1zz=0.0d0

	        cg=1.0/(9.0*pi*p)
	        cgk=(1-4*p*p)*(8*p*p-3)
	        cge=16*p*p*(2*p*p-1)
		fk=-(1.0d0-4.0d0*p*p)*(3.0d0-8.0d0*p*p)/(9.0d0*pi*p)
	        fkpp=-fk/(2*p)-cg*2*p*(44*p*p-16)-2.0/(3*pi)
	        fkzz=-fk/(2*p)-cg*2*p*(20*p*p- 4)+2.0/(3*pi)
		fe=16.0d0*p*(2.0d0*p*p-1.0d0)/(9.0d0*pi)
	        fepp=-fe/(2*p)+cg*4*p*(-4+22*p*p)+2.0/(3*pi)
	        fezz=-fe/(2*p)+cg*4*p*(-4+10*p*p)-2.0/(3*pi)

		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

	        t=sqrt(4*p*p-1)
	        f2=(k1+p*p*(3*p*p)*k0-0.25*(1+6*p*p)*t)/(2.0*pi)
	        f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/pi
	        f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/pi

		k=1.0d0/(2.0d0*p)
	        kpp=-1/(4*p*p)
	        kzz=-1/(4*p*p)
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=1
		calc_epi=0
C Case 'D'
	elseif(z.lt.1.0d0-p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2

		f0=p*p
		f0pp=2*p
		f0zz=0.0d0

		f1=0.0d0
		f1pp=0.0d0
		f1zz=0.0d0

		ci=2.0d0/(9.0d0*pi*sqrt(1.0d0-a))
		cik=(1.0d0-5.0d0*z*z+f0+a*b)
		cie=(z*z+7.0d0*f0-4.0d0)*(1.0d0-a)
		cip=-3.0d0*(p+z)/(p-z)
		fk=ci*cik
		fkpp=+ci*2*p*(1+2*(p*p-z*z))+ci*cik*(p-z)/(1-a)
		fkzz=-ci*2*z*(5+2*(p*p-z*z))-ci*cik*(p-z)/(1-a)
		fe=ci*cie
		fepp=+ci*(14*p*(1-a)-2*(p-z)*(z*z+7*f0-4))+ci/(1-a)*(p-z)*cie
		fezz=+ci*( 2*z*(1-a)+2*(p-z)*(z*z+7*f0-4))-ci/(1-a)*(p-z)*cie
		fp=ci*cip
		fppp=+ci*3.0*(+(2*z)/a-(p+z)/(1-a))
		fpzz=+ci*3.0*(-(2*p)/a+(p+z)/(1-a))

		f2=0.5d0*f0*(f0+2.0d0*z*z)
	        f2pp=p*(f0+2*z*z)+f0*p
	        f2zz=2*f0*z

		k=sqrt(4.0d0*z*p/(1.0d0-a))
	        kpp=2*z*(1+p*p-z*z)/((1-a)*(1-a)*k)
	        kzz=2*p*(1-p*p+z*z)/((1-a)*(1-a)*k)
		n=(a-b)/a
	        npp=+4*z*(p+z)/((p-z)*a)
	        nzz=-4*p*(p+z)/((p-z)*a)

	        calc_eke=1
	        calc_epi=1

C Case 'E'
	elseif(abs(z-1.0d0+p).le.SMALL) then
		f0=p*p
		f0pp=2*p
		f0zz=0.0d0

		f1=2.0d0/(3.0d0*pi)*acos(1.0d0-2.0d0*p)-4.0d0/(9.0d0*pi)*(3.0d0+2.0d0*p-8.0d0*f0)*sqrt(p*(1.0d0-p))
	        f1pp=(8.0/pi)*p*sqrt(p*(1-p))
	        f1zz=-f1pp/3.0
		
		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=0.5d0*f0*(f0+2.0d0*z*z)
		f2pp=p*(f0+2*z*z)+f0*p
		f2zz=2*f0*z

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=0
		calc_epi=0
C Case 'F'
	elseif(z.lt.1.0d0+p-SMALL) then
		a=(p-z)**2
		b=(p+z)**2
		k0=acos((p*p+z*z-1.0d0)/(2.0d0*p*z))
		k1=acos((1.0d0+z*z-p*p)/(2.0d0*z))

	        t=sqrt(z*z-0.25*(1+z*z-p*p)**2)
	        f0=(p*p*k0+k1-t)/pi
	        f0pp=(2*p*k0)/pi
	        f0zz=(-2*p*sin(k0))/pi

		f1=0.0d0
	        f1pp=0.0d0
	        f1zz=0.0d0

		cg=1.0d0/(9.0d0*pi*sqrt(p*z))
		cgk=((1.0d0-b)*(2.0d0*b+a-3.0d0)-3.0d0*(p+z)*(p-z)*(b-2.0d0))
		cge=4.0d0*p*z*(z*z+7*p*p-4.0d0)
		cgp=-3.0d0*(p+z)/(p-z)

		fk=cg*cgk
	        fkpp=-fk/(2*p)-cg*2*(p*p*(12*p+21*z)+z*(z*z-4)+2*p*(5*z*z-6))
	        fkzz=-fk/(2*z)-cg*2*p*(-4+7*p*p+10*p*z+3*z*z)
		fe=cg*cge
	        fepp=-fe/(2*p)+cg*4*z*(-4+21*p*p+  z*z)
	        fezz=-fe/(2*z)+cg*4*p*(-4+ 7*p*p+3*z*z)
		fp=cg*cgp
	        fppp=-fp/(2*p)+cg*6*z/a
	        fpzz=-fp/(2*z)-cg*6*p/a

	        t=sqrt((1-a)*(b-1))
	        f2=(k1+p*p*(p*p+2*z*z)*k0-0.25*(1+5*p*p+z*z)*t)/(2.0*pi)
	        f2pp=(2*p*(p*p+z*z)*k0-4*z*p*p*sin(k0))/pi
	        f2zz=(-2*p*(p*p+z*z)*sin(k0)+z*p*p*(2*k0+sin(2*k0)))/pi

	        k=sqrt((1-a)/(4*p*z))
	        kpp=-(1+p*p-z*z)/(8*k*p*p*z)
	        kzz=-(1-p*p+z*z)/(8*k*p*z*z)
	        n=(a-1)/a
	        npp=2/(a*(p-z))
	        nzz=-npp

	        calc_eke=1
	        calc_epi=1
C Case 'G'
	else
		f0=0.0d0
		f0pp=0.0d0
		f0zz=0.0d0

		f1=0.0d0
		f1pp=0.0d0
		f1zz=0.0d0

		fk=0.0d0
		fkpp=0.0d0
		fkzz=0.0d0
		fe=0.0d0
		fepp=0.0d0
		fezz=0.0d0
		fp=0.0d0
		fppp=0.0d0
		fpzz=0.0d0

		f2=0.0d0
		f2pp=0.0d0
		f2zz=0.0d0

		k=0.0d0
		kpp=0.0d0
		kzz=0.0d0
		n=0.0d0
		npp=0.0d0
		nzz=0.0d0

		calc_eke=0
		calc_epi=0

	endif

	df=w0*f0+w1*f1+w2*f2
	diff(1)=w0*f0pp+w1*f1pp+w2*f2pp
	diff(2)=w0*f0zz+w1*f1zz+w2*f2zz
	diff(3)=w0g1*f0+w1g1*f1+w2g1*f2
	diff(4)=w0g2*f0+w1g2*f1+w2g2*f2

	if(calc_eke.ne.0) then
		q=(1.0d0-k)*(1.0d0+k)
		ek=carlson_elliptic_rf(0.0d0,q,1.0d0)
		df=df+w1*fk*ek
		if(fp.ne.0) then
			t=fp/(2*n*(n-1))
		else
			t=0
		endif
	        diff(1)=diff(1)+w1*ek*(fkpp-(fk+fe)*kpp/k+t*npp)
	        diff(2)=diff(2)+w1*ek*(fkzz-(fk+fe)*kzz/k+t*nzz)
	        diff(3)=diff(3)+w1g1*fk*ek
	        diff(4)=diff(4)+w1g2*fk*ek
		ee=ek-k*k*carlson_elliptic_rd(0.0d0,q,1.0d0)/3.0d0
		df=df+w1*fe*ee
		if(fp.ne.0) then
			t1=fp*k/((n-k*k)*(k*k-1))
			t2=fp/(2*(k*k-n)*(n-1))
		else
			t1=0.0d0
			t2=0.0d0
		endif
	        diff(1)=diff(1)+w1*ee*(fepp+(fk/(k*(1-k*k))+fe/k+t1)*kpp+t2*npp)
	        diff(2)=diff(2)+w1*ee*(fezz+(fk/(k*(1-k*k))+fe/k+t1)*kzz+t2*nzz)
	        diff(3)=diff(3)+w1g1*fe*ee
	        diff(4)=diff(4)+w1g2*fe*ee
		if(calc_epi.ne.0) then
			ep=ek+n*carlson_elliptic_rj(0.0d0,q,1.0d0,1.0d0-n)/3.0d0
			df=df+w1*fp*ep
	                diff(3)=diff(3)+w1g1*fp*ep
	                diff(4)=diff(4)+w1g2*fp*ep
		endif
	endif

	diff(2)=diff(2)*signz

	ntiq_diff=df

	return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program ntiqtest
	integer i
	real*8	ntiq_diff,n,p,z,diff(4)
        p = 0.1d0
        do 10 i=0,1000
          z = sqrt(((1d0+2d0*p)*(i-500.0)/500.)**2)
     	  n=ntiq_diff(p,z,0.2d0,0.3d0,diff)
    	  write(6,*) n,diff(1),diff(2),diff(3),diff(4)
10	continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	     

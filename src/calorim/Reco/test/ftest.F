
        subroutine ftest(nnn)        

        double precision tutu(3,3),b(3)
	integer		ir(2000)

c       0.*x + 1.*y + 2.*z = 3.
c      -1.*x - 1.*y + 1.*z =-1.
c       2.*x + 1.*y + 3.*z = 6.
c 
         tutu(1,1)= 0.
         tutu(2,1)=-1.
         tutu(3,1)= 2.

         tutu(1,2)= 1.
         tutu(2,2)=-1.
         tutu(3,2)= 1.

         tutu(1,3)= 2.
         tutu(2,3)= 1.
         tutu(3,3)= 3.
         
         b(1) = 3.
         b(2) =-1.
         b(3) = 6.
 
         nf=3
         neq=1
         ifail=0
 
         call deqinv(nf,tutu,nf,ir,ifail,neq,b)
 
         write (*,*) ' ifail  ',ifail
         write (*,*) ' x vect ',b(1),b(2),b(3)
        
        return
        end
        
	subroutine bsig2b(np,nf,b,sig2,wb)
c
	common /ifail/	ifail
	integer		ir(2000)
	dimension	wb(*)
	dimension	b(*)
	dimension	sig2(np)
c
	do i=1,nf
	do j=i,nf
		s=0.
		do k=1,np
		ki=k+np*(i-1)
		kj=k+np*(j-1)
			s=s+b(ki)*b(kj)*sig2(k)
		enddo
		ij=i+nf*(j-1)
		ji=j+nf*(i-1)
		wb(ij)=s
		if(i.ne.j) wb(ji)=s
	enddo
	enddo
c	
	if(nf.eq.1) then
		wb(1)=1./wb(1)
		ifail = 0
	else
		call rinv(nf,wb,nf,ir,ifail)
	endif
c	
	return
	end
c

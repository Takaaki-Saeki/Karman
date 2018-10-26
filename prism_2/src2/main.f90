!**********************
!*角柱周りの二次元剥離流れ*
!**********************
!
!
!---------------
  module cmpcnd
!---------------
!
     implicit none
!
	 double precision::omegav,errorv,dt
	 double precision::re=70.0,cfl=0.2,omegap=1.0,errorp=1.0e-04
	 integer,parameter::nlast=5000,nlp=250,maxitp=100
!  
  end module cmpcnd
!
!----------------
  module grid2d
!----------------
!
     implicit none
!
     integer::mx,i1,i2,my,j1,j2,j3,j4
	 integer,parameter::mdx=500,mdy=300
     double precision::dx,dy
     double precision,dimension(-2:mdx,-2:mdy)::x,y	 
!
  end module grid2d
!
!-----------------
  module others
!-----------------
!
     implicit none
!
     integer::i,j,icent,jcent,jcent1,jcent2,n,nstep,nbegin,itr,itrp
	 double precision::time,resp,cp1top,cp1btm,cp1fore,cp1back,cp2top,cp2btm,cp2fore,cp2back
     double precision::cp1,cp2,cp3,cp4,cd1,cl1,cd2,cl2,dp,res,ux,uy,vx,vy
!
  end module others
!
!-----------------
  module flowvp
!-----------------
!
     implicit none
!
	 integer,parameter::mdx=500,mdy=300
	 double precision,dimension(-1:mdx,-1:mdy)::u,v,p
!
  end module flowvp
!
!---------------
  module worksp
!---------------
!
     implicit none
!
     integer,parameter::mdx=500,mdy=300
     double precision,dimension(-1:mdx,-1:mdy)::rhs,urhs,vrhs  
!
  end module worksp
!
!================
  program main
!================
!
     use cmpcnd
     implicit none
!
     open(unit=60,file='flow1.data',status='replace')
	 open(unit=50,file='flow2.data',status='replace')
	 open(unit=40,file='flow3.data',status='replace')
	 open(unit=70,file='hist.data',form='formatted',status='replace')
! 計算条件はmodule部分ですでに設定
!
! 計算格子の設定
     call setgrd
! 流れ場を解く
     call slvflw
!
     close(60)
	 close(50)
	 close(40)
     close(70)
  end program main
!  
!====================
  subroutine setflow
!====================
!
! 計算条件の設定
!
     use cmpcnd
!
     implicit none
!
!
     errorp=1.0e-04   ! SOR parameters for p(圧力計算に必要)
!
  end subroutine setflow
!
!===================
  subroutine setgrd
!===================
!
! 計算格子の設定
!
     use grid2d
	 use others
!
     implicit none
!
! set x-grid parameters
     mx=401
     i1=96
     i2=106
     dx=1.0/float(i2-i1)
! set y-grid parameters
     my=201
     j1=67
     j2=77
     j3=124
     j4=134
     dy=1.0/float(j2-j1)
! set xy-grid parameters
     icent=(i1+i2)/2
     jcent=(96+106)/2
     jcent1=(j1+j2)/2
     jcent2=(j3+j4)/2
     do i=1,mx
        do j=1,my
           x(i,j)=dx*float(i-icent)
           y(i,j)=dy*float(j-jcent)
        end do
     end do
!
  end subroutine setgrd
!
!===================
  subroutine slvflw
!===================
!
! 流れを解く
!
     use cmpcnd
     use grid2d	 
     use flowvp
	 use others
!
     implicit none
!
!
! set time step size
     dt=cfl*amin1(dx,dy)
! set initial conditions
     call intcnd
	 call bcforp
	 call bcforv
! time marching
     write(6,*) ' step / res(p) at itr. / CD / CL / Cp1 / Cp2'
	 write(70,*) '    step      res        CD      CL   Cp1   Cp2'
! nステップからn+1ステップを計算することをnlast回繰り返す
     do n=1,nlast
	    nstep=n+nbegin
		time=time+dt
! solve poisson for parameter
        call poiseq             ! 圧力を更新
		call bcforp             ! 圧力の境界条件を設定
! update u,v
        call veloeq             ! 速度の更新
		call bcforv             ! 速度の境界条件の設定
! caluculate CD and CL
        cd1=0.0
        cd2=0.0
		do j=j1,j2-1
		   cp1fore=(2.*p(i1,j)+2.*p(i1,j+1))/2.     ! Cpを積分
		   cp1back=(2.*p(i2,j)+2.*p(i2,j+1))/2.     ! Cp=2pである
		   cd1=cd1+(cp1fore-cp1back)*dy
        end do
        do j=j3, j4-1
           cp2fore=(2.*p(i1,j)+2.*p(i1,j+1))/2.
           cp2back=(2.*p(i2,j)+2.*p(i2,j+1))/2.
           cd2=cd2+(cp2fore-cp2back)*dy
        end do
        cl1=0.0
        cl2=0.0
		do i=i1,i2-1
		   cp1btm=(2.*p(i,j1)+2.*p(i+1,j1))/2.
		   cp1top=(2.*p(i,j2)+2.*p(i+1,j2))/2.
           cp2btm=(2.*p(i,j3)+2.*p(i+1,j3))/2.
           cp2top=(2.*p(i,j4)+2.*p(i+1,j4))/2.
		   cl1=cl1+(cp1btm-cp1top)*dx
           cl2=cl2+(cp2btm-cp2top)*dx
		end do
! 結果をnlpステップおきにモニター
        if((n/nlp)*nlp==n) then
		   cp1=2.0*p(i2+i2-i1,j1)
		   cp2=2.0*p(i2+i2-i1,j2)
           cp3=2.0*p(i2+i2-i1,j3)
           cp4=2.0*p(i2+i2-i1,j4)
		   write(6,601) nstep,resp
		   write(70,602) nstep,resp
		   do i=1,mx
		      do j=1,my
                 write(40,*) p(i,j)
			     write(60,*) u(i,j), v(i,j)
			  end do
		   end do
        end if
     end do
! write final results
     do i=1,mx
	    do j=1,my
		   p(i,j)=2.0*p(i,j)
		end do
		
	 end do
!	 write(60,*) re,cfl,dt,nlast,time
!	 write(60,*) mx,i1,i2,my,j1,j2
	 do i=1,mx
	    do j=1,my
	       write(50,*) x(i,j),y(i,j)
!	       write(60,*) u(i,j),v(i,j)
!		   write(40,*) p(i,j)
		end do
	 end do	
!
601  format(i6,e12.3,i6,4f12.3)
602  format(i10,e15.4,4f15.4)
!
  end subroutine slvflw
!
!===================
  subroutine intcnd
!===================
!
     use cmpcnd
	 use grid2d
	 use flowvp
	 use others
!
     implicit none
!
!
     nbegin=0
	 time=0.0
!
     do i=1,mx
	    do j=1,my
		   u(i,j)=1.0
		   v(i,j)=0.0     !一様条件を至るところ(壁面上は除く)で与える
		   p(i,j)=0.0
		end do
     end do
!
  end subroutine intcnd
!
!=====================
  subroutine bcforp
!=====================
!
! 圧力の境界条件の設定 
     use cmpcnd
	 use grid2d
	 use flowvp
	 use others
!
     implicit none
!
! (1)inflow condition(i=1)
	 do j=1,my
	    p(1,j)=0.0
	 end do
! (2)downstream consition(i=mx)
     do j=1,my
	    p(mx,j)=0.0
	 end do
! (3)bottom condition(j=1)
     do i=1,mx
	    p(i,1)=0.0
	 end do
! (4)bottom condition(j=my)
     do i=1,mx
	    p(i,my)=0.0
	 end do
! (5)wall condition
     p(i1,j1)=p(i1-1,j1-1)
	 p(i1,j2)=p(i1-1,j2+1)
	 p(i2,j1)=p(i2+1,j1-1)
	 p(i2,j2)=p(i2+1,j2+1)
     p(i1,j4)=p(i1-1,j4+1)
     p(i1,j3)=p(i1-1,j3-1)
     p(i2,j4)=p(i2+1,j4+1)
     p(i2,j3)=p(i2+1,j3-1)
	 do j=j1+1,j2-1
	    p(i1,j)=p(i1-1,j)
	 end do
	 do j=j1+1,j2-1
	    p(i2,j)=p(i2+1,j)
	 end do
     do j=j3+1,j4-1
        p(i1,j)=p(i1-1,j)
     end do
     do j=j3+1,j4-1
        p(i2,j)=p(i2+1,j)
     end do
	 do i=i1+1,i2-1
	    p(i,j1)=p(i,j1-1)
	 end do
	 do i=i1+1,i2-1
	    p(i,j2)=p(i,j2+1)
	 end do
     do i=i1+1,i2-1
        p(i,j3)=p(i,j3-1)
     end do
     do i=i1+1,i2-1
        p(i,j4)=p(i,j4+1)
     end do
!
  end subroutine bcforp
!
!
!====================
  subroutine bcforv
!====================
!
! 速度の境界条件の設定
!
     use cmpcnd
	 use grid2d
	 use flowvp
	 use others
!
     implicit none
!
! (1)inflow condition(i=1)
     do j=1,my
	    u(1,j)=1.0
		v(1,j)=0.0
		u(0,j)=1.0
		v(0,j)=0.0
	 end do
! (2)downstream condition(i=mx)
     do j=1,my
	    u(mx,j)=2.0*u(mx-1,j)-u(mx-2,j)
		v(mx,j)=2.0*v(mx-1,j)-v(mx-2,j)
		u(mx+1,j)=2.0*u(mx,j)-u(mx-1,j)
		v(mx+1,j)=2.0*v(mx,j)-v(mx-1,j)
     end do
! (3)bottom condition(j=1)
	 do i=1,mx
	    u(i,1)=2.0*u(i,2)-u(i,3)
		v(i,1)=2.0*v(i,2)-v(i,3)
		u(i,0)=2.0*u(i,1)-u(i,2)
		v(i,0)=2.0*v(i,1)-v(i,2)
	 end do
! (4)bottom condition(j=my)
     do i=1,mx
	    u(i,my)=2.0*u(i,my-1)-u(i,my-2)
		v(i,my)=2.0*v(i,my-1)-v(i,my-2)
		u(i,my+1)=2.0*u(i,my)-u(i,my-1)
		v(i,my+1)=2.0*v(i,my)-v(i,my-1)
     end do
! (5)wall condition
     do i=i1,i2
	    do j=j1,j2
		   u(i,j)=0.0
		   v(i,j)=0.0
		end do
     end do
     do i=i1,i2
        do j=j3,j4
           u(i,j)=0.0
           v(i,j)=0.0
        end do
     end do
!
  end subroutine bcforv
!
!===============================
  subroutine poiseq
!===============================
!
! 緩和法によるポアソン方程式の解法
!
     use cmpcnd
     use grid2d
     use flowvp
     use worksp
	 use others
!
     implicit none
!
! (1)compute RHS(右辺の計算)
     do i=2,mx-1
	    do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then  ! 正方形内部は計算しない
		      ux=(u(i+1,j)-u(i-1,j))/(2.0*dx)
		      uy=(u(i,j+1)-u(i,j-1))/(2.0*dy)
		      vx=(v(i+1,j)-v(i-1,j))/(2.0*dx)
		      vy=(v(i,j+1)-v(i,j-1))/(2.0*dy)
		      rhs(i,j)=(ux+vy)/dt-(ux**2+2.0*uy*vx+vy**2)
		   end if
		end do
     end do
!	 
! (2)iterations
     do itr=1,maxitp
	    res=0.0	 
!    relaxation
        do i=2,mx-1
	       do j=2,my-1
		      if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then  ! 正方形内部は計算しない
		         dp=(p(i+1,j)+p(i-1,j))/dx**2+(p(i,j+1)+p(i,j-1))/dy**2-rhs(i,j)
		         dp=dp/(2.0/dx**2+2.0/dy**2)-p(i,j)  ! p(i,j)の修正すべき量
		         res=res+dp**2
		         p(i,j)=p(i,j)+omegap*dp ! 修正すべき量に緩和係数をかけて修正
			  end if	 
		   end do
        end do
!    set boundary condition	 
	    call bcforp
!
        res=sqrt(res/float(mx*my))  ! 残差(これが0または十分小さいなら収束)
!
        if(res<errorp) exit
     end do
!     
	 resp=res
     itrp=itr
!
  end subroutine poiseq
!
!======================
  subroutine veloeq
!======================
!
! 速度場を解く
!
     use cmpcnd
	 use grid2d
	 use flowvp
	 use worksp
	 use others
!
     implicit none
!
!
! (1)pressure gradient
     do i=2,mx-1
	    do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then
		      urhs(i,j)=-(p(i+1,j)-p(i-1,j))/(2.0*dx)
		      vrhs(i,j)=-(p(i,j+1)-p(i,j-1))/(2.0*dy)
		   end if
		end do
     end do
! (2)viscous term
     do i=2,mx-1
	    do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then
		      urhs(i,j)=urhs(i,j)                                         &
		                +(u(i+1,j)-2.0*u(i,j)+u(i-1,j))/(re*dx**2)        &
					    +(u(i,j+1)-2.0*u(i,j)+u(i,j-1))/(re*dy**2)
              vrhs(i,j)=vrhs(i,j)                                         &
 		                +(v(i+1,j)-2.0*v(i,j)+v(i-1,j))/(re*dx**2)        &
					    +(v(i,j+1)-2.0*v(i,j)+v(i,j-1))/(re*dy**2)      
           end if						
        end do
     end do
! (3)Advection term in x-direction
! kawamuraスキームは左右2点の計5点の値が必要。物体内部や計算領域外部のものも必要であり、適当に外挿するしかない。
     do j=j1+1,j2-1
	    u(i1+1,j)=2.0*u(i1,j)-u(i1-1,j)
		u(i2-1,j)=2.0*u(i2,j)-u(i2+1,j)
		v(i1+1,j)=2.0*v(i1,j)-v(i1-1,j)
		v(i2-1,j)=2.0*v(i2,j)-v(i2+1,j)
	 end do
     do j=j3+1,j4-1
        u(i1+1,j)=2.0*u(i1,j)-u(i1-1,j)
        u(i2-1,j)=2.0*u(i2,j)-u(i2+1,j)
        v(i1+1,j)=2.0*v(i1,j)-v(i1-1,j)
        v(i2-1,j)=2.0*v(i2,j)-v(i2+1,j)
     end do
	 do i=2,mx-1
	    do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3))  then
		      urhs(i,j)=urhs(i,j)                                                                        &
		                -u(i,j)*(-u(i+2,j)+8.0*(u(i+1,j)-u(i-1,j))+u(i-2,j))/(12.0*dx)                   & 
					    -abs(u(i,j))*(u(i+2,j)-4.0*u(i+1,j)+6.0*u(i,j)-4.0*u(i-1,j)+u(i-2,j))/(4.0*dx)   
              vrhs(i,j)=vrhs(i,j)                                                                        &
		                -u(i,j)*(-v(i+2,j)+8.0*(v(i+1,j)-v(i-1,j))+v(i-2,j))/(12.0*dx)                   &
		                -abs(u(i,j))*(v(i+2,j)-4.0*v(i+1,j)+6.0*v(i,j)-4.0*v(i-1,j)+v(i-2,j))/(4.0*dx)  
           end if						
        end do
     end do
! (4)Advecion term in y-direction
     do i=i1+1,i2-1
	    u(i,j1+1)=2.0*u(i,j1)-u(i,j1-1)
		u(i,j2-1)=2.0*u(i,j2)-u(i,j2+1)
		v(i,j1+1)=2.0*v(i,j1)-v(i,j1-1)
		v(i,j2-1)=2.0*v(i,j2)-v(i,j2+1)
	 end do
     do i=i1+1,i2-1
        u(i,j4-1)=2.0*u(i,j4)-u(i,j4+1)
        u(i,j3+1)=2.0*u(i,j3)-u(i,j3-1)
        v(i,j4-1)=2.0*v(i,j4)-u(i,j4+1)
        v(i,j3+1)=2.0*v(i,j3)-u(i,j3-1)
     end do
     do i=2,mx-1
        do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then
		      urhs(i,j)=urhs(i,j)                                                                        &
		                -v(i,j)*(-u(i,j+2)+8.0*(u(i,j+1)-u(i,j-1))+u(i,j-2))/(12.0*dy)                   &
		                -abs(v(i,j))*(u(i,j+2)-4.0*u(i,j+1)+6.0*u(i,j)-4.0*u(i,j-1)+u(i,j-2))/(4.0*dy)	  	   
              vrhs(i,j)=vrhs(i,j)                                                                        &
		                -v(i,j)*(-v(i,j+2)+8.0*(v(i,j+1)-v(i,j-1))+v(i,j-2))/(12.0*dy)                   &
		                -abs(v(i,j))*(v(i,j+2)-4.0*v(i,j+1)+6.0*v(i,j)-4.0*v(i,j-1)+v(i,j-2))/(4.0*dy)  
           end if						
        end do
     end do
! (5)update
     do i=2,mx-1
	    do j=2,my-1
		   if(i<i1.or.i>i2.or.j<j1.or.j>j4.or.(j>j2.and.j<j3)) then
		      u(i,j)=u(i,j)+dt*urhs(i,j)
		      v(i,j)=v(i,j)+dt*vrhs(i,j)
		   end if
		end do
     end do
!
  end subroutine veloeq
 
	 

  
	 
		   
		   
     
	 
		













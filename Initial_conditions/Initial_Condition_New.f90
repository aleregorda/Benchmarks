program GM
 use WorldBuilder

 implicit none

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 integer, parameter :: nodi_elem=4
 double precision, parameter :: T_aria=273.15D0
 double precision, parameter :: z=0.D0
 double precision, parameter :: depth=1.D-04
 double precision, parameter :: gravity=9.81D0 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VARIABILI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 integer :: i,j,jj,k,ii,elementi,nodi_x,nodi_y,nodi,length,elem_x,elem_y,marc_elem,nmats,nx,indice,&
 elem_x1,elem_x2,elem_x3

 double precision :: xm,ym,temp,delta_x,delta_y,x_min,x_max,y_min,y_max,h,elem_xd,elem_yd,&
 temperature,composition,delta_x1,delta_x2,delta_x3,var_grid

 integer, dimension(:), allocatable :: markers
 double precision, dimension(:), allocatable :: x,y,xmin_elem,ymin_elem,deltax,deltay
 integer, dimension(:,:), allocatable :: nodo,connessioni
 character(len=:),allocatable :: file_name
 character(len=50) :: filename
 character(len=50) :: output_dir = "../../doc/manual/"//C_NULL_CHAR
 character(len=1) :: input,rdm,inputG
 logical :: random
 logical(1) :: has_output_dir = .false.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     INPUT
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 random=.false.
 print *,''
 print *,'Utilizzare griglia esistente? (y/n)'
 read *,input
 if (input.eq.'y') goto 1000
 
 print *,'Inserire nome file per le condizioni iniziali'
 read *,filename
 length=LEN(TRIM(filename))+20
 allocate(character(length) :: file_name)
 file_name="./Initial_conditions/"//TRIM(ADJUSTL(filename))//C_NULL_CHAR

100 print *,'Inserire X minima (in km)'
 read *,x_min
 x_min=x_min*1.D03
 print *,'Inserire X massima (in km)'
 read *,x_max
 x_max=x_max*1.D03
 if(x_max.le.x_min) then
   print *,'La X massima deve essere maggiore della X minima'
   goto 100
 endif
200 print *,'Inserire Y minima (in km)'
 read *,y_min
 y_min=y_min*1.D03
 print *,'Inserire Y massima (in km)'
 read *,y_max
 y_max=y_max*1.D03
 if(y_max.le.y_min) then
   print *,'La Y massima deve essere maggiore della Y minima'
   goto 200
 endif
300 print *,'Inserire spessore aria (in km)'
 read *,h
 h=h*1.D03
350 print *,'Scegliere se griglia regolare (1) o variabile (2)'
 read *,inputG
 if(inputG.ne.'1'.and.inputG.ne.'2')then
   print *,'Scelta non corretta'
   goto 350
 endif
 if(inputG.eq.'1')goto 400
 if(inputG.eq.'2')goto 700
400 print *,'Scegliere se inserire il numero (1) o la dimensione (2) degli elementi'
 read *,input
 if(input.ne.'1'.and.input.ne.'2')then
   print *,'Scelta non corretta'
   goto 400
 endif
 if(input.eq.'1')then
   print *,'Inserire elementi lungo X'
   read *,elem_x
   elem_x=INT(elem_x)
   print *,'Inserire elementi lungo Y'
   read *,elem_y
   elem_y=INT(elem_y)
   elementi=elem_x*elem_y
   delta_x=((x_max-x_min)/elem_x)
   delta_y=((y_max-y_min)/elem_y)
 elseif(input.eq.'2')then
500 print *,'Inserire dimensione elementi lungo X (in km)'
   read *,delta_x
   delta_x=delta_x*1.D03
   elem_xd=(x_max-x_min)/delta_x
   if(elem_xd.eq.INT(elem_xd))then
     elem_x=INT(elem_xd)
   else
     print *,'Dimensione non corretta'
     goto 500
   endif
600 print *,'Inserire dimensione elementi lungo Y (in km)'
   read *,delta_y
   delta_y=delta_y*1.D03
   elem_yd=(y_max-y_min)/delta_y
   if(elem_yd.eq.INT(elem_yd))then
     elem_y=INT(elem_yd)
   else
     print *,'Dimensione non corretta'
     goto 600
   endif
   elementi=((x_max-x_min)/delta_x)*((y_max-y_min)/delta_y)
 endif
 goto 900
700 print *,'Inserire dimensione minima elementi lungo X (in km)'
 read *,delta_x3
 delta_x3=delta_x3*1.D03
 delta_x2=5.D0*delta_x3
 delta_x1=2.D0*delta_x2
 var_grid=x_max/5.D0
 elem_x1=INT(var_grid/delta_x1)
 elem_x2=INT(var_grid/delta_x2)
 elem_x3=INT(var_grid/delta_x3)
 elem_x=elem_x1*2+elem_x2*2+elem_x3
800 print *,'Inserire dimensione elementi lungo Y (in km)'
 read *,delta_y
 delta_y=delta_y*1.D03
 elem_yd=(y_max-y_min)/delta_y
 if(elem_yd.eq.INT(elem_yd))then
   elem_y=INT(elem_yd)
 else
   print *,'Dimensione non corretta'
   goto 800
 endif
 elem_y=(y_max-y_min)/delta_y
 elementi=elem_x*elem_y
900 print *,'Inserire numero marcatori per elemento (in caso di griglia variabile considerare elemento minore)'
 read *,marc_elem
 marc_elem=ABS(INT(marc_elem))
 print *,'Utilizzare una distribuzione random dei marcatori? (y/n)'
 read *,rdm
 if(rdm.eq.'y')random=.true.
 print *,'Inserire numero di materiali (escluso aria)'
 read *,nmats
 nmats=INT(nmats)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     CREAZIONE GRIGLIA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 print *,''
 print *,'****************************************************************************************'
 print *,'---------------------------- CREATING GRID AND MARKERS FILES ---------------------------'
 print *,'****************************************************************************************'
 print *,''

 nodi_x=elem_x+1
 nodi_y=elem_y+1
 nodi=nodi_x*nodi_y

 print *,'Creating grid constituted by:'
 print 1001,nodi
 print 1002,elementi
 print *,''

 allocate(nodo(nodi_x,nodi_y))
 allocate(x(nodi))
 allocate(y(nodi))
 allocate(connessioni(elementi,nodi_elem))
 allocate(xmin_elem(elementi))
 allocate(ymin_elem(elementi))
 allocate(deltax(elementi))
 allocate(deltay(elementi))
 allocate(markers(elementi))

 do j=1,nodi_y
   do i=1,nodi_x
     indice=nodi_x*(j-1)+i
     if(j.eq.1)then
       y(indice)=y_min
     else
       y(indice)=y(indice-nodi_x)+delta_y
     endif
     if(i.eq.1)then
       x(indice)=x_min
     else
       if(inputG.eq.'1')then
         x(indice)=x(indice-1)+delta_x
         delta_x3=delta_x
       else
         if(x(indice-1).lt.var_grid.or.x(indice-1).ge.var_grid*4.D0)then
           x(indice)=x(indice-1)+delta_x1
         elseif(x(indice-1).lt.var_grid*2.D0.or.x(indice-1).ge.var_grid*3.D0)then
           x(indice)=x(indice-1)+delta_x2
         else
           x(indice)=x(indice-1)+delta_x3
         endif
       endif
     endif
     nodo(i,j)=indice
   enddo
 enddo

 do j=1,elem_y
   do i=1,elem_x
     indice=elem_x*(j-1)+i
     connessioni(indice,1)=nodo(i,j)
     connessioni(indice,2)=nodo(i+1,j)
     connessioni(indice,3)=nodo(i+1,j+1)
     connessioni(indice,4)=nodo(i,j+1)
     xmin_elem(indice)=MINVAL(x(connessioni(indice,:)))
     ymin_elem(indice)=MINVAL(y(connessioni(indice,:)))
     deltax(indice)=MAXVAL(x(connessioni(indice,:)))-MINVAL(x(connessioni(indice,:)))
     deltay(indice)=MAXVAL(y(connessioni(indice,:)))-MINVAL(y(connessioni(indice,:)))
     if(.not.random)then
       markers(indice)=INT(SQRT(DBLE(marc_elem*(deltax(indice)/delta_x3))))**2
     else
       markers(indice)=marc_elem*(deltax(indice)/delta_x3)
     endif
   enddo
 enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     SCRITTURA FILE GRIGLIA+TEMPERATURA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 call create_world(cworld, file_name, has_output_dir, output_dir)

 print *,'Writing grid file'
 print *,''
 open(unit=1,file='Griglia',status='unknown') 
 rewind(1)

 if(inputG.eq.'1')then
   write(1,*)'Fixed'
   write(1,*)0,0,0,0.D0
 else
   write(1,*)'Variable'
   write(1,*)elem_x1,elem_x2,elem_x3,var_grid
 endif
 write(1,*)'Spessore aria'
 write(1,*)h
 write(1,*)'Nodi'
 write(1,*)nodi,nodi_x,nodi_y
 do j=1,nodi_y
   do i=1,nodi_x
     indice=nodi_x*(j-1)+i
     if(y_max-h-y(indice).ge.y_min)then
       call temperature_2d(cworld,x(indice),depth,y_max-h-y(indice),gravity,temperature)
       write(1,'(i10,1x,2(f15.4,1x),f15.4,1x,i4)')indice,x(indice),y(indice),temperature
     else
       write(1,'(i10,1x,2(f15.4,1x),f15.4,1x,i4)')indice,x(indice),y(indice),T_aria
     endif
   enddo
 enddo
   
 write(1,*)'Elementi'
 write(1,*)elementi,elem_x,elem_y
 do j=1,elem_y
   do i=1,elem_x
     indice=elem_x*(j-1)+i
     write(1,'(5(i10,1x))')indice,nodo(i,j),nodo(i+1,j),nodo(i+1,j+1),nodo(i,j+1)
   enddo
 enddo
 close(1)

 deallocate(nodo)
 deallocate(connessioni)
 deallocate(x)
 deallocate(y)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     SCRITTURA FILE MARCATORI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
 print 1003,SUM(markers)
 open(unit=1,file='Marcatori',status='unknown') 
 rewind(1)
 write(1,*)SUM(markers)
 k=0
 do i=1,elementi
   if(random)then
     do j=1,markers(i)
       call random_number(xm)
       call random_number(ym)
       xm=(xm*deltax(i))+xmin_elem(i)
       ym=(ym*deltay(i))+ymin_elem(i)
       k=k+1
       if(y_max-h-ym.ge.y_min)then
         do ii=1,nmats
           call composition_2d(cworld,xm,depth,y_max-h-ym,ii,composition)
           if(composition.eq.1.D0)then
             write(1,'(i10,2(f15.4,1x),i10)')k,xm,ym,ii
             exit
           endif
         enddo
       else
         write(1,'(i10,2(f15.4,1x),i10)')k,xm,ym,nmats+1
       endif
     enddo
   else
     nx=INT(SQRT(DBLE(markers(i))))
     do j=1,nx
     do jj=1,nx
       xm=xmin_elem(i)+(deltax(i)/(nx+1))*j
       ym=ymin_elem(i)+(deltay(i)/(nx+1))*jj
       k=k+1
       if(y_max-h-ym.ge.y_min)then
         do ii=1,nmats
           call composition_2d(cworld,xm,depth,y_max-h-ym,ii,composition)
           if(composition.eq.1.D0)then
             write(1,'(i10,2(f15.4,1x),i10)')k,xm,ym,ii
             exit
           endif
         enddo
       else
         write(1,'(i10,2(f15.4,1x),i10)')k,xm,ym,nmats+1
       endif
     enddo
     enddo
   endif
 enddo
 close(1)

 deallocate(xmin_elem)
 deallocate(ymin_elem)

 call release_world(cworld)

 print *,''
 print *,'Done!'
 print *,''
 print *,'****************************************************************************************'
 print *,''

1000 continue

1001 format(' -',i8,' nodes')
1002 format(' -',i8,' elements')
1003 format(' Writing markers file (',i8,' markers)')

end program

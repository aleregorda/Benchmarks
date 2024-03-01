   if method==3: # aspect ALE

      print('computing vmesh for bc with Q1 elements')

      #first compute the L2 projection of the mesh velocity

      A_fs = np.zeros(((nelx+1)*ndim,(nelx+1)*ndim),dtype=np.float64)
      b_fs = np.zeros((nelx+1)*ndim,dtype=np.float64)
      A_el =np.zeros((4,4),dtype=np.float64)
      b_el =np.zeros((4),dtype=np.float64)
      iconFS=np.zeros((2,nelx),dtype=np.int32)
      umesh_fs = np.zeros((nelx+1),dtype=np.float64)
      vmesh_fs = np.zeros((nelx+1),dtype=np.float64)
      x_FS = np.zeros((nelx+1),dtype=np.float64)
      y_FS = np.zeros((nelx+1),dtype=np.float64)

      counter=0
      for j in range(0,nely):
          for i in range(0,nelx):
              if j==nely-1: # if element is top row

                 # compute normal to element edge
                 x1=xV[iconV[3,counter]] ; y1=yV[iconV[3,counter]]
                 x2=xV[iconV[2,counter]] ; y2=yV[iconV[2,counter]]
                 n_x=-y2+y1
                 n_y= x2-x1
                 nnorm=np.sqrt(n_x**2+n_y**2)
                 n_x/=nnorm
                 n_y/=nnorm
                 d12=np.sqrt((x1-x2)**2+(y1-y2)**2)
                 jcob=d12/2

                 #build connectivity array on the fly
                 iconFS[0,i]=i+0 
                 iconFS[1,i]=i+1 

                 x_FS[iconFS[0,i]]=xV[iconV[3,counter]]
                 x_FS[iconFS[1,i]]=xV[iconV[2,counter]]
                 y_FS[iconFS[0,i]]=yV[iconV[3,counter]]
                 y_FS[iconFS[1,i]]=yV[iconV[2,counter]]


                 A_el[:,:]=0.
                 b_el[:]=0.
                 for iq in [-1,1]:
                     rq=iq/np.sqrt(3.)
                     weightq=1. ; JxW=weightq*jcob
                     N1= 0.5*(1-rq) 
                     N2= 0.5*(1+rq) 
                     xq=N1*xV[iconV[3,counter]]+N2*xV[iconV[2,counter]]
                     yq=N1*yV[iconV[3,counter]]+N2*yV[iconV[2,counter]]
                     uq=N1*u[iconV[3,counter]]+N2*u[iconV[2,counter]]
                     vq=N1*v[iconV[3,counter]]+N2*v[iconV[2,counter]]
                     bob=uq*n_x+vq*n_y
                     #print (xq,bob,'aaa')
                     A_el[0,0]+=N1*N1*JxW ; A_el[0,2]+=N1*N2*JxW 
                     A_el[1,1]+=N1*N1*JxW ; A_el[1,3]+=N1*N2*JxW 
                     A_el[2,0]+=N2*N1*JxW ; A_el[2,2]+=N2*N2*JxW 
                     A_el[3,1]+=N2*N1*JxW ; A_el[3,3]+=N2*N2*JxW 
                     b_el[0]+=bob*N1*n_x*JxW
                     b_el[1]+=bob*N1*n_y*JxW
                     b_el[2]+=bob*N2*n_x*JxW
                     b_el[3]+=bob*N2*n_y*JxW
                 #end for iq
      
                 # assemble matrix A_fs and rhs b_fs
                 for k1 in range(0,2):
                     for i1 in range(0,ndofV):
                         ikk=ndofV*k1          +i1
                         m1 =ndofV*iconFS[k1,i]+i1
                         for k2 in range(0,2):
                             for i2 in range(0,ndofV):
                                 jkk=ndofV*k2          +i2
                                 m2 =ndofV*iconFS[k2,i]+i2
                                 A_fs[m1,m2] += A_el[ikk,jkk]
                             #end for
                         #end for
                         b_fs[m1]+=b_el[ikk]
                     #end for i1
                 #end for k1
              #end if
              counter+=1
          #end for
      #end for

      print("     -> A_fs (m,M) %.4e %.4e " %(np.min(A_fs),np.max(A_fs)))
      print("     -> b_fs (m,M) %.4e %.4e " %(np.min(b_fs),np.max(b_fs)))

      sparse_matrix=sps.csr_matrix(A_fs)
      sol=sps.linalg.spsolve(sparse_matrix,b_fs)
      umesh_fs,vmesh_fs=np.reshape(sol[0:2*(nelx+1)],(nelx+1,2)).T

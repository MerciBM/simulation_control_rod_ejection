# Capitalization of the code
# 2026-09-01
# Pr Bertrand MERCIER bertrand.mercier@cea.fr

import numpy as np
from os import path

if __name__ == '__main__':
    main_path = path.abspath(__file__)
    src_path = path.dirname(main_path)
    folder_path = path.dirname(src_path)
    input_path = path.join(folder_path, 'input')
    output_path = path.join(folder_path, 'output')
    # This is the only input file. It contains the starting critical states
    # it is the output of another program that is not included here
    path_sortie_eigensh3D = path.join(input_path, 'sortie eigensh3D.txt')
    # These are the outputs.
    # The first contains the final state
    path_sortie_yy = path.join(output_path, 'sortieYY.txt')
    # The second one contains the evolution
    path_sortie = path.join(output_path, 'sortie.txt')
    # This is a 3D version  
    # The control rod is extracted step by step
    Sigmaf=0.24
    nue=2.
    Sigma0=0.472565
    # this value has been slected in such a way that the initial state
    # is a critical state
    Diff=10.
    h=20.
    gama=Diff/h**2
    print("gama=",gama)
    vit=220000.
    # thermal neutrons speed (cm/s)
    dt=.00001
    beta=.0045
    mu=.08
    ZZ=49/73*np.loadtxt(path_sortie_eigensh3D)
    # contains the initial critical state
    XX = np.zeros([9,9,9], dtype=float)
    i12=0
    for i1 in range(9) :
        for i2 in range (9) :
            for i3 in range(9):
                XX[i1,i2,i3]=ZZ[i12,i3]/9.
            i12=i12+1
    RHO=np.zeros([9,9,9], dtype=float)
    SIG0=Sigma0*np.ones([9,9,9], dtype=float)
    for i3 in range(8):
        SIG0[4,4,8-i3]=4.8
    # 
    # XX[i,j,k] is the power in layer k of assembly n°ij, in MW
    # 
    # Wcore is then the initial power of the core in MW
    # 
    Wcore=np.sum(XX)
    print(" Wcore=",Wcore)
    # moderator temperature (it will not change)
    Tmod=300.
    hs=7.
    hsij=hs/729.
    # fuel -> moderator transfer coefficient in one layer of assembly i,j  
    TF=np.zeros([9,9,9], dtype=float)
    TF=Tmod+XX/hsij
    TF0=TF
    print("average =",np.mean(TF))
    # normalized population of precursors
    CQ=beta*vit*nue*Sigmaf*XX/mu 
    # fuel heat capacity of one layer of one assembly MJ/°
    mc=0.005
    # à t=0.
    # Doppler pcm/°
    alfaD=-3.
    alfad=alfaD/100000.
    t=0.
    nit=90
    tab=np.zeros([100,9])
    for it in range(9000) :
        AA=XX-hsij*(TF-Tmod)
        Tfa=np.mean(TF)
        TFnew=TF+dt*AA/mc
        RHOnew=RHO+alfad*(TFnew-TF)
        YY = np.zeros([9,9,9], dtype=float)
        t=t+dt
            #  matrix vector product
        for k in range(9):
            for j in range(9):
                for i in range(9):
                    Cdiag=6.
                    XW=0.
                    if i>0 :
                        XW=XX[i-1,j,k]
                    XS=0.
                    if j>0 :
                        XS=XX[i,j-1,k]    
                    XE=0.
                    if i<8 :
                        XE=XX[i+1,j,k]
                    XN=0.
                    if j<8 :
                        XN=XX[i,j+1,k]  
                    XL=0.
                    if k>0 :
                        XL=XX[i,j,k-1]
                    else :
                        Cdiag=5.
                    XU=0.
                    if k<8 :
                        XU=XX[i,j,k+1]
                    else :
                        Cdiag=5.
                    YY[i,j,k]=Cdiag*XX[i,j,k]-XW-XS-XE-XN-XL-XU
        #   first step of the splitting method
        XX=XX-dt*vit*gama*YY
        # cross sections in the central fuel assembly are modified 
        # to take into account control rod extraction
        if it<9000 :
            k=int(it/1000.)
            theta=(it-k*1000.)/1000.
            SIG0[4,4,k]=theta*Sigma0+(1.-theta)*4.8
        # second step of the splitting method
        for k in range(9):
            for j in range(9):
                for i in range(9):
                    a=1.-dt*vit*(nue*Sigmaf*(1.-beta)-SIG0[i,j,k]*(1.-RHOnew[i,j,k]))
                    b=-mu*dt
                    c=-dt*beta*vit*nue*Sigmaf
                    d=1.+mu*dt
                    det=a*d-c*b
                    xx=(d*XX[i,j,k]-b*CQ[i,j,k])/det
                    cc=(-c*XX[i,j,k]+a*CQ[i,j,k])/det
                    XX[i,j,k]=xx
                    CQ[i,j,k]=cc
        TF=TFnew
        RHO=RHOnew
        Wcore=np.sum(XX)
        if it==int(it/nit)*nit :
            jt=int(it/nit)
            tab[jt,0]=t-dt
            tab[jt,1]=XX[2,2,2]
            tab[jt,2]=TF[2,2,2]
            tab[jt,3]=1.E5*RHO[2,2,2]
            tab[jt,4]=XX[4,4,4]
            tab[jt,5]=TF[4,4,4]
            tab[jt,6]=1.E5*RHO[4,4,4]
            tab[jt,7]=CQ[4,4,4]
            tab[jt,8]=Wcore
            print("t=",t-dt," XX[4,4,4]=",XX[4,4,4],"RHO[4,4,4]",1E5*RHO[4,4,4])
    i12=0
    for i1 in range(9) :
        for i2 in range (9) :
            for i3 in range(9):
                ZZ[i12,i3]=XX[i1,i2,i3]
            i12=i12+1
    np.savetxt(path_sortie_yy,ZZ,fmt='%10.5f',delimiter=" ")
    np.savetxt(path_sortie,tab,fmt='%10.5f',delimiter=" ")
    Wcore=np.sum(XX)
    print(" Wcore=",Wcore)
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
    path_sortie_eigensh4 = path.join(input_path, 'sortieEIGENSH4.txt')
    # These are the outputs.
    # The first one contains the reactivity
    path_sortie_rho = path.join(output_path, 'sortieRHO.txt')
    # The second contains the final state
    path_sortie_yy = path.join(output_path, 'sortieYY.txt')
    # The third one contains the evolution
    path_sortie = path.join(output_path, 'sortie.txt')


    # We shall proceed to a local reactivity injection
    #  rho0 such that (1-rho0)*10=1
    # 
    rho0=1.-1./10.
    # drho will be injected at each time step for 10 000 time steps
    drho=rho0/10000
    print(" drho=",drho)
    lam0=0.195773934819386
    # lam0 is the eigenvalue of matrix A
    KINF0=1.01
    gs=(KINF0-1)/lam0
    print(" gs=",gs)
    Sigmaf=0.24
    nue=2.
    Sigma0=nue*Sigmaf/KINF0
    gama=gs*Sigma0
    print("gama=",gama)
    vit=220000.
    # thermal neutrons speed (cm/s)
    dt=.00001
    # timestep (s)
    beta=.006
    # delayed neutrons fraction
    mu=.08
    # decay time constant of precursors (s^-1)
    XX=0.97*np.loadtxt(path_sortie_eigensh4)
    # XX is the initial power map (with inserted rod)
    RHO=np.zeros([9,9], dtype=float)
    RHO[4,4]=0.
    SIG0=0.9948*Sigma0*np.ones([9,9], dtype=float)
    SIG0[4,4]=Sigma0*10.
    # Assembly 4,4 is the central assembly where the rod is initially inserted
    # at t=0 rho = 0 in all 81 assemplies
    # 
    # XX[i,j] is the power (MW) of assembly n°ij (MW)
    # Since the sum is about equal to  40 
    # It means that Wcore = 40 MW.
    # 
    Wcore=np.sum(XX)
    print(" Wcore=",Wcore)
    # moderator temperature (°Celsius)
    Tmod=300.
    hs=7.
    #  fuel -> moderator transfer coefficient (MW/K)
    hsij=hs/81.
    #  same but for one assembly
    TF=np.zeros([9,9], dtype=float)
    # Fuel temperature of assemblies
    TF=Tmod+XX/hsij
    print(" starting point ")
    print(XX)
    print(" initial Tfuel ")
    print (TF)
    TF0=TF
    print("average =",np.mean(TF))
    CQ=beta*vit*nue*Sigmaf*XX/mu
    print("initial concentration of precursors")
    print(CQ)
    # fuel heat capacity for one assembly MJ/°
    mc=0.01
    # Doppler coefficient pcm/°
    alfaD=-3.
    alfad=alfaD/100000.
    t=0.
    # nb. of timesteps between outputs
    nit=200
    tab=np.zeros([100,9])
    # tab will contain the outputs
    for it in range(20000) :
        AA=XX-hsij*(TF-Tmod)
        Tfa=np.mean(TF)
        TFnew=TF+dt*AA/mc
        RHOnew=RHO+alfad*(TFnew-TF)
        # local reactivity will decrease due to the Doppler effect
        YY = np.zeros([9,9], dtype=float)
        # YY is an intermediate array
        t=t+dt
            # matrix A times vector product
        for j in range(9):
            for i in range(9):
                XW=0.
                if i>0 :
                    XW=XX[i-1,j]
                XS=0.
                if j>0 :
                    XS=XX[i,j-1]    
                XE=0.
                if i<8 :
                    XE=XX[i+1,j]
                XN=0.
                if j<8 :
                    XN=XX[i,j+1]  
                YY[i,j]=4.*XX[i,j]-XW-XS-XE-XN
        XX=XX-dt*vit*gama*YY
        for j in range(9):
            for i in range(9):
            # coefficients of the 2x2 matrix to be inverted
                a=1.-dt*vit*(nue*Sigmaf*(1.-beta)-SIG0[i,j]*(1.-RHOnew[i,j]))
                b=-mu*dt
                c=-dt*beta*vit*nue*Sigmaf
                d=1.+mu*dt
                det=a*d-c*b
                xx=(d*XX[i,j]-b*CQ[i,j])/det
                cc=(-c*XX[i,j]+a*CQ[i,j])/det
                XX[i,j]=xx
                CQ[i,j]=cc
        TF=TFnew
        RHO=RHOnew
        if it<10000 :
            RHO[4,4]=RHO[4,4]+drho
        Wcore=np.sum(XX)
        if it==int(it/nit)*nit :
            jt=int(it/nit)
            tab[jt,0]=t-dt
            tab[jt,1]=XX[2,2]
            tab[jt,2]=TF[2,2]
            tab[jt,3]=1.E5*RHO[2,2]
            tab[jt,4]=XX[4,4]
            tab[jt,5]=TF[4,4]
            tab[jt,6]=1.E5*RHO[4,4]
            tab[jt,7]=CQ[4,4]
            tab[jt,8]=Wcore
            print("t=",t-dt," XX[4,4]=",XX[4,4],"RHO[4,4]",1E5*RHO[4,4])
    print("XX=",XX)
    print("RHO=",RHO*1E5)
    print("Tmod=",Tmod)
    np.savetxt(path_sortie_rho,RHO,fmt='%10.5f',delimiter=" ")
    np.savetxt(path_sortie_yy,XX,fmt='%10.5f',delimiter=" ")
    np.savetxt(path_sortie,tab,fmt='%10.5f',delimiter=" ")
    Wcore=np.sum(XX)
    print(" Wcore=",Wcore)
# version 2023_07_31 
RDT = function(y, J, pois = rep(10,5), t_pois = rep(1:5,each=10), nchain=2, iters=15000, burnin=5000, thin=3, ReTrain=F, mod=NULL){
  
  cat("In preparation.\n")
  
  # MCMC
  ss = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I/2,b+sum(x^2)/2)
    return(s)
  }
  
  ss2 = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I,b+sum(x))
    return(s)
  }
  
  tdoor = function(door){
    door = array(0,dim=c(length()))
  }
  
  Prob_RD = function(delta,tau,a1,ga,theta1,theta2){
    N = length(theta1)
    J = length(tau)+1
    
    # RD
    trd = (a1*theta1+exp(theta2+ga)-delta)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_trd = matrix(trd, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp(m_trd - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob_rd = prob / c_Prob
    
    return(Prob_rd)
  }
  
  Prob_T = function(tau,w1,theta3){
    
    N = length(theta3)
    J = length(tau)+1
    
    # T
    tt = matrix(theta3+w1,N,J-1,byrow=F)
    mw2 = matrix(tau,N,J-1,byrow=T)
    door = exp(-(tt-mw2)*(tt>(mw2)))
    door[which(is.na(door)==T)] = 1
    
    return(door)
  }
  
  Prob_RDT = function(Prob_rd,door,N,J){
    
    # RDT
    d = (door[,1:(J-1)] - cbind(0,door[,1:(J-2)])) # 11, 12, ..., 1J-1 
    d2 = cbind((matrix(1,N,J-1) - door[,(J-1):1]),1) # 1J
    Prob_rdt = NULL
    for(j in 1:(J-1)){
      Prob_t = d[,j:1]
      tmp =  matrix((Prob_rd[,1:j] * Prob_t),N,j,byrow=F)
      Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    }
    tmp = matrix((Prob_rd * d2),N,J,byrow=F)
    Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    
    return(Prob_rdt)
  }
  
  Prob_RD_theta = function(delta,tau,a1,ga,theta1,theta2,t_pois,R=T){
    
    I = length(delta)
    J = length(tau)+1
    
    tt1 = ifelse(rep(R,I),theta1,theta1[t_pois])
    
    # RD
    trd = (a1*tt1+exp(theta2+ga)-delta)
    m_I = matrix(0:(J-1), I, J, byrow=T)
    m_trd = matrix(trd, I, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, I, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp(m_trd - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), I,J, byrow=F)
    Prob_rd = prob / c_Prob
    
    return(Prob_rd)
  }
  
  Prob_T_theta = function(tau,w1,theta3){
    
    I = length(w1)
    J = length(tau)+1
    
    # T
    tt = matrix(theta3+w1,I,J-1,byrow=F)
    mw2 = matrix(tau,I,J-1,byrow=T)
    door = exp(-(tt-mw2)*(tt>(mw2)))
    door[which(is.na(door)==T)] = 1
    
    return(door)
  }
  
  Prob_RDT_theta = function(Prob_rd,door,I,J){
    
    # RDT
    d = (door[,1:(J-1)] - cbind(0,door[,1:(J-2)])) # 11, 12, ..., 1J-1 
    d2 = cbind((matrix(1,I,J-1) - door[,(J-1):1]),1) # 1J
    Prob_rdt = NULL
    for(j in 1:(J-1)){
      Prob_t = d[,j:1]
      tmp =  matrix((Prob_rd[,1:j] * Prob_t),I,j,byrow=F)
      Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    }
    tmp = matrix((Prob_rd * d2),I,J,byrow=F)
    Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    
    return(Prob_rdt)
  }
  
  pdelta = function(delta,Sd,yy,Prob_rdt){
    
    N = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=N, ncol=J, byrow=F) == matrix(0:(J-1), nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    ((-1/2)*((delta)^2/Sd)+sum(log(sum_Prob)))
  }
  
  pa1 = function(a1,Sa1,yy,Prob_rdt){
    
    N = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=N, ncol=J, byrow=F) == matrix(0:(J-1), nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    (-a1/Sa1+sum(log(sum_Prob)))
  }
  
  pa2 = function(a2,Sa2,yy,Prob_rdt){
    
    N = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=N, ncol=J, byrow=F) == matrix(0:(J-1), nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    (-0.5*((a2^2)/Sa2)+sum(log(sum_Prob)))
  }
  
  pw1 = function(w1,Sw,yy,Prob_rdt){
    
    N = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=N, ncol=J, byrow=F) == matrix(0:(J-1), nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    -0.5*((w1)^2/Sw)+sum(log(sum_Prob))
  }
  
  
  ptau = function(tau,St,yy,Prob_rdt,J){
    
    Prob = 0
    for(j in 0:(J-1)){
      m_y = (y == j)
      Prob = Prob + Prob_rdt[,,j+1] * m_y
    }
    
    (-sum(tau^2)/(2*St)+sum(log(Prob)))
  }
  
  ptheta1 = function(theta1,S,yy,Prob_rdt,theta3){
    
    I = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=I, ncol=J, byrow=F) == matrix(0:(J-1), nrow=I, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    tt = c(theta1, theta3)
    (-0.5*t(tt)%*%solve(S)%*%tt+sum(log(sum_Prob)))
  }
  
  ptheta2 = function(theta2,Stt2,yy,Prob_rdt){
    
    I = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=I, ncol=J, byrow=F) == matrix(0:(J-1), nrow=I, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    (-0.5*sum(theta2^2)/(Stt2)+sum(log(sum_Prob)))
  }

  ptheta3 = function(theta3,S,yy,Prob_rdt,theta1){
    
    I = nrow(Prob_rdt)
    J = ncol(Prob_rdt)
    
    m_y = matrix(yy, nrow=I, ncol=J, byrow=F) == matrix(0:(J-1), nrow=I, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    tt = c(theta1, theta3)
    (-0.5*t(tt)%*%solve(S)%*%tt+sum(log(sum_Prob)))
  }
  
  plambda = function(l,t1,t3,la,v0,d0){
    D = diag(l)
    la = as.vector(la)
    tmp = diag(1,length(l))
    tmp[upper.tri(tmp)] = la
    LA = t(tmp) 
    S = solve(LA) %*% D %*% t(solve(LA))
    
    N = length(t3)
    theta = t(cbind(t1,t3))
    (sum((-v0-1)*log(l)-d0/l)+
        -(N/2)*log(det(S))-0.5*psych::tr(t(theta)%*%solve(S)%*%theta))
  }
  
  pla = function(la,t1,t3,l,v){
    D = diag(l)
    N = length(t3)
    la = as.vector(la)
    tmp = diag(1,length(l))
    tmp[upper.tri(tmp)] = la
    LA = t(tmp) 
    
    S = solve(LA) %*% D %*% t(solve(LA))
    Sl = v*rep(l[2:length(l)],1:(length(l)-1))
    
    tt = t(cbind(t1,t3))
    
    (-0.5*sum(la^2/Sl)
      -(N/2)*log(det(S))-0.5*psych::tr(t(tt)%*%solve(S)%*%tt))
  }
  
  pS2 = function(S2,t2,v0,d0){
    n = length(t2)
    (-v0-1)*log(S2)-d0/(S2)-(n/2)*log(S2)-0.5*sum(t2^2)/S2
  }
  
  Sigmat = function(l, la){
    D = diag(l)
    la = as.vector(la)
    tmp = diag(1,length(l))
    tmp[upper.tri(tmp)] = la
    LA = t(tmp) 
    solve(LA) %*% D %*% t(solve(LA))
  }
  
  # Estimate
  eap = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,node,c]/nchain 
    }
    out = apply(tmp, 1, mean)
    return(out)
  } 
  
  eap1d = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[node,c]/nchain 
    }
    out = mean(tmp)
    return(out)
  } 
  
  MCMC = function(x){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,,c]/nchain 
    }
    return(tmp)
  } 
  
  # R Square
  RR = function(x){
    items = dim(x)[1]
    len = dim(x)[2]
    nchain = dim(x)[3]
    B = NULL
    W = NULL
    for(i in 1:items){
      B_tmp = NULL
      W_tmp = NULL
      for(c in 1:nchain){
        B_tmp = c(B_tmp, mean(x[i,,c]))
        W_tmp = c(W_tmp, var(x[i,,c]))
      }
      B = c(B, len*var(B_tmp))
      W = c(W, mean(W_tmp))
    }
    Var = (len-1)*W/len + B/len
    sqrt(Var/W)
  }
  
  RR_m = function(x){
    items = dim(x)[1]
    Dim = dim(x)[2]
    len = dim(x)[3]
    nchain = dim(x)[4]
    B = NULL
    W = NULL
    for(d in 1:Dim){
      B_tmp2 = NULL
      W_tmp2 = NULL
      for(i in 1:items){
        B_tmp = NULL
        W_tmp = NULL
        for(c in 1:nchain){
          B_tmp = c(B_tmp, mean(x[i,d,,c]))
          W_tmp = c(W_tmp, var(x[i,d,,c]))
        }
        B_tmp2 = c(B_tmp2, len*var(B_tmp))
        W_tmp2 = c(W_tmp2, mean(W_tmp))
      }
      B = cbind(B,B_tmp2)
      W = cbind(W,W_tmp2)
    }
    Var = (len-1)*W/len + B/len
    Out = sqrt(Var/W)
    colnames(Out) = paste('Dim',1:Dim)
    
    return(Out)
  }
  
  r = 0
  N = nrow(y)
  K = J
  items = ncol(y)
  
  ## Initial value
  
  if(ReTrain){
    Est = mod$Estimate
    delta = array(Est$Delta, dim = c(items, iters, nchain))
    a1 = array(Est$Alpha, dim = c(items, iters, nchain))
    a2 = array(Est$Gamma, dim = c(items, iters, nchain))
    w1 = array(Est$Omega, dim = c(items, iters, nchain))
    
    tau =array(Est$Tau,dim = c(J-1, iters, nchain))
    
    theta1 = array(Est$ThetaR,dim = c(N, length(pois), iters, nchain))
    theta2 = array(Est$ThetaD,dim = c(N, iters, nchain))
    theta3 = array(Est$ThetaT,dim = c(N, iters, nchain))
    
    M = BNSP::chol(Est$Sigma11)
    
    A = array(t(M$L)[upper.tri(t(M$L))],dim=c(sum(1:length(pois)),iters,nchain))
    L = array(diag(M$D),dim=c(length(pois)+1,iters,nchain))
    Sa1 = rep(1,nchain)
    Sa2 = rep(1,nchain)
    Sd = rep(1,nchain)
    Sw1 = rep(1,nchain)
    St = rep(1,nchain)
    Stt2 = matrix(Est$Sigma22,iters,nchain)
    
    P_rd = array(0,dim=c(N,items,J,nchain))
    P_t = array(0,dim=c(N,items,J-1,nchain))
    P_rdt = array(0,dim=c(N,items,J,nchain))
  }else{
    ### Chain
    delta = array(scale(apply(y,2,sum)),dim = c(items, iters, nchain))
    a1 = array(1,dim = c(items, iters, nchain))
    a2 = array(0,dim = c(items, iters, nchain))
    w1 = array(0,dim = c(items, iters, nchain))
    
    tau =array(seq(-2,2,length.out=J-1),dim = c(J-1, iters, nchain))
    
    theta1 = array(scale(apply(y,1,sum)),dim = c(N, length(pois), iters, nchain))
    if(length(pois)>1){
      for(l in 1:length(pois)){
        theta1[,l,,] = scale(apply(y[,t_pois==l],1,sum))
      }
    }
    theta2 = array(0,dim = c(N, iters, nchain))
    theta3 = array(0,dim = c(N, iters, nchain))
    
    A = array(rep(0,sum(1:length(pois))),dim=c(sum(1:length(pois)),iters,nchain))
    L = array(rep(1,(length(pois)+1)),dim=c(length(pois)+1,iters,nchain))
    Sa1 = rep(1,nchain)
    Sa2 = rep(1,nchain)
    Sd = rep(1,nchain)
    Sw1 = rep(1,nchain)
    St = rep(1,nchain)
    Stt2 = matrix(1,iters,nchain)

    P_rd = array(0,dim=c(N,items,J,nchain))
    P_t = array(0,dim=c(N,items,J-1,nchain))
    P_rdt = array(0,dim=c(N,items,J,nchain))
  }
  
  # Parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(nchain, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  cat('Estimation process begins: \n')
  
  ## Estimation
  tmp = foreach(c = 1:nchain, .combine = list, .multicombine=TRUE) %dopar%{
    
    # All copyright / intellectual property rights of the Software belong to the Licensor, and they are protected by Intellectual.
    # If you use this software in your work, you shall add this work to your reference list.
    # Leng, CH., Huang, HY. & Yao, G. A Social Desirability Item Response Theory Model: Retrieve–Deceive–Transfer. Psychometrika 85, 56–74 (2020). https://doi.org/10.1007/s11336-019-09689-y
    
    for(i in 1:items){
      P_rd[,i,,c] = Prob_RD(delta[i,1,c],tau[,1,c],a1[i,1,c],a2[i,1,c],theta1[,t_pois[i],1,c],theta2[,1,c])
      P_t[,i,,c] = Prob_T(tau[,1,c],w1[i,1,c],theta3[,1,c])
      P_rdt[,i,,c] = Prob_RDT(P_rd[,i,,c],P_t[,i,,c],N,J)
    }
    
    for(iter in 1:(iters-1)){
      for(i in 1:items){
        ### a1
        proposal = rlnorm(1,log(a1[i,iter,c]),0.01)
        tmp_P_rd = Prob_RD(delta[i,iter,c],tau[,iter,c],proposal,a2[i,iter,c],theta1[,t_pois[i],iter,c],theta2[,iter,c])
        tmp_P_rdt = Prob_RDT(tmp_P_rd,P_t[,i,,c],N,J)
        R = exp(log(proposal)+
                  pa1(proposal,Sa1[c],y[,i],tmp_P_rdt)
                -log(a1[i,iter,c])-
                  pa1(a1[i,iter,c],Sa1[c],y[,i],P_rdt[,i,,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        a1[i,iter+1,c] = ifelse(accept,proposal,a1[i,iter,c])
        if(accept){
          P_rd[,i,,c] = tmp_P_rd
          P_rdt[,i,,c] = tmp_P_rdt
        }
        
        ### a2
        proposal = rnorm(1,(a2[i,iter,c]),1)
        tmp_P_rd = Prob_RD(delta[i,iter,c],tau[,iter,c],a1[i,iter+1,c],proposal,theta1[,t_pois[i],iter,c],theta2[,iter,c])
        tmp_P_rdt = Prob_RDT(tmp_P_rd,P_t[,i,,c],N,J)
        R = exp(pa2(proposal,Sa2[c],y[,i],tmp_P_rdt)-
                  pa2(a2[i,iter,c],Sa2[c],y[,i],P_rdt[,i,,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        a2[i,iter+1,c] = ifelse(accept,proposal,a2[i,iter,c])
        if(accept){
          P_rd[,i,,c] = tmp_P_rd
          P_rdt[,i,,c] = tmp_P_rdt
        }
        
        ### delta
        proposal = rnorm(1,(delta[i,iter,c]),0.5)
        tmp_P_rd = Prob_RD(proposal,tau[,iter,c],a1[i,iter+1,c],a2[i,iter+1,c],theta1[,t_pois[i],iter,c],theta2[,iter,c])
        tmp_P_rdt = Prob_RDT(tmp_P_rd,P_t[,i,,c],N,J)
        R = exp(pdelta(proposal,Sd[c],y[,i],tmp_P_rdt)-
                  pdelta(delta[i,iter,c],Sd[c],y[,i],P_rdt[,i,,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        delta[i,iter+1,c] = ifelse(accept,proposal,delta[i,iter,c])
        if(accept){
          P_rd[,i,,c] = tmp_P_rd
          P_rdt[,i,,c] = tmp_P_rdt
        }
        
        ### w1
        proposal = rnorm(1,(w1[i,iter,c]),1)
        tmp_P_t = Prob_T(tau[,iter,c],proposal,theta3[,iter,c])
        tmp_P_rdt = Prob_RDT(P_rd[,i,,c],tmp_P_t,N,J)
        R = exp(pw1(proposal,Sd[c],y[,i],tmp_P_rdt)-
                  pw1(w1[i,iter,c],Sd[c],y[,i],P_rdt[,i,,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        w1[i,iter+1,c] = ifelse(accept,proposal,w1[i,iter,c])
        if(accept){
          P_t[,i,,c] = tmp_P_t
          P_rdt[,i,,c] = tmp_P_rdt
        }
      }
      
      {
        ### tau
        proposal = mvtnorm::rmvnorm(1,tau[1:(J-2),iter,c],diag(0.05,J-2))
        proposal = (c(proposal,-sum(proposal)))
        tmp_P_rd = array(0,dim=c(N,items,J))
        tmp_P_t = array(0,dim=c(N,items,J-1))
        tmp_P_rdt = array(0,dim=c(N,items,J))
        for(i in 1:items){
          tmp_P_rd[,i,] = Prob_RD(delta[i,iter+1,c],proposal,a1[i,iter+1,c],a2[i,iter+1,c],theta1[,t_pois[i],iter,c],theta2[,iter,c])
          tmp_P_t[,i,] = Prob_T(proposal,w1[i,iter+1,c],theta3[,iter,c])
          tmp_P_rdt[,i,] = Prob_RDT(tmp_P_rd[,i,],tmp_P_t[,i,],N,J)
        }
        R = exp(ptau(proposal,St[c],y,tmp_P_rdt,J)-
                  ptau(tau[,iter,c],St[c],y,P_rdt[,,,c],J))
        R[is.na(R)] = 0
        accept = rbinom(1,1,min(1,R))
        tau[,iter+1,c] = ifelse(rep(accept,J-1),proposal,tau[,iter,c])
        if(accept){
          P_rd[,,,c] = tmp_P_rd
          P_t[,,,c] = tmp_P_t
          P_rdt[,,,c] = tmp_P_rdt
        }
      }
      
      {
        ### Vairance of item parameters
        Sd[c] = ss(delta[,iter+1,c],3,10)
        Sa1[c] = ss2((a1[,iter+1,c]),3,8)
        Sa2[c] = ss((a2[,iter+1,c]),3,10)
        Sw1[c] = ss((w1[,iter+1,c]),3,10)
        St[c] = ss(tau[,iter+1,c],3,10)
        
        ### Variance of theta2
        b3 = rgamma(1, shape=5.5, rate = 2*Stt2[iter,c]/(2+Stt2[iter,c]))
        proposal = rlnorm(1,log(Stt2[iter,c]),0.5)
        R = exp(pS2(proposal,theta2[,iter,c],3,b3)+log(proposal)-
                  (pS2(Stt2[iter,c],theta2[,iter,c],3,b3)+log(Stt2[iter,c])))
        R[is.na(R)] = 0
        accept = rbinom(1,1,min(1,R))
        Stt2[iter+1,c] = ifelse(accept,proposal,Stt2[iter,c]) 
        
        ### L
        proposal2 = L[,iter,c]
        for(l in 1:(length(pois)+1)){
          ### Lambda1, 3
          b3 = rgamma(1, shape=5.5, rate = 2*L[l,iter,c]/(2+L[l,iter,c]))
          temp = proposal2
          proposal = rlnorm(1,log(L[l,iter,c]),0.5)
          proposal2[l] = proposal
          R = exp(plambda(proposal2,theta1[,,iter,c],theta3[,iter,c],A[,iter,c],3,b3)+log(proposal2[l])-
                    (plambda(temp,theta1[,,iter,c],theta3[,iter,c],A[,iter,c],3,b3)+log(temp[l])))
          R[is.na(R)] = 0
          accept = rbinom(1,1,min(1,R))
          proposal2[l] = ifelse(accept,proposal,L[l,iter,c]) 
        }
        L[,iter+1,c] = proposal2
        
        ### A
        proposal2 = A[,iter,c]
        for(j in 1:sum(1:length(pois))){
          temp = proposal2
          proposal = rnorm(1,A[j,iter,c],0.2)
          proposal2[j] = proposal
          R = exp(pla(proposal2,theta1[,,iter,c],theta3[,iter,c],L[,iter+1,c],5)-
                    pla(temp,theta1[,,iter,c],theta3[,iter,c],L[,iter+1,c],5))
          if(is.na(R)==T){R=0}
          accept = rbinom(1,1,min(1,R))
          proposal2[j] = ifelse(accept, proposal, A[j,iter,c])
        }
        A[,iter+1,c] = proposal2
        
        Sigma = Sigmat(L[,iter+1,c],A[,iter+1,c])
      }
      
      for(n in 1:N){
        ### theta1
        proposal2 = theta1[n,,iter,c]
        for(l in 1:length(pois)){
          tmp = proposal2
          proposal = rnorm(1,theta1[n,l,iter,c],0.05)
          proposal2[l] = proposal 
          tmp_P_rd = Prob_RD_theta(delta[t_pois==l,iter+1,c],tau[,iter+1,c],a1[t_pois==l,iter+1,c],a2[t_pois==l,iter+1,c],proposal,theta2[n,iter,c],t_pois,R=T)
          tmp_P_rdt = Prob_RDT_theta(tmp_P_rd,P_t[n,t_pois==l,,c],pois[l],J)
          R = exp(ptheta1(proposal2,Sigma,y[n,t_pois==l],tmp_P_rdt,theta3[n,iter,c])-
                    ptheta1(tmp,Sigma,y[n,t_pois==l],P_rdt[n,t_pois==l,,c],theta3[n,iter,c]))
          if(is.na(R)==T){R=0}
          accept = rbinom(1,1,min(1,R))
          proposal2[l] = ifelse(accept, proposal, theta1[n,l,iter,c])
          if(accept){
            P_rd[n,t_pois==l,,c] = tmp_P_rd
            P_rdt[n,t_pois==l,,c] = tmp_P_rdt
          }
        }
        theta1[n,,iter+1,c] = proposal2
        
        ### theta2
        proposal = rnorm(1,(theta2[n,iter,c]),1)
        tmp_P_rd = Prob_RD_theta(delta[,iter+1,c],tau[,iter+1,c],a1[,iter+1,c],a2[,iter+1,c],theta1[n,,iter+1,c],proposal,t_pois,R=F)
        tmp_P_rdt = Prob_RDT_theta(tmp_P_rd,P_t[n,,,c],items,J)
        R = exp(ptheta2(proposal,Stt2[iter+1,c],y[n,],tmp_P_rdt)-
                  ptheta2(theta2[n,iter,c],Stt2[iter+1,c],y[n,],P_rdt[n,,,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        theta2[n,iter+1,c] = ifelse(accept,proposal,theta2[n,iter,c])
        if(accept){
          P_rd[n,,,c] = tmp_P_rd
          P_rdt[n,,,c] = tmp_P_rdt
        }
        
        ### theta3
        proposal = rnorm(1,theta3[n,iter,c],1)
        tmp_P_t = Prob_T_theta(tau[,iter+1,c],w1[,iter+1,c],proposal)
        tmp_P_rdt = Prob_RDT_theta(P_rd[n,,,c],tmp_P_t,items,J)
        R = exp(ptheta3(proposal,Sigma,y[n,],tmp_P_rdt,theta1[n,,iter+1,c])-
                  ptheta3(theta3[n,iter,c],Sigma,y[n,],P_rdt[n,,,c],theta1[n,,iter+1,c]))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        theta3[n,iter+1,c] = ifelse(accept,proposal,theta3[n,iter,c])
        if(accept){
          P_t[n,,,c] = tmp_P_t
          P_rdt[n,,,c] = tmp_P_rdt
        }
      }
      
    }
    out = (list(tau = tau[,,c],
               a1 = a1[,,c],
               a2 = a2[,,c],
               delta = delta[,,c],
               w1 = w1[,,c],
               theta1 = theta1[,,,c],
               theta2 = theta2[,,c],
               theta3 = theta3[,,c],
               Stt2 = Stt2[,c],
               L = L[,,c],
               A = A[,,c]))
            
    return(out)
    
  }
  
  stopCluster(my.cluster)
  
  cat("Estimation process ends. \n")
  cat("Outputing the estimation results. \n")
  
  for(c in 1:nchain){
    a1[,,c] = tmp[[c]]$a1
    a2[,,c] = tmp[[c]]$a2
    delta[,,c] = tmp[[c]]$delta
    tau[,,c] = tmp[[c]]$tau
    w1[,,c] = tmp[[c]]$w1
    theta1[,,,c] = tmp[[c]]$theta1
    theta2[,,c] = tmp[[c]]$theta2
    theta3[,,c] = tmp[[c]]$theta3
    Stt2[,c] = tmp[[c]]$Stt2
    L[,,c] = tmp[[c]]$L
    A[,,c] = tmp[[c]]$A
  }
  
  # EAP
  node = seq(burnin+thin,iters,by=thin)

  # Estimate
  a1_out = eap(a1, node); 
  a2_out = eap(a2, node); 
  delta_out = eap(delta, node); 
  w1_out = eap(w1, node); 
  tau_out = eap(tau, node)
  
  theta1_out = matrix(0,N,length(pois))
  for(l in 1:length(pois)){
    theta1_out[,l] = eap(theta1[,l,,], node); 
  }
  theta2_out = eap(theta2, node); 
  theta3_out = eap(theta3, node); 
  
  Stt2_out = eap1d(Stt2, node)  
  
  A_out = rep(0, sum(1:length(pois)))
  for(l in 1:length(A_out)){
    A_out[l] = eap1d(A[l,,], node)
  }
  L_out = eap(L, node)
  Sigma = Sigmat(L_out,A_out)
  
  # R Square
  a1_RR = RR(a1[,node,])
  a2_RR = RR(a2[,node,])
  delta_RR = RR(delta[,node,])
  w1_RR = RR(w1[,node,])
  tau_RR = RR(tau[,node,])
  if(length(pois)==1){
    theta1_RR = RR(theta1[,,node,])
  }else{
    theta1_RR = RR_m(theta1[,,node,])
  }
  theta2_RR = RR(theta2[,node,])
  theta3_RR = RR(theta3[,node,])
  
  # Output
  out = list(
    Estimate = list(
      Alpha = a1_out,
      Gamma = a2_out,
      Delta = delta_out,
      Omega = w1_out,
      Tau = tau_out,
      ThetaR = theta1_out,
      ThetaD = theta2_out,
      ThetaT = theta3_out,
      Sigma11 = Sigma,
      Sigma22 = Stt2_out
    ),
    MCMC = list(
      Alpha = a1,
      Gamma = a2,
      Delta = delta,
      Omega = w1,
      Tau = tau,
      ThetaR = theta1,
      ThetaD = theta2,
      ThetaT = theta3,
      Sigma11 = Sigma,
      Sigma22 = Stt2
    ),
    RSquare = list(
      Alpha = a1_RR,
      Gamma = a2_RR,
      Delta = delta_RR,
      Omega = w1_RR,
      Tau = tau_RR,
      ThetaR = theta1_RR,
      ThetaD = theta2_RR,
      ThetaT = theta3_RR
    ),
    Node = node
  )
  
  cat('Done. \n')
  return(out)
}

plotRDTMCMC = function(mod, x, i, node=mod$Node){
  X = mod$MCMC[[x]]
  ylim = c(min(X[i,node,]),max(X[i,node,]))
  
  for(c in 1:(dim(X)[3])){
    plot(node,X[i,node,c],'l',col=c,ylim=ylim,ylab=""); par(new=T)
  }
  par(new=F)
  title(ylab=paste(x,i,sep='_'))
}

plotRDTMCMC_M = function(mod, x, d, i, node=mod$Node){
  X = mod$MCMC[[x]]
  ylim = c(min(X[i,d,node,]),max(X[i,d,node,]))
  
  for(c in 1:(dim(X)[4])){
    plot(node,X[i,d,node,c],'l',col=c,ylim=ylim,ylab=""); par(new=T)
  }
  par(new=F)
  title(ylab=paste(x,i,sep='_'))
}

# DIC
DIC = function(y, mod, nchain, t_pois, iters=15000, burnin=5000, thin=10){
  # ppDIC
  logL = function(yy,delta,tau,a1,ga,w1,theta1,theta2,theta3){
    
    N = length(theta1)
    J = length(tau)+1
    
    # RD
    trd = (a1*theta1+exp(theta2+ga)-delta)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_trd = matrix(trd, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp(m_trd - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob_rd = prob / c_Prob
    
    # T
    tt = matrix(theta3+w1,N,J-1,byrow=F)
    mw2 = matrix(tau,N,J-1,byrow=T)
    door = exp(-(tt-mw2)*(tt>(mw2)))
    door[which(is.na(door)==T)] = 1
    
    # RDT
    d = (door[,1:(J-1)] - cbind(0,door[,1:(J-2)])) # 11, 12, ..., 1J-1
    d2 = cbind((matrix(1,N,J-1) - door[,(J-1):1]),1) # 1J
    Prob_rdt = NULL
    for(j in 1:(J-1)){
      Prob_t = d[,j:1]
      tmp =  matrix((Prob_rd[,1:j] * Prob_t),N,j,byrow=F)
      Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    }
    tmp = matrix((Prob_rd * d2),N,J,byrow=F)
    Prob_rdt = cbind(Prob_rdt, apply(tmp, 1, sum))
    
    # Y
    m_y = matrix(yy, nrow=N, ncol=J, byrow=F) == matrix(0:(J-1), nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob_rdt * m_y,1,sum)
    
    return(sum(log(sum_Prob)))
  }
  
  DIC = LogL = NULL
  {
    S = mod$Node
    N = nrow(y)
    items = ncol(y)
  }
  
  d = tau = a = g = w = t1 = t2 = t3 = 0
  for(c in 1:nchain){
    d = d + mod$MCMC$Delta[, mod$Node,c] / nchain
    tau = tau + mod$MCMC$Tau[, mod$Node,c] / nchain
    a = a + mod$MCMC$Alpha[, mod$Node,c] / nchain
    g = g + mod$MCMC$Gamma[, mod$Node,c] / nchain
    w = w + mod$MCMC$Omega[, mod$Node,c] / nchain
    t1 = t1 + mod$MCMC$ThetaR[,, mod$Node,c] / nchain
    t2 = t2 + mod$MCMC$ThetaD[, mod$Node,c] / nchain
    t3 = t3 + mod$MCMC$ThetaT[, mod$Node,c] / nchain
  }
  t1 = array(t1, dim=c(N,length(unique(t_pois)),length(S)))
  
  Est = mod$Estimate
  
  ppdic = 0
  for(s in 1:length(S)){
    for(i in 1:items){
      ppdic = ppdic + logL(y[,i],d[i,s],tau[,s],a[i,s],g[i,s],w[i,s],t1[,t_pois[i],s],t2[,s],t3[,s])
    }
  }
  ppdic = ppdic / length(S)
  
  logl = 0
  for(i in 1:items){
    logl = logl + logL(y[,i],Est$Delta[i],Est$Tau,Est$Alpha[i],Est$Gamma[i],Est$Omega[i],
                       Est$ThetaR[,t_pois[i]],Est$ThetaD,Est$ThetaT)
  }
  
  pdic = 2*(logl-ppdic)
  DIC = c(DIC, -2*logl+2*pdic)
  LogL = c(LogL, -2*logl)
  
  out = list(LogL = LogL, DIC=DIC)
  return(out)
}

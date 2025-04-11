
library(pracma)
library(MASS)
library(mvtnorm)
library(bbmle)
library(mle.tools)
library(rworldmap)
library(geosphere)
library(plotly)
library(numDeriv)



c_H<-function(u,v,h){
  return(0.5*(u^(2*h)+v^(2*h)-abs(u-v)^(2*h)))
 }



cov_OUH_zeta_2<-function(T,n,b,h) {
    D<-T/n
    index<-0
    G<-function(u,v){return(exp(-b*(u+v))*abs(v-u+D*index)^(2*h))}
    H<-function(u){return(exp(-b*u)*(D*(index+1)-u)^(2*h))}
    int_2<-rep(0,n)
    int_1<-rep(0,n)
    for(i in 1:n){
    int_2[i]<-integral2(G,0,D,0,D,reltol = 1e-6)$Q
    int_1[i]<-integrate(H,0,D)$val
    index<-i 
    }
    K<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:i){
              K[i,j]<-0.5*((1-exp(-b*D))/b)*(int_1[i]+int_1[j])-0.5*int_2[i-j+1]
              K[j,i]<-K[i,j]
        }
    }
    return(K)
 }






cov_OUH_zeta<-function(T,n,r,b_1,b_2,h_1,h_2){
	D<-T/n
	index<-0
    G<-function(u,v){return(exp(-b_1*u-b_2*v)*abs(u-v+D*index)^(h_1+h_2))}
    H_1<-function(u){return(exp(-b_1*u)*(D*(index+1)-u)^(h_1+h_2))}
    H_2<-function(v){return(exp(-b_2*v)*(D*(index+1)-v)^(h_1+h_2))}
    int_2<-rep(0,2*n-1)
    #int_2<-rep(0,n)
    int_1_1<-rep(0,n)
    int_1_2<-rep(0,n)
    for(i in 1:(2*n-1)){
     index<--(n-1)+i-1	
     int_2[i]<-integral2(G,0,D,0,D,reltol = 1e-6)$Q 
    }
    index<-0
    for(i in 1:n){
    int_1_1[i]<-integrate(H_1,0,D)$val
    int_1_2[i]<-integrate(H_2,0,D)$val
    index<-i
    }
    K<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:n){
              
              K[i,j]<-0.5*r*(((1-exp(-b_2*D))/b_2)*int_1_1[i]
              +((1-exp(-b_1*D))/b_1)*int_1_2[j]
              -int_2[j-i+n])
              
        }
    }
    return(K)
}



cov_OUH<-function(T,n,r_12,r_13,r_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3) {
    D<-T/n
    K<-matrix(rep(0,(9*n*n)),3*n,3*n)
    S_B_1<-matrix(rep(0,n*n),n,n)
    S_B_2<-matrix(rep(0,n*n),n,n)
    S_B_3<-matrix(rep(0,n*n),n,n)
    for(i in 1:n){
       for(j in 1:i){
        	 S_B_1[i,j]<-exp((j-i)*b_1*D)
             S_B_2[i,j]<-exp((j-i)*b_2*D)
             S_B_3[i,j]<-exp((j-i)*b_3*D)
             } 
       }
    K[1:n,1:n]<-(s_1*S_B_1%*%cov_OUH_zeta_2(T,n,b_1,h_1)%*%t(s_1*S_B_1))
    K[(n+1):(2*n),(n+1):(2*n)]<-(s_2*S_B_2%*%cov_OUH_zeta_2(T,n,b_2,h_2)%*%t(s_2*S_B_2))
    K[(2*n+1):(3*n),(2*n+1):(3*n)]<-(s_3*S_B_3%*%cov_OUH_zeta_2(T,n,b_3,h_3)%*%t(s_3*S_B_3))
    K[1:n,(n+1):(2*n)]<-(s_1*S_B_1%*%cov_OUH_zeta(T,n,r_12,b_1,b_2,h_1,h_2)%*%t(s_2*S_B_2))
    K[1:n,(2*n+1):(3*n)]<-(s_1*S_B_1%*%cov_OUH_zeta(T,n,r_13,b_1,b_3,h_1,h_3)%*%t(s_3*S_B_3))
    K[(n+1):(2*n),(2*n+1):(3*n)]<-(s_2*S_B_2%*%cov_OUH_zeta(T,n,r_23,b_2,b_3,h_2,h_3)%*%t(s_3*S_B_3))
    K[(2*n+1):(3*n),(n+1):(2*n)]<-t(K[(n+1):(2*n),(2*n+1):(3*n)])
    K[(n+1):(2*n),1:n]<-t(K[1:n,(n+1):(2*n)])
    K[(2*n+1):(3*n),1:n]<-t(K[1:n,(2*n+1):(3*n)])
    return(K)
 }







cov_OUH_zeta_3<-function(T,n,b,h) {
    index<-0
    G<-function(u,v){return(exp(-b*(u+v))*abs(v-u+index)^(2*h))}
    H<-function(u){return(exp(-b*u)*((index+1)-u)^(2*h))}
    int_2<-rep(0,n)
    int_1<-rep(0,n)
    for(i in 1:n){
    int_2[i]<-integral2(G,0,1,0,1,reltol = 1e-6)$Q
    int_1[i]<-integrate(H,0,1)$val
    index<-i 
    }
    K<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:i){
              K[i,j]<-0.5*((1-exp(-b))/b)*(int_1[i]+int_1[j])-0.5*int_2[i-j+1]
              K[j,i]<-K[i,j]
        }
    }

    return(K)
 }




simulation_vel<-function(T,s,b,h,mu,index,sim){
    n<-index[length(index)]
    D<-T/n
    K_2<-cov_OUH_zeta_2(T,n,b,h)    
    z_2<-rep(0,n)
    S_B<-matrix(rep(0,n*n),n,n)
    K_1<-matrix(rep(0,n*n),n,n)
    for(i in 1:n){
      for(j in 1:i){
          K_1[i,j]<-s*s*(c_H(j*D,i*D,h)-exp(-b*D)*(c_H(j*D,(i-1)*D,h)+c_H(i*D,(j-1)*D,h))+exp(-2*b*D)*c_H((i-1)*D,(j-1)*D,h))
          K_1[j,i]<-K_1[i,j]
      }
    }
   index_aux<-0
   H<-function(u){return(exp(-b*u)*((index_aux+1)*D-u)^(2*h))}
   int_2<-rep(0,2*n)
   int_1<-rep(0,n)
   for(i in 1:n){
    int_1[i]<-integrate(H,0,D)$val
    index_aux<-i 
   }
   G<-function(u){return(exp(b*(u-D))*abs(index_aux*D-u)^(2*h))}
   for(i in (-n+1):n){
    index_aux<-i
    int_2[n+i]<-integrate(G,0,D)$val 
   }
   K_21<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:n){
              K_21[i,j]<-( -0.5*s*s*b*(  ( (j*D)^(2*h)       )*( (1-exp(-b*D))/b )+int_1[i]-int_2[j-i+1+n]          )
                           +0.5*s*s*b*exp(-1*b*D)*(  ( ( (j-1)*D )^(2*h) )*( (1-exp(-b*D))/b )+int_1[i]-int_2[j-i+n]           )     )
       }
    }
 
  for(i in 1:n) {
        for(j in 1:i){
             S_B[i,j]<-exp((j-i)*b*D)
       }
    }
  M<-(s*S_B%*%K_2%*%t(s*S_B))
  M_2<-M[index,index]
  M_1<-M[-index,-index]
  M_21<-M[index,-index]
  M_12<-M[-index,index]
  K_12<-matrix(rep(0,n*n),n,n)  
  K_12<-t(K_21)   
  vel_aux<-rep(0,n)
  vel<-rep(0,n)
  mu_com<-rep(0,n)
  #---------------Simulation----
  mu_vel_sim<-matrix(rep(0,n*2*sim),2*sim,n)
  for(i in 1:sim){
  mu_vel_sim[2*i-1,-index]<-mvrnorm( 1,rep(mu[1],n-length(index))+M_12%*%solve(M_2)%*%as.vector(mu[-1]-mu[1]), M_1-M_12%*%solve(M_2)%*%M_21)
  mu_vel_sim[2*i-1,index]<-mu[-1]         #--Es sin el punto inicial, para imprimir debemos poner el punto inicial 
  z_2<--b*solve(S_B)%*%(as.vector(mu_vel_sim[2*i-1,]-mu[1])) #---con los cambios creo que es sin exp(-b*D)
  z_1<-mvrnorm( 1,K_12%*%solve(s*s*b*b*K_2)%*%as.vector(z_2), K_1-K_12%*%solve(s*s*b*b*K_2)%*%K_21)
  vel_aux<-z_1+z_2
  mu_vel_sim[2*i,1]<-vel_aux[1]
  for(j in 2:n){
    mu_vel_sim[2*i,j]<-vel_aux[j]+exp(-b*D)*mu_vel_sim[2*i,j-1]
  }

 }
  return(mu_vel_sim)
}





zeta_2<-function(D,n,s_1,s_2,b_1,b_2,h_1,h_2,rho){  #i=1,j=2. cov y vs w
   index_aux<-0
   H<-function(u){return(exp(-b_1*u)*((index_aux+1)*D-u)^(h_1+h_2))}
   int_2<-rep(0,2*n)
   int_1<-rep(0,n)
   for(i in 1:n){
    int_1[i]<-integrate(H,0,D)$val
    index_aux<-i 
   }
   index_aux<-0
   G<-function(u){return(exp(-b_1*u)*abs(index_aux*D+u)^(h_1+h_2))}
   int_2<-rep(0,2*n)
   for(i in 1:(2*n)){
     index_aux<--n+i-1	
     int_2[i]<-integrate(G,0,D)$val
    }
   K_21<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:n){
              K_21[i,j]<- 0.5*rho*(  ( (j*D)^(h_1+h_2)       )*( (1-exp(-b_1*D))/b_1 )+int_1[i]-int_2[j-i+n+1]          )
                         
       }
    }
 return(K_21)
 }



K<-function(s,t,h_1,h_2,rho)
{
	return(rho*0.5*(t^(h_1+h_2)+s^(h_1+h_2)-abs(t-s)^(h_1+h_2)))
	
}



simulation_vel_3D<-function(T,r,s,b,h,mu_dat,index,sim){
    n<-index[length(index)]
    D<-T/n
    K_W<-matrix(rep(0,9*n*n),3*n,3*n)
    for(i in 1:n){
      for(j in 1:n){
          K_W[i,j]<-K(i*D,j*D,h[1],h[1],1)  #W_11,W_11  
          
          K_W[n+i,n+j]<-K(i*D,j*D,h[2],h[2],1)  #W_21,W_21
          
          K_W[2*n+i,2*n+j]<-K(i*D,j*D,h[3],h[3],1)   #W_31,W_31
          
          K_W[i,n+j]<-K(i*D,j*D,h[1],h[2],r[1])  #W_11,W_21
          
          K_W[n+i,j]<-K(i*D,j*D,h[1],h[2],r[1]) #W_21,W_11  
                
          K_W[i,2*n+j]<-K(i*D,j*D,h[1],h[3],r[2])  #W_11,W_31   
               
          K_W[2*n+i,j]<-K(i*D,j*D,h[1],h[3],r[2]) #W_31,W_11
          
          K_W[n+i,2*n+j]<-K(i*D,j*D,h[2],h[3],r[3])  #W_21,W_31
          
          K_W[2*n+i,n+j]<-K(i*D,j*D,h[2],h[3],r[3])  #W_31,W_21
      }
    }
    
    
    K_mu_W<-matrix(rep(0,9*n*n),3*n,3*n)
    
    
    S_B_1<-matrix(rep(0,n*n),n,n)
    S_B_2<-matrix(rep(0,n*n),n,n)
    S_B_3<-matrix(rep(0,n*n),n,n)
    for(i in 1:n){
       for(j in 1:i){
        	 S_B_1[i,j]<-exp((j-i)*b[1]*D)
             S_B_2[i,j]<-exp((j-i)*b[2]*D)
             S_B_3[i,j]<-exp((j-i)*b[3]*D)
             } 
       }
    
    
    
    #---zeta_2 da y vs w
    
    
    K_mu_W[c((1):(n)),c(1:n)]<-s[1]*S_B_1%*%zeta_2(D,n,s[1],s[1],b[1],b[1],h[1],h[1],1) #z_1,2,z_1,1
    K_mu_W[c((1):(n)),c((n+1):(2*n))]<-s[1]*S_B_1%*%zeta_2(D,n,s[1],s[2],b[1],b[2],h[1],h[2],r[1])   #z_1,2,z_2,1
    K_mu_W[c((1):(n)),c((2*n+1):(3*n))]<-s[1]*S_B_1%*%zeta_2(D,n,s[1],s[3],b[1],b[3],h[1],h[3],r[2])   #z_1,2,z_3,1     
    K_mu_W[c((n+1):(2*n)),c(1:n)]<-s[2]*S_B_2%*%zeta_2(D,n,s[2],s[1],b[2],b[1],h[2],h[1],r[1]) #z_2,2,z_1,1
    K_mu_W[c((n+1):(2*n)),c((n+1):(2*n))]<-s[2]*S_B_2%*%zeta_2(D,n,s[2],s[2],b[2],b[2],h[2],h[2],1)   #z_2,2,z_2,1
    K_mu_W[c((n+1):(2*n)),c((2*n+1):(3*n))]<-s[2]*S_B_2%*%zeta_2(D,n,s[2],s[3],b[2],b[3],h[2],h[3],r[3])   #z_2,2,z_3,1   
    K_mu_W[c((2*n+1):(3*n)),c(1:n)]<-s[3]*S_B_3%*%zeta_2(D,n,s[3],s[1],b[3],b[1],h[3],h[1],r[2]) #z_3,2,z_1,1
    K_mu_W[c((2*n+1):(3*n)),c((n+1):(2*n))]<-s[3]*S_B_3%*%zeta_2(D,n,s[3],s[2],b[3],b[2],h[3],h[2],r[3])   #z_3,2,z_2,1
    K_mu_W[c((2*n+1):(3*n)),c((2*n+1):(3*n))]<-s[3]*S_B_3%*%zeta_2(D,n,s[3],s[3],b[3],b[3],h[3],h[3],1)   #z_3,2,z_3,1 
    
   K_mu<-cov_OUH(T,n,r[1],r[2],r[3],s[1],s[2],s[3],b[1],b[2],b[3],h[1],h[2],h[3]) 
   K_mu_i<-solve(K_mu) 
    
    
   z_1<-rep(0,3*n)
   z_2<-rep(0,3*n)
   vel_d<-matrix(rep(0,sim*3*n),sim,3*n)
   vel<-matrix(rep(0,sim*3*n),sim,3*n)
   mu<-rep(0,3*(n+1))
   mu[c(1,1+index,n+2,n+2+index,2*n+3,2*n+3+index)]<-mu_dat 
   K_mu_dat<-K_mu[c(index,n+index,2*n+index),c(index,n+index,2*n+index)]
   K_mu_dat_inv<-solve(K_mu[c(index,n+index,2*n+index),c(index,n+index,2*n+index)])  
   K_res_mu_dat<-K_mu[-c(index,n+index,2*n+index),-c(index,n+index,2*n+index)]
   K_res_mu_dat_and_mu_dat<-K_mu[-c(index,n+index,2*n+index),c(index,n+index,2*n+index)]
   mean_mu_aux<-c(rep(mu[1],n-length(index)),rep(mu[n+1+1],n-length(index)),rep(mu[2*(n+1)+1],n-length(index)))+K_res_mu_dat_and_mu_dat%*%K_mu_dat_inv%*%(mu_dat[-c(1,length(index)+1+1,2*(length(index)+1)+1)]-c(rep(mu[1],length(index)),rep(mu[n+1+1],length(index)),rep(mu[2*(n+1)+1],length(index))))   
   K_mu_aux<-K_res_mu_dat-K_res_mu_dat_and_mu_dat%*%K_mu_dat_inv%*%t(K_res_mu_dat_and_mu_dat) 
   K_W_given_mu<-(K_W-t(K_mu_W)%*%K_mu_i%*%K_mu_W+t(K_W-t(K_mu_W)%*%K_mu_i%*%K_mu_W))/2
  ev<-eigen(K_W_given_mu)
  Diag<-diag(sqrt(abs(ev$values)))
  V<-ev$vectors
  L<-V%*%Diag
  K_mean<-t(K_mu_W)%*%K_mu_i 
  v_initial<-rep(0,3)
  v_initial_d<-rep(0,3)
  v_initial_d[1]<-exp(-b[1]*D)*v_initial[1]
  v_initial_d[2]<-exp(-b[2]*D)*v_initial[2]
  v_initial_d[3]<-exp(-b[3]*D)*v_initial[3]
   for(q in 1:sim){ 
   mu[-c(1,1+index,n+2,n+2+index,2*n+3,2*n+3+index)]<-mvrnorm( 1,mean_mu_aux,K_mu_aux)
  mean_mu<-mu[-c(1,n+2,2*n+3)]-c(rep(mu[1],n),rep(mu[n+2],n),rep(mu[2*n+3],n))
  K_mean_c<-K_mean%*%mean_mu
  #---simular z_1 y z_2
  simulation_W<-rep(0,3*(n+1))
  simulation_W[c(((2):(n+1)),((n+3):(2*n+2)),((2*n+4):(3*n+3)))]<-L%*%rnorm( 3*n, 0,1)+ K_mean_c 
  for(i in 1:n)
      {
      	z_2[i]<-b[1]*(1-exp(-b[1]*D))*mu[1]-b[1]*(mu[i+1]-exp(-b[1]*D)*mu[i])
      	z_2[n+i]<-b[2]*(1-exp(-b[2]*D))*mu[(n+1)+1]-b[2]*(mu[(n+1)+i+1]-exp(-b[2]*D)*mu[(n+1)+i])
      	z_2[2*n+i]<-b[3]*(1-exp(-b[3]*D))*mu[2*(n+1)+1]-b[3]*(mu[2*(n+1)+i+1]-exp(-b[3]*D)*mu[2*(n+1)+i])
        z_1[i]<-s[1]*(simulation_W[i+1]-exp(-b[1]*D)*simulation_W[i])
      	z_1[n+i]<-s[2]*(simulation_W[(n+1)+i+1]-exp(-b[2]*D)*simulation_W[(n+1)+i])
      	z_1[2*n+i]<-s[3]*(simulation_W[2*(n+1)+i+1]-exp(-b[3]*D)*simulation_W[2*(n+1)+i])
      } 
  
  #------------
  vel_d[q,]<-z_1+z_2
  vel[q,1]<-vel_d[q,1]+v_initial_d[1]
  vel[q,n+1]<-vel_d[q,n+1]+v_initial_d[2]
  vel[q,2*n+1]<-vel_d[q,2*n+1]+v_initial_d[3]
  for(j in 1:(n-1)) {
            vel[q,j+1]<-vel_d[q,j+1]+exp(-b[1]*D)*vel[q,j]
            vel[q,n+j+1]<-vel_d[q,n+j+1]+exp(-b[2]*D)*vel[q,n+j]
            vel[q,2*n+j+1]<-vel_d[q,2*n+j+1]+exp(-b[3]*D)*vel[q,2*n+j]      
        }
  }  
 return(vel)
}













#------Log-likelihood in index--------#
                   



sig<-function(x){
  if(x<0){return(-1)}
  if(x>=0){return(1)}
}


degree_km<-function(lon,lat){
k<-length(lon)
change<-matrix(rep(0,2*k),2,k) 
change[1,1]<-0
change[2,1]<-0
for(i in 2:k){
change[1,i]<-sig(lon[i]-lon[1])*distm(c(lon[1], lat[1]), c(lon[i], lat[1]), fun = distHaversine)
change[2,i]<-sig(lat[i]-lat[1])*distm(c(lon[1], lat[1]), c(lon[1], lat[i]), fun = distHaversine)
}
  return(change)
}



simulation<-function(s,b,D,K,mu_ini){
    n<-length(K[1,])
    sim<-mvrnorm( 1,rep(0,n), K )
    mu<-rep(0,n+1)
    mu[1]<-mu_ini
    for(i in 1:n){
      mu[i+1]=(1-exp(-b*D))*mu[1]+s*sim[i]+exp(-b*D)*mu[i]
     }
    return(mu)
}

log_verosimilitude<-function(T,b,h,mu){
n<-length(mu)-1
D<-T/n
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*b*D)
  } 
}

M<-(S_B%*%K%*%t(S_B))
s<-(t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n
S<-s*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(n+1)],rep(mu[1],n),S,log=TRUE))
}



 
log_verosimilitude_practical<-function(T,b,h,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n                     
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*D*b)
  } 
}
M<-(S_B%*%K%*%t(S_B))
S<-M[index,index]
g<-(t(as.vector(mu[2:(k+1)]-mu[1]))%*%solve(S)%*%as.vector(mu[2:(k+1)]-mu[1]) )[1,1]/k
#g<--2*dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE)/k+2*dmvnorm(rep(0,k),rep(0,k),S,log=TRUE)/k
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),g*S,log=TRUE))
}



prediction<-function(T,s,b,h,mu,sim,step){
n<-length(mu)-1    #----sim numero de simulaciones
D<-T/n         #----step numero de pasos a predecir
K<-cov_OUH_zeta_2(T,n+step,b,h)
y<-rep(0,n)
for(i in 1:n){
    y[i]<--((1-exp(-b*D))*mu[1]-mu[i+1]+exp(-b*D)*mu[i])/s
  }

K_11<-K[1:n,1:n]
K_12<-K[1:n,(n+1):(n+step)]
K_21<-t(K[1:n,(n+1):(n+step)])
K_22<-K[(n+1):(n+step),(n+1):(n+step)]
K_11_inv<-solve(K_11)

y_predic<-rmvnorm(sim,K_21%*%K_11_inv%*%y,K_22-K_21%*%K_11_inv%*%K_12)

mu_predic<-matrix(rep(0,step*sim),sim,step)


for(i in 1:sim){
 mu_aux<-c(mu,rep(0,step))
 for(j in 1:step){
  mu_aux[n+j+1]<-(1-exp(-b*D))*mu_aux[1]+exp(-b*D)*mu_aux[n+j]+s*y_predic[i,j]
 }
 mu_predic[i,]<-mu_aux[(n+2):(n+step+1)]
}

return(mu_predic)

}

log_vero<-function(x){
-log_verosimilitude(T,b=x[1],h=x[2],mu=mu)
}





norm<-function(x){
	
	sqrt(sum(x*x))
}



coeff<-function(z_12,z_13,z_23,h_1,h_2,h_3){
gamma_11<-gamma(h_1+h_1+1)*sin(pi/2*(h_1+h_1))
gamma_12<-gamma(h_1+h_2+1)*sin(pi/2*(h_1+h_2))
gamma_13<-gamma(h_1+h_3+1)*sin(pi/2*(h_1+h_3))
gamma_23<-gamma(h_2+h_3+1)*sin(pi/2*(h_2+h_3))
gamma_22<-gamma(h_2+h_2+1)*sin(pi/2*(h_2+h_2))
gamma_33<-gamma(h_3+h_3+1)*sin(pi/2*(h_3+h_3))
r<-rep(0,3)
r[1]<-z_12*sqrt(gamma_11*gamma_22)/gamma_12
r[2]<-z_13*sqrt(gamma_11*gamma_33)/gamma_13
r[3]<-(z_23*sqrt(1-z_13^2)*sqrt(1-z_12^2)+z_12*z_13)*sqrt(gamma_22*gamma_33)/gamma_23
return(r)
}




log_likelihood_3d<-function(T,n,telemetry,z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3){
mu_1<-telemetry[1:(n+1)]
mu_2<-telemetry[(n+2):(2*n+2)]
mu_3<-telemetry[(2*n+3):(3*n+3)]	
r=coeff(z_12,z_13,z_23,h_1,h_2,h_3)
val<-0
val<-dmvnorm( c(mu_1[-1],mu_2[-1],mu_3[-1]),c(rep(mu_1[1],n),rep(mu_2[1],n),rep(mu_3[1],n)), cov_OUH(T,n,r[1],r[2],r[3],s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3),T)
return(val)
}



log_likelihood_3d_prac<-function(T,n,telemetry,z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3,index){
k<-length(index)
mu_1<-telemetry[1:(k+1)]
mu_2<-telemetry[(k+2):(2*k+2)]
mu_3<-telemetry[(2*k+3):(3*k+3)]	
r=coeff(z_12,z_13,z_23,h_1,h_2,h_3)
val<-0
K<-cov_OUH(T,n,r[1],r[2],r[3],s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)
K_2<-K[c(index,n+index,2*n+index),c(index,n+index,2*n+index)]
val<-dmvnorm( c(mu_1[-1],mu_2[-1],mu_3[-1]),c(rep(mu_1[1],k),rep(mu_2[1],k),rep(mu_3[1],k)), K_2,T)
return(val)
}




log_likelihood_3d_prac_r<-function(T,n,telemetry,r_12,r_13,r_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3,index){
k<-length(index)
mu_1<-telemetry[1:(k+1)]
mu_2<-telemetry[(k+2):(2*k+2)]
mu_3<-telemetry[(2*k+3):(3*k+3)]	
val<-0
K<-cov_OUH(T,n,r_12,r_13,r_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)
K_2<-K[c(index,n+index,2*n+index),c(index,n+index,2*n+index)]
val<-dmvnorm( c(mu_1[-1],mu_2[-1],mu_3[-1]),c(rep(mu_1[1],k),rep(mu_2[1],k),rep(mu_3[1],k)), K_2,T)
return(val)
}




simulation_3d<-function(mu_ini_1,mu_ini_2,mu_ini_3,z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3){
r=coeff(z_12,z_13,z_23,h_1,h_2,h_3)
telemetry<-rep(0,3*(n+1))
telemetry[c(2:(n+1),(n+3):(2*n+2),(2*n+4):(3*n+3))]<-mvrnorm( 1,c(rep(mu_ini_1,n),rep(mu_ini_2,n),rep(mu_ini_3,n)), cov_OUH(T,n,r[1],r[2],r[3],s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3),T)
telemetry[c(1,n+2,2*n+3)]<-c(mu_ini_1,mu_ini_2,mu_ini_3)
return(telemetry)
}





install.packages("pracma")
install.packages("MASS")
install.packages("mvtnorm")
install.packages("optimr")
install.packages("bbmle")
install.packages("mle.tools")
install.packages("plot3D")
install.packages("rworldmap")
install.packages("geosphere")
install.packages("numDeriv")


library(pracma)
library(MASS)
library(mvtnorm)
library(optimr)
library(bbmle)
library(mle.tools)
library(plot3D)
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

#-------Simulacion inferencia 3D--------



#-----------------------Simulaciones 3D------------------------------------------------#


#-----Parametros movimiento 1
s_1<-2
b_1<-7
h_1<-0.75
mu_i_1<-0
v_i_1<-0
#------Parametros movimiento 2
s_2<-3
b_2<-10
h_2<-0.5
mu_i_2<-0
v_i_2<-0


#------Parametros movimiento 3
s_3<-2.5
b_3<-5
h_3<-0.25
mu_i_3<-0
v_i_3<-0


z_12<-0.2
z_13<--0.5
z_23<-0.7

#max



#-----Varibles globales
T<-200
n<-400
t<-seq(0,T,length=n)
D<-T/n






#----------Imagen funcion de correlacion-----

r=coeff(z_12,z_13,z_23,h_1,h_2,h_3)


M<-cov_OUH(T,n,r[1],r[2],r[3],s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)


k<-50


par(mfrow=c(2,3))


col<-c(1:(n/k))


ymin=min(M[1:n,1:n])
ymax=max(M[1:n,1:n])


plot(t,M[1:n,k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['1,1'](s,t)),cex.lab=1.1,col=col[1])
for(i in c(2:(n/k))){
lines(t,M[1:n,i*k],lwd=2,col=col[i])
}


#a<-(k/2)*c(1:(n/k))
 
#legend("topleft", legend=c(paste("s=",a)),
   #    col=c(1:(n/k)), lty=1, cex=1.1)




ymin=min(M[(n+1):(2*n),1:n])
ymax=max(M[(n+1):(2*n),1:n])



plot(t,M[(n+1):(2*n),k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['1,2'](s,t)),cex.lab=1.1,,col=col[1])
for(i in c(2:(n/k))){
lines(t,M[(n+1):(2*n),i*k],lwd=2,col=col[i])
}






ymin=min(M[(2*n+1):(3*n),1:n])
ymax=max(M[(2*n+1):(3*n),1:n])



plot(t,M[(2*n+1):(3*n),k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['1,3'](s,t)),cex.lab=1.1,col=col[1])


for(i in c(2:(n/k))){
lines(t,M[(2*n+1):(3*n),i*k],lwd=2,col=col[i])
}




ymin=min(M[(n+1):(2*n),(n+1):(2*n)])
ymax=max(M[(n+1):(2*n),(n+1):(2*n)])



plot(t,M[(n+1):(2*n),n+k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['2,2'](s,t)),cex.lab=1.1,col=col[1])


for(i in c(2:(n/k))){
lines(t,M[(n+1):(2*n),n+i*k],lwd=2,col=col[i])
}





ymin=min(M[(n+1):(2*n),(2*n+1):(3*n)])
ymax=max(M[(n+1):(2*n),(2*n+1):(3*n)])

plot(t,M[(2*n+1):(3*n),n+k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['2,3'](s,t)),cex.lab=1.1,col=col[1])

for(i in c(2:(n/k))){
lines(t,M[(2*n+1):(3*n),n+i*k],lwd=2,col=col[i])
}







ymin=min(M[(2*n+1):(3*n),(2*n+1):(3*n)])
ymax=max(M[(2*n+1):(3*n),(2*n+1):(3*n)])


plot(t,M[(2*n+1):(3*n),2*n+k],type="l",lwd=2,xlim=c(0,T+0.5),ylim=c(ymin,ymax),ylab=expression(R['3,3'](s,t)),cex.lab=1.1,col=col[1])


for(i in c(2:(n/k))){
lines(t,M[(2*n+1):(3*n),2*n+i*k],lwd=2,col=col[i])
}


#-------Fin imagen funcion de correlacion----




#-----Simulaciones

telemetry<-simulation_3d(mu_i_1,mu_i_2,mu_i_3,z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)


mu_1<-telemetry[1:(n+1)]
mu_2<-telemetry[(n+2):(2*n+2)]
mu_3<-telemetry[(2*n+3):(3*n+3)]


log_likelihood_3d(c(mu_1[-1],mu_2[-1],mu_3[-1]),z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)



#---Imagen e imagen con angulo-----

plot_ly(x=mu_1,y=mu_2,z=mu_3,type = 'scatter3d',mode = 'lines+markers')
x<-rep(0,n)
y<-rep(0,n)
z<-rep(0,n)
w_1<-rep(0,n)
w_2<-rep(0,n)
w_3<-rep(0,n)
aux_x<-0
aux_y<-0
aux_z<-0
for(i in 1:n){
	aux_x<-mu_1[i+1]-mu_1[i]
	aux_y<-mu_2[i+1]-mu_2[i]
	aux_z<-mu_3[i+1]-mu_3[i]
	x[i]<-aux_x/norm(c(aux_x,aux_y))
	y[i]<-aux_y/norm(c(aux_x,aux_y))
	w_3[i]<-aux_x/norm(c(aux_x,aux_z))
	z[i]<-aux_z/norm(c(aux_x,aux_z))
	w_1[i]<-aux_y/norm(c(aux_y,aux_z))
	w_2[i]<-aux_z/norm(c(aux_y,aux_z))
}

par(mfrow=c(1,3))

plot(x,y,main="Angulo (x,y)")
lines(c(0,0),c(-1,1),col="green")
lines(c(-1,1),c(0,0),col="green")

plot(w_3,z,main="Angulo (x,z)")
lines(c(0,0),c(-1,1),col="green")
lines(c(-1,1),c(0,0),col="green")


plot(w_1,w_2,main="Angulo (y,z)")
lines(c(0,0),c(-1,1),col="green")
lines(c(-1,1),c(0,0),col="green")



#-------------Inferencia coordenada a coordenada------------#


#--------mu_1

#------------Change to lon data---------------
log_vero<-function(x){
-log_verosimilitude(T,x[1],x[2],mu_1)
}


#----------------Imprimirmos la perfil de h y b para tener una idea de donde puede el EMV de h


R_2<-1
b_a<-b_1/5
b_b<-b_1*5



block<-75
reg_b<-seq(b_a*R_2,b_b*R_2,length=block)
reg_h<-seq(0.1,0.9,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<-log_vero(c(reg_b[i],reg_h[j]))
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")



max_mu_1<-rep(0,3)
mu<-mu_1
maximos_mu_1<-opm(c(10,0.7),fn=log_vero, lower=c(3,0.5), upper=c(15,0.95),method="L-BFGS-B")
 max_mu_1[2]<-maximos_mu_1$p1
 max_mu_1[3]<-maximos_mu_1$p2
 K_2<-cov_OUH_zeta_2(T,n,max_mu_1[2],max_mu_1[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_1[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_mu_1[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)

p<-exp(-log_vero(c(max_mu_1[2],max_mu_1[3]))/3)
res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")
mypoints <- trans3d(max_mu_1[2],max_mu_1[3],p,pmat=res)
points(mypoints, pch=19, col="red",cex=1)



#---------mu_2

#------------Change to lon data---------------
log_vero<-function(x){
-log_verosimilitude(T,x[1],x[2],mu_2)
}


#----------------Imprimirmos la perfil de h y b para tener una idea de donde puede el EMV de h


R_2<-1
b_a<-b_2/5
b_b<-b_2*5



block<-75
reg_b<-seq(b_a*R_2,b_b*R_2,length=block)
reg_h<-seq(0.1,0.9,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<-log_vero(c(reg_b[i],reg_h[j]))
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")



max_mu_2<-rep(0,3)
mu<-mu_2
maximos_mu_2<-opm(c(10,0.7),fn=log_vero, lower=c(3,0.5), upper=c(15,0.95),method="L-BFGS-B")
 max_mu_2[2]<-maximos_mu_2$p1
 max_mu_2[3]<-maximos_mu_2$p2
 K_2<-cov_OUH_zeta_2(T,n,max_mu_2[2],max_mu_2[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_2[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_mu_2[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)

p<-exp(-log_vero(c(max_mu_2[2],max_mu_2[3]))/3)
res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")
mypoints <- trans3d(max_mu_2[2],max_mu_2[3],p,pmat=res)
points(mypoints, pch=19, col="red",cex=1)

#-------mu_3------


#------------Change to lon data---------------
log_vero<-function(x){
-log_verosimilitude(T,x[1],x[2],mu_3)
}


#----------------Imprimirmos la perfil de h y b para tener una idea de donde puede el EMV de h


R_2<-1
b_a<-b_3/5
b_b<-b_3*5



block<-75
reg_b<-seq(b_a*R_2,b_b*R_2,length=block)
reg_h<-seq(0.1,0.9,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<-log_vero(c(reg_b[i],reg_h[j]))
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")



max_mu_3<-rep(0,3)
mu<-mu_3
maximos_mu_3<-opm(c(10,0.3),fn=log_vero, lower=c(1.5,0.05), upper=c(15,0.5),method="L-BFGS-B")
 max_mu_3[2]<-maximos_mu_3$p1
 max_mu_3[3]<-maximos_mu_3$p2
 K_2<-cov_OUH_zeta_2(T,n,max_mu_3[2],max_mu_3[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_3[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_mu_3[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)

p<-exp(-log_vero(c(max_mu_3[2],max_mu_3[3]))/3)
res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")
mypoints <- trans3d(max_mu_3[2],max_mu_3[3],p,pmat=res)
points(mypoints, pch=19, col="red",cex=1)



#--------valores iniciales-----#

max_mu_1
max_mu_2
max_mu_3

write.csv(telemetry,"telemetry.csv")
write.csv(c(max_mu_1,max_mu_2,max_mu_3),"initial.csv")



#-----Lectura datos----

telemetry<-data.frame(read.csv("telemetry.csv"))
mu_1<-telemetry[1:(n+1),2]
mu_2<-telemetry[(n+2):(2*n+2),2]
mu_3<-telemetry[(2*n+3):(3*n+3),2]

initial<-data.frame(read.csv("initial.csv"))
max_mu_1<-initial[1:3,2]
max_mu_2<-initial[4:6,2]
max_mu_3<-initial[7:9,2]


z_12<-runif(1,-1,1)
z_13<-runif(1,-1,1)
z_23<-runif(1,-1,1)

log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),z_12,z_13,z_23,max_mu_1[1],max_mu_2[1],max_mu_3[1],max_mu_1[2],max_mu_2[2],max_mu_3[2],max_mu_1[3],max_mu_2[3],max_mu_3[3])





#------------Change to lon data---------------
log_likelihood<-function(x){
-log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12])
}


#----------------Imprimirmos la perfil de h y b para tener una idea de donde puede el EMV de h


maximos<-rep(0,12)
maximos<-opm(c(z_12,z_13,z_23,max_mu_1[1],max_mu_2[1],max_mu_3[1],max_mu_1[2],max_mu_2[2],max_mu_3[2],max_mu_1[3],max_mu_2[3],max_mu_3[3]),fn=log_likelihood, lower=c(-0.98,-0.98,-0.98,0.2,0.5,0.5,0.5,0.5,0.5,0.51,0.51,0.01), upper=c(0.98,0.98,0.98,5,7,7,12,12,10,0.95,0.95,0.49),method="L-BFGS-B")

max<-rep(0,12)

max[1]<-maximos$p1
max[2]<-maximos$p2
max[3]<-maximos$p3
max[4]<-maximos$p4
max[5]<-maximos$p5
max[6]<-maximos$p6
max[7]<-maximos$p7
max[8]<-maximos$p8
max[9]<-maximos$p9
max[10]<-maximos$p10
max[11]<-maximos$p11
max[12]<-maximos$p12

write.csv(max,"max_val.csv")
max<-data.frame(read.csv("max_val.csv"))$x



maximo<--log_likelihood(max)


#--------------------Profile likelihood---------

#------------Change to lon data---------------
log_likelihood_1<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),x,max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12])
}

log_likelihood_2<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],x,max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_3<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],x,max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_4<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],x,max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_5<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],x,max[6],max[7],max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_6<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],x,max[7],max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_7<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],x,max[8],max[9],max[10],max[11],max[12])
}
log_likelihood_8<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],max[7],x,max[9],max[10],max[11],max[12])
}
log_likelihood_9<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],x,max[10],max[11],max[12])
}
log_likelihood_10<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],x,max[11],max[12])
}
log_likelihood_11<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],x,max[12])
}
log_likelihood_12<-function(x){
log_likelihood_3d(T,n,c(mu_1,mu_2,mu_3),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],x)
}


par(mfrow=c(4,3))

max


val_right<-rep(0,12)
val_left<-rep(0,12)

#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.22,0.5,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_1(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_1,max[1]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[1],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[1],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[1],sqrt(norm_1))
val1<-dnorm(valc1_1,max[1],sqrt(norm_1))

val_left[1]<-valc1_1
val_right[1]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[12])),xlab=expression(z[12]),ylab=expression(paste("Profile Likelihood ",z[12])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[1],max[1]),c(0,12),col="red",lwd=2)
lines(c(z_12,z_12),c(0,12),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,12),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,12),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------

#-----------rho_13 #------(-0.1,0.1) 
m<-50
reg_aux<-seq(-0.56,-0.3,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_2(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_2,max[2]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[2],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[2],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[2],sqrt(norm_1))
val1<-dnorm(valc1_1,max[2],sqrt(norm_1))

val_left[2]<-valc1_1
val_right[2]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[13])),xlab=expression(z[13]),ylab=expression(paste("Profile Likelihood ",z[13])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[2],max[2]),c(0,13),col="red",lwd=2)
lines(c(z_13,z_13),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,13),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,13),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)



#----------


#-----------rho_23 #----(-0.05,0.15)
m<-50
reg_aux<-seq(0.45,0.75,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_3(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_3,max[3]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[3],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[3],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[3],sqrt(norm_1))
val1<-dnorm(valc1_1,max[3],sqrt(norm_1))

val_left[3]<-valc1_1
val_right[3]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[23])),xlab=expression(z[23]),ylab=expression(paste("Profile Likelihood ",z[23])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[3],max[3]),c(0,15),col="red",lwd=2)
lines(c(z_23,z_23),c(0,15),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,15),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,15),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------



#-----------s_1  #----(0.05,0.065)
m<-50
reg_aux<-seq(1.8,2.5,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_4(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_4,max[4]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[4],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[4],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[4],sqrt(norm_1))
val1<-dnorm(valc1_1,max[4],sqrt(norm_1))

val_left[4]<-valc1_1
val_right[4]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[1])),xlab=expression(sigma[1]),ylab=expression(paste("Profile Likelihood ",sigma[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[4],max[4]),c(0,20),col="red",lwd=2)
lines(c(s_1,s_1),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,20),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,20),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------




#-----------s_2  #------(0.029,0.038)
m<-50
reg_aux<-seq(2.25,3.25,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_5(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_5,max[5]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[5],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[5],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[5],sqrt(norm_1))
val1<-dnorm(valc1_1,max[5],sqrt(norm_1))

val_left[5]<-valc1_1
val_right[5]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[2])),xlab=expression(sigma[2]),ylab=expression(paste("Profile Likelihood ",sigma[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[5],max[5]),c(0,20),col="red",lwd=2)
lines(c(s_2,s_2),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,20),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,20),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------s_3   #------(0.028,0.037)
m<-50
reg_aux<-seq(2.1,2.9,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_6(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_6,max[6]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[6],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[6],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[6],sqrt(norm_1))
val1<-dnorm(valc1_1,max[6],sqrt(norm_1))

val_left[6]<-valc1_1
val_right[6]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[3])),xlab=expression(sigma[3]),ylab=expression(paste("Profile Likelihood ",sigma[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[6],max[6]),c(0,250),col="red",lwd=2)
lines(c(s_3,s_3),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------


#-----------b_1  #----------(9,12)
m<-50
reg_aux<-seq(5,11,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_7(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_7,max[7]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[7],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[7],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[7],sqrt(norm_1))
val1<-dnorm(valc1_1,max[7],sqrt(norm_1))

val_left[7]<-valc1_1
val_right[7]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[1])),xlab=expression(beta[1]),ylab=expression(paste("Profile Likelihood ",beta[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[7],max[7]),c(0,250),col="red",lwd=2)
lines(c(b_1,b_1),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------



#-----------b_2  #-----(8.5,11)    
m<-50
reg_aux<-seq(0.75,6.75,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_8(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_8,max[8]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[8],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[8],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[8],sqrt(norm_1))
val1<-dnorm(valc1_1,max[8],sqrt(norm_1))

val_left[8]<-valc1_1
val_right[8]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[2])),xlab=expression(beta[2]),ylab=expression(paste("Profile Likelihood ",beta[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[8],max[8]),c(0,250),col="red",lwd=2)
lines(c(b_2,b_2),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------



#-----------b_3  #-----(0.6,1.2)
m<-50
reg_aux<-seq(0.52,6.52,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_9(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_9,max[9]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[9],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[9],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[9],sqrt(norm_1))
val1<-dnorm(valc1_1,max[9],sqrt(norm_1))

val_left[9]<-valc1_1
val_right[9]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[3])),xlab=expression(beta[3]),ylab=expression(paste("Profile Likelihood ",beta[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[9],max[9]),c(0,250),col="red",lwd=2)
lines(c(b_3,b_3),c(0,13),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_1   #--(0.89,0.94)
m<-50
reg_aux<-seq(0.77,0.83,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_10(reg_aux[i])-maximo
}


#hessian(log_likelihood_10,max[10])
norm_1<-1/(-hessian(log_likelihood_10,max[10]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[10],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[10],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[10],sqrt(norm_1))
val1<-dnorm(valc1_1,max[10],sqrt(norm_1))

val_left[10]<-valc1_1
val_right[10]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[1])),xlab=expression(H[1]),ylab=expression(paste("Profile Likelihood ",H[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[10],max[10]),c(0,250),col="red",lwd=2)
lines(c(h_1,h_1),c(0,250),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_2   #---(0.8,0.88)
m<-50
reg_aux<-seq(0.52,0.62,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_11(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_11,max[11]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[11],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[11],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[11],sqrt(norm_1))
val1<-dnorm(valc1_1,max[11],sqrt(norm_1))

val_left[11]<-valc1_1
val_right[11]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[2])),xlab=expression(H[2]),ylab=expression(paste("Profile Likelihood ",H[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[11],max[11]),c(0,250),col="red",lwd=2)
lines(c(h_2,h_2),c(0,250),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_3  #------(0.22,0.47)
m<-50
reg_aux<-seq(0.08,0.25,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_12(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_12,max[12]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[12],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[12],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[12],sqrt(norm_1))
val1<-dnorm(valc1_1,max[12],sqrt(norm_1))

val_left[12]<-valc1_1
val_right[12]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[3])),xlab=expression(H[3]),ylab=expression(paste("Profile Likelihood ",H[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[12],max[12]),c(0,250),col="red",lwd=2)
lines(c(h_3,h_3),c(0,250),col="yellow",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#---------









par(mfrow=c(4,3))

f<-qchisq(0.99,1)

root_1_l<-rep(0,12)
root_1_r<-rep(0,12)


#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.25,0.47,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_1(reg_aux[i])-maximo
}


like_ratio_1<-function(x){
	-2*(log_likelihood_1(x)-maximo)-f
}


root_1_l[1]<-uniroot(like_ratio_1,lower=0.22,upper=max[1],tol=1e-9)$root
root_1_r[1]<-uniroot(like_ratio_1,lower=max[1],upper=0.5,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[12])),xlab=expression(z[12]),ylab=expression(paste("Likelihood Ratio Statistic ",z[12])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.22,0.5),c(f,f),col="blue",lwd=2)
lines(c(max[1],max[1]),c(0,15),col="red",lwd=2)
lines(c(z_12,z_12),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[1],root_1_l[1]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[1],root_1_r[1]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------rho_13 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.52,-0.31,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_2(reg_aux[i])-maximo
}


like_ratio_2<-function(x){
	-2*(log_likelihood_2(x)-maximo)-f
}


root_1_l[2]<-uniroot(like_ratio_2,lower=-0.56,upper=max[2],tol=1e-9)$root
root_1_r[2]<-uniroot(like_ratio_2,lower=max[2],upper=-0.3,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[13])),xlab=expression(z[13]),ylab=expression(paste("Likelihood Ratio Statistic ",z[13])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.56,-0.3),c(f,f),col="blue",lwd=2)
lines(c(max[2],max[2]),c(0,15),col="red",lwd=2)
lines(c(z_13,z_13),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[2],root_1_l[2]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[2],root_1_r[2]),c(0,15),col="green",lwd=2)

#-------------------------------------------------



#-----------rho_23 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.5,0.69,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_3(reg_aux[i])-maximo
}


like_ratio_3<-function(x){
	-2*(log_likelihood_3(x)-maximo)-f
}


root_1_l[3]<-uniroot(like_ratio_3,lower=0.45,upper=max[3],tol=1e-9)$root
root_1_r[3]<-uniroot(like_ratio_3,lower=max[3],upper=0.71,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[23])),xlab=expression(z[23]),ylab=expression(paste("Likelihood Ratio Statistic ",z[23])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.5,0.7),c(f,f),col="blue",lwd=2)
lines(c(max[3],max[3]),c(0,15),col="red",lwd=2)
lines(c(z_23,z_23),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[3],root_1_l[3]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[3],root_1_r[3]),c(0,15),col="green",lwd=2)

#-------------------------------------------------



#-----------s_1 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(1.95,2.4,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_4(reg_aux[i])-maximo
}


like_ratio_4<-function(x){
	-2*(log_likelihood_4(x)-maximo)-f
}


root_1_l[4]<-uniroot(like_ratio_4,lower=1.8,upper=max[4],tol=1e-9)$root
root_1_r[4]<-uniroot(like_ratio_4,lower=max[4],upper=2.4,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[1])),xlab=expression(sigma[1]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(1.8,2.4),c(f,f),col="blue",lwd=2)
lines(c(max[4],max[4]),c(0,15),col="red",lwd=2)
lines(c(s_1,s_1),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[4],root_1_l[4]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[4],root_1_r[4]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------s_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(2.5,3.1,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_5(reg_aux[i])-maximo
}


like_ratio_5<-function(x){
	-2*(log_likelihood_5(x)-maximo)-f
}


root_1_l[5]<-uniroot(like_ratio_5,lower=2.2,upper=max[5],tol=1e-9)$root
root_1_r[5]<-uniroot(like_ratio_5,lower=max[5],upper=3.2,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[2])),xlab=expression(sigma[2]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(2.2,3.2),c(f,f),col="blue",lwd=2)
lines(c(max[5],max[5]),c(0,15),col="red",lwd=2)
lines(c(s_2,s_2),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[5],root_1_l[5]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[5],root_1_r[5]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------s_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(2.25,2.8,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_6(reg_aux[i])-maximo
}


like_ratio_6<-function(x){
	-2*(log_likelihood_6(x)-maximo)-f
}


root_1_l[6]<-uniroot(like_ratio_6,lower=2.1,upper=max[6],tol=1e-9)$root
root_1_r[6]<-uniroot(like_ratio_6,lower=max[6],upper=2.9,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[3])),xlab=expression(sigma[3]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(2.1,2.9),c(f,f),col="blue",lwd=2)
lines(c(max[6],max[6]),c(0,15),col="red",lwd=2)
lines(c(s_3,s_3),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[6],root_1_l[6]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[6],root_1_r[6]),c(0,15),col="green",lwd=2)

#-------------------------------------------------





#-----------b_1 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(5.5,10.5,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_7(reg_aux[i])-maximo
}


like_ratio_7<-function(x){
	-2*(log_likelihood_7(x)-maximo)-f
}


root_1_l[7]<-uniroot(like_ratio_7,lower=3,upper=max[7],tol=1e-9)$root
root_1_r[7]<-uniroot(like_ratio_7,lower=max[7],upper=13,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[1])),xlab=expression(beta[1]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(3,13),c(f,f),col="blue",lwd=2)
lines(c(max[7],max[7]),c(0,15),col="red",lwd=2)
lines(c(b_1,b_1),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[7],root_1_l[7]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[7],root_1_r[7]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------b_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(1.9,5.5,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_8(reg_aux[i])-maximo
}


like_ratio_8<-function(x){
	-2*(log_likelihood_8(x)-maximo)-f
}


root_1_l[8]<-uniroot(like_ratio_8,lower=1.5,upper=max[8],tol=1e-9)$root
root_1_r[8]<-uniroot(like_ratio_8,lower=max[8],upper=6.2,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[2])),xlab=expression(beta[2]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(1.5,6.2),c(f,f),col="blue",lwd=2)
lines(c(max[8],max[8]),c(0,15),col="red",lwd=2)
lines(c(b_2,b_2),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[8],root_1_l[8]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[8],root_1_r[8]),c(0,15),col="green",lwd=2)

#-------------------------------------------------


#-----------b_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(1.65,5.35,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_9(reg_aux[i])-maximo
}


like_ratio_9<-function(x){
	-2*(log_likelihood_9(x)-maximo)-f
}


root_1_l[9]<-uniroot(like_ratio_9,lower=1.65,upper=max[9],tol=1e-9)$root
root_1_r[9]<-uniroot(like_ratio_9,lower=max[9],upper=5.35,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[3])),xlab=expression(beta[3]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(1.65,5.35),c(f,f),col="blue",lwd=2)
lines(c(max[9],max[9]),c(0,15),col="red",lwd=2)
lines(c(b_3,b_3),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[9],root_1_l[9]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[9],root_1_r[9]),c(0,15),col="green",lwd=2)

#-------------------------------------------------






#-----------H_1 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.78,0.82,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_10(reg_aux[i])-maximo
}


like_ratio_10<-function(x){
	-2*(log_likelihood_10(x)-maximo)-f
}


root_1_l[10]<-uniroot(like_ratio_10,lower=0.76,upper=max[10],tol=1e-9)$root
root_1_r[10]<-uniroot(like_ratio_10,lower=max[10],upper=0.83,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[1])),xlab=expression(H[1]),ylab=expression(paste("Likelihood Ratio Statistic ",H[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.76,0.83),c(f,f),col="blue",lwd=2)
lines(c(max[10],max[10]),c(0,15),col="red",lwd=2)
lines(c(h_1,h_1),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[10],root_1_l[10]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[10],root_1_r[10]),c(0,15),col="green",lwd=2)

#-------------------------------------------------





#-----------H_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.535,0.602,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_11(reg_aux[i])-maximo
}


like_ratio_11<-function(x){
	-2*(log_likelihood_11(x)-maximo)-f
}


root_1_l[11]<-uniroot(like_ratio_11,lower=0.52,upper=max[11],tol=1e-9)$root
root_1_r[11]<-uniroot(like_ratio_11,lower=max[11],upper=0.62,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[2])),xlab=expression(H[2]),ylab=expression(paste("Likelihood Ratio Statistic ",H[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.52,0.62),c(f,f),col="blue",lwd=2)
lines(c(max[11],max[11]),c(0,20),col="red",lwd=2)
lines(c(h_2,h_2),c(0,20),col="yellow",lwd=2)
lines(c(root_1_l[11],root_1_l[11]),c(0,20),col="green",lwd=2)
lines(c(root_1_r[11],root_1_r[11]),c(0,20),col="green",lwd=2)

#-------------------------------------------------


#-----------H_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.10,0.23,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_12(reg_aux[i])-maximo
}


like_ratio_12<-function(x){
	-2*(log_likelihood_12(x)-maximo)-f
}


root_1_l[12]<-uniroot(like_ratio_12,lower=0.10,upper=max[12],tol=1e-9)$root
root_1_r[12]<-uniroot(like_ratio_12,lower=max[12],upper=0.23,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[3])),xlab=expression(H[3]),ylab=expression(paste("Likelihood Ratio Statistic ",H[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.10,0.23),c(f,f),col="blue",lwd=2)
lines(c(max[12],max[12]),c(0,15),col="red",lwd=2)
lines(c(h_3,h_3),c(0,15),col="yellow",lwd=2)
lines(c(root_1_l[12],root_1_l[12]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[12],root_1_r[12]),c(0,15),col="green",lwd=2)

#-------------------------------------------------------------

#----Intervalos 99% confianza---


val_left
val_right


root_1_l
root_1_r



#------Etiquetas--------
x<-c(1:10)
plot(x,x,col="white")
legend("topleft", legend=c("Profile Likelihood","Normal Approximation","True Value","MLE", "99% Confidence Interval Normal Approximation"),
       col=c("black","blue","yellow","red","green"), lty=1,cex=1.3)



plot(x,x,col="white")
legend("topleft", legend=c("Likelihood Ratio Statistic",expression(chi['1,0.99']^2),"True Value","MLE", "99% Confidence Interval Chi-Squared"),
       col=c("black","blue","yellow","red","green"), lty=1,cex=1.3)


plot(x,x,col="white")
legend("topleft", legend=c("Likelihood Ratio Statistic",expression(chi['1,0.99']^2),"MLE", "99% Confidence Interval Chi-Squared"),
       col=c("black","blue","red","green"), lty=1,cex=1.3)

plot(x,x,col="white")
legend("topleft", legend=c("Profile Likelihood","Normal Approximation","MLE", "99% Confidence Interval Normal Approximation"),
       col=c("black","blue","red","green"), lty=1,cex=1.3)

#------------------


#--------Imagenes de la funcion de covarianza---------











#--------------Fin inferencia simulacion-------


#---------------Caso practico bat

telemetry_data<-data.frame(read.csv("Bat53D.csv"))
k<-length(telemetry_data$lon)
mu_1<-telemetry_data$lon
mu_2<-telemetry_data$lat
mu_3<-telemetry_data$alt


telemetry_2D<-degree_km(mu_1,mu_2)/1000

mu_1_r<-telemetry_2D[1,]
mu_2_r<-telemetry_2D[2,]
mu_3_r<-mu_3/1000






simulation_vel<-simulation_vel(T,max[4],max[7],max[10],mu_1,index,2)



plot_ly(x=mu_1_r,y=mu_2_r,z=mu_3_r,type = 'scatter3d',mode = 'lines+markers')
time<-telemetry_data$time
index<-time[-1]*2


k<-length(index)

min(mu_1)
min(mu_2)
max(mu_1)
max(mu_2)



T<-index[k]/2
D<-0.5
n<-T/D
t<-seq(0,T,length=(n+1))


#-------------Inferencia coordenada a coordenada------------#

#--------mu_1


block<-30
reg_b<-seq(1,20,length=block)
reg_h<-seq(0.1,0.99,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<--log_verosimilitude_practical(T,reg_b[i],reg_h[j],mu_1,index)
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/30),theta = 180, phi = 0,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")




log_vero<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],mu_1,index)
}


max_mu_1<-rep(0,3)
maximos_mu_1<-opm(c(10,0.7),fn=log_vero, lower=c(0.01,0.5), upper=c(30,0.99),method="L-BFGS-B")
 #max_mu_1[2]<-maximos_mu_1$p1
 max_mu_1[3]<-maximos_mu_1$p2


#----------Imprimirmos la perfil de b-----------------#

log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,max_mu_1[3],mu_1,index)
}


reg_b<-seq(0.1,15,length=100)
reg_b_ver<-rep(0,length(reg_b))
for(i in 1:length(reg_b)){
reg_b_ver[i]<--log_vero_pract(reg_b[i])
}
plot(reg_b,reg_b_ver,xlab="b",ylab="L_p(b)",type="l",col="blue",main="Likelihood beta")


max_mu_1[2]<-runif(1,5,10)


 K_2<-cov_OUH_zeta_2(T,n,max_mu_1[2],max_mu_1[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_1[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_mu_1[1]<-sqrt((t(as.vector(mu_1[2:(k+1)]-mu_1[1]))%*%solve(M)%*%as.vector(mu_1[2:(k+1)]-mu_1[1]) )[1,1]/k)




#-------------Inferencia coordenada a coordenada------------#

#--------mu_2


block<-30
reg_b<-seq(5,25,length=block)
reg_h<-seq(0.7,0.9,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<--log_verosimilitude_practical(T,reg_b[i],reg_h[j],mu_2,index)
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/30),theta = 45, phi = 0,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")




log_vero<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],mu_2,index)
}




max_mu_2<-rep(0,3)
maximos_mu_2<-opm(c(10,0.7),fn=log_vero, lower=c(2,0.5), upper=c(20,0.95),method="L-BFGS-B")
 #max_mu_1[2]<-maximos_mu_1$p1
 max_mu_2[3]<-maximos_mu_2$p2



#----------Imprimirmos la perfil de b-----------------#

log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,max_mu_2[3],mu_2,index)
}


reg_b<-seq(0.1,15,length=100)
reg_b_ver<-rep(0,length(reg_b))
for(i in 1:length(reg_b)){
reg_b_ver[i]<--log_vero_pract(reg_b[i])
}
plot(reg_b,reg_b_ver,xlab="b",ylab="L_p(b)",type="l",col="blue",main="Likelihood beta")


max_mu_2[2]<-runif(1,5,10)


 K_2<-cov_OUH_zeta_2(T,n,max_mu_2[2],max_mu_2[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_2[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_mu_2[1]<-sqrt((t(as.vector(mu_2[2:(k+1)]-mu_2[1]))%*%solve(M)%*%as.vector(mu_2[2:(k+1)]-mu_2[1]) )[1,1]/k)


#-------------Inferencia coordenada a coordenada------------#

#--------mu_3



block<-30
reg_b<-seq(1,8,length=block)
reg_h<-seq(0.5,0.97,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<--log_verosimilitude_practical(T,reg_b[i],reg_h[j],mu_3_r,index)
 }
}

res<-persp(reg_b,reg_h,exp(-profile_bh/30),theta = 270, phi = 0,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")




log_vero<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],mu_3_r,index)
}




max_mu_3<-rep(0,3)
maximos_mu_3<-opm(c(10,0.5),fn=log_vero, lower=c(0.01,0.1), upper=c(20,0.9),method="L-BFGS-B")
 max_mu_3[2]<-maximos_mu_3$p1
 max_mu_3[3]<-maximos_mu_3$p2
max_mu_3

#----------Imprimirmos la perfil de b-----------------#

 
 
 
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,max_mu_3[3],mu_3,index)
}


reg_b<-seq(0.1,15,length=100)
reg_b_ver<-rep(0,length(reg_b))
for(i in 1:length(reg_b)){
reg_b_ver[i]<--log_vero_pract(reg_b[i])
}
plot(reg_b,reg_b_ver,xlab="b",ylab="L_p(b)",type="l",col="blue",main="Likelihood beta")


 
 K_2<-cov_OUH_zeta_2(T,n,max_mu_3[2],max_mu_3[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_mu_3[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_mu_3[1]<-sqrt((t(as.vector(mu_3_r[2:(k+1)]-mu_3_r[1]))%*%solve(M)%*%as.vector(mu_3_r[2:(k+1)]-mu_3_r[1]) )[1,1]/k)



#-----Correlation

plot(mu_1,mu_2)
z_12<-runif(1,0,1)

plot(mu_1,mu_3)
z_13<-runif(1,0,1)

plot(mu_2,mu_3)
z_23<-runif(1,0,1)



#------maximos----
max_mu_1
max_mu_2
max_mu_3





log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),0,0,0,max_mu_1[1],max_mu_2[1],max_mu_3[1],max_mu_1[2],max_mu_2[2],max_mu_3[2],max_mu_1[3],max_mu_2[3],max_mu_3[3],index)


log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),0.3,0.3,0.3,0.2,0.2,0.2,11,10,1.5,0.94,0.92,0.4,index)

log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),-0.5,-0.5,-0.5,0.005,0.005,0.005,1,1,0.1,0.6,0.6,0.1,index)

#------------Change to lon data---------------
log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}




maximos<-rep(0,12)
maximos<-opm(c(0,0,0,max_mu_1[1],max_mu_2[1],max_mu_3[1],max_mu_1[2],max_mu_2[2],max_mu_3[2],max_mu_1[3],max_mu_2[3],max_mu_3[3]),fn=log_likelihood_prac, lower=c(-0.5,-0.5,-0.5,0.005,0.005,0.005,1,1,0.1,0.6,0.6,0.1), upper=c(0.3,0.3,0.3,0.2,0.2,0.2,11,10,1.5,0.96,0.92,0.4),method="L-BFGS-B")




log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),-0.1646701,-0.01833508,0.03301999,0.05364833,0.02982051,0.03182533,x,8.921156,0.8937402,0.9210071,0.8504374,0.3092366,index),-10^7)
}

aux_max<-opm(10,fn=log_likelihood_prac, lower=8, upper=11,method="L-BFGS-B")




log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),-0.1646701,0,0,0.05364833,0.02982051,0.03182533,10.14239,8.921156,0.8937402,0.9210071,0.8504374,0.3092366,index)

max<-rep(0,12)

max[1]<--0.1646701
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05364833
max[5]<-0.02982051
max[6]<-0.03182533
max[7]<-10.14239
max[8]<-8.921156
max[9]<-0.8937402
max[10]<-0.9210071
max[11]<-0.8504374
max[12]<-0.3092366

write.csv(max,"max_val_Bat53D.csv")

maximo<--log_likelihood_prac(max)


log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],x[1],max[5],max[6],max[7],x[2],max[9],x[3],max[11],max[12],index),-10^7)
}

aux_max_2<-opm(c(max[4],max[8],max[10]),fn=log_likelihood_prac, lower=c(0.053,8.7,0.915), upper=c(0.057,9.5,0.925),method="L-BFGS-B")


max[1]<--0.1646701
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05363331
max[5]<-0.02982051
max[6]<-0.03182533
max[7]<-10.14239
max[8]<-9.189999
max[9]<-0.8937402
max[10]<-0.9189999
max[11]<-0.8504374
max[12]<-0.3092366


maximo<--aux_max_2$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],x[1],max[6],x[2],max[8],max[9],max[10],x[3],max[12],index),-10^7)
}

aux_max_3<-opm(c(max[5],max[7],max[11]),fn=log_likelihood_prac, lower=c(0.027,9.5,0.8), upper=c(0.035,11.2,0.88),method="L-BFGS-B")

max[1]<--0.1646701
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05363331
max[5]<-0.03048197
max[6]<-0.03182533
max[7]<-10.15364
max[8]<-9.189999
max[9]<-0.8937402
max[10]<-0.9189999
max[11]<-0.8493179
max[12]<-0.3092366

maximo<--aux_max_3$value







log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],x[1],x[2],max[9],max[10],max[11],max[12],index),-10^7)
}

aux_max_4<-opm(c(max[7],max[8]),fn=log_likelihood_prac, lower=c(9.5,8.5), upper=c(11.2,10.5),method="L-BFGS-B")


max[1]<--0.1646701
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05363331
max[5]<-0.03048197
max[6]<-0.03182533
max[7]<-10.29486
max[8]<-9.452026
max[9]<-0.8937402
max[10]<-0.9189999
max[11]<-0.8493179
max[12]<-0.3092366


maximo<--aux_max_4$value






log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],x[1],x[2],max[6],x[3],x[4],max[9],x[5],x[6],max[12],index),-10^7)
}

aux_max_5<-opm(c(max[4],max[5],max[7],max[8],max[10],max[11]),fn=log_likelihood_prac, lower=c(0.047,0.027,9.5,8.5,0.89,0.8), upper=c(0.062,0.036,11.2,10.5,0.94,0.88),method="L-BFGS-B")


max[1]<--0.1646701
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05558302
max[5]<-0.03184931
max[6]<-0.03182533
max[7]<-10.29501
max[8]<-9.452318
max[9]<-0.8937402
max[10]<-0.9227176
max[11]<-0.8561966
max[12]<-0.3092366


maximo<--aux_max_5$value




log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],max[2],max[3],max[4],max[5],max[6],x[2],x[3],max[9],max[10],max[11],max[12],index),-10^7)
}

aux_max_6<-opm(c(max[1],max[7],max[8]),fn=log_likelihood_prac, lower=c(-0.30,9.5,8.5), upper=c(-0.08,11.2,10.5),method="L-BFGS-B")


max[1]<--0.171322
max[2]<--0.01833508
max[3]<-0.03301999
max[4]<-0.05558302
max[5]<-0.03184931
max[6]<-0.03182533
max[7]<-10.44006
max[8]<-9.718063
max[9]<-0.8937402
max[10]<-0.9227176
max[11]<-0.8561966
max[12]<-0.3092366


maximo<--aux_max_6$value





log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}

aux_max_7<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=c(0,0.1,0.15,0.065,0.038,0.037,12,11,1.2,0.94,0.88,0.47),method="L-BFGS-B")


max[1]<--0.171376
max[2]<--0.01704026
max[3]<-0.03436776
max[4]<-0.05617801
max[5]<-0.03256014
max[6]<-0.03186451
max[7]<-10.44024
max[8]<-9.718014
max[9]<-0.8971786
max[10]<-0.9223305
max[11]<-0.8556108
max[12]<-0.3082678


maximo<--aux_max_7$value





log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],x[1],x[2],max[9],max[10],max[11],max[12],index),-10^7)
}

aux_max_8<-opm(c(max[7],max[8]),fn=log_likelihood_prac, lower=c(9.5,8.5), upper=c(12,11),method="L-BFGS-B")


max[1]<--0.171376
max[2]<--0.01704026
max[3]<-0.03436776
max[4]<-0.05617801
max[5]<-0.03256014
max[6]<-0.03186451
max[7]<-10.60129
max[8]<-9.9671
max[9]<-0.8971786
max[10]<-0.9223305
max[11]<-0.8556108
max[12]<-0.3082678


maximo<--aux_max_8$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}

aux_max_9<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=c(0,0.1,0.15,0.065,0.038,0.037,12,11,1.2,0.94,0.88,0.47),method="L-BFGS-B")


max[1]<--0.1837625
max[2]<--0.01857889
max[3]<-0.03532642
max[4]<-0.05837285
max[5]<-0.03363158
max[6]<-0.0318744
max[7]<-10.60578
max[8]<-9.975434
max[9]<-0.899218
max[10]<-0.9261341
max[11]<-0.85837
max[12]<-0.309301


maximo<--aux_max_9$value



log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],x[1],x[2],max[9],max[10],max[11],max[12],index),-10^7)
}

aux_max_10<-opm(c(max[7],max[8]),fn=log_likelihood_prac, lower=c(9.5,8.5), upper=c(12,11),method="L-BFGS-B")

max[1]<--0.1837625
max[2]<--0.01857889
max[3]<-0.03532642
max[4]<-0.05837285
max[5]<-0.03363158
max[6]<-0.0318744
max[7]<-10.75103
max[8]<-10.24041
max[9]<-0.899218
max[10]<-0.9261341
max[11]<-0.85837
max[12]<-0.309301


maximo<--aux_max_10$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}

aux_max_11<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=c(0,0.1,0.15,0.065,0.038,0.037,12,11,1.2,0.94,0.88,0.47),method="L-BFGS-B")


max[1]<--0.1839092
max[2]<--0.01868052
max[3]<-0.03530391
max[4]<-0.05903609
max[5]<-0.03439293
max[6]<-0.03185918
max[7]<-10.75105
max[8]<-10.24044
max[9]<-0.8991071
max[10]<-0.9259927
max[11]<-0.8584903
max[12]<-0.3093874


maximo<--aux_max_11$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],x[1],x[2],max[9],max[10],max[11],max[12],index),-10^7)
}

aux_max_12<-opm(c(max[7],max[8]),fn=log_likelihood_prac, lower=c(9.5,8.5), upper=c(12,11),method="L-BFGS-B")



max[1]<--0.1839092
max[2]<--0.01868052
max[3]<-0.03530391
max[4]<-0.05903609
max[5]<-0.03439293
max[6]<-0.03185918
max[7]<-10.902
max[8]<-10.48451
max[9]<-0.8991071
max[10]<-0.9259927
max[11]<-0.8584903
max[12]<-0.3093874


maximo<--aux_max_12$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}

aux_max_13<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=c(0,0.1,0.15,0.065,0.038,0.037,12,11,1.2,0.94,0.88,0.47),method="L-BFGS-B")


max[1]<--0.1861503
max[2]<--0.01881601
max[3]<-0.03537574
max[4]<-0.06062584
max[5]<-0.03566942
max[6]<-0.03186014
max[7]<-10.90253
max[8]<-10.48546
max[9]<-0.8980585
max[10]<-0.9283528
max[11]<-0.8631826
max[12]<-0.3092829


maximo<--aux_max_13$value


log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],max[2],max[3],x[2],x[3],max[6],x[4],x[5],max[9],x[6],x[7],max[12],index),-10^7)
}

aux_max_14<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=
c(0,0.1,0.15,0.7,0.7,0.037,90,140,1.2,0.98,0.95,0.47)
,method="L-BFGS-B")

max[1]<--0.2335114
max[2]<--0.01881601
max[3]<-0.03537574
max[4]<-0.4357578
max[5]<-0.433507
max[6]<-0.03186014
max[7]<-66.12875
max[8]<-120
max[9]<-0.8980585
max[10]<-0.9557297
max[11]<-0.8989568
max[12]<-0.3092829


maximo<--aux_max_14$value



log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}


aux_max_15<-opm(c(max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12]),fn=log_likelihood_prac, lower=c(-0.30,-0.1,-0.05,0.05,0.029,0.028,9,8.5,0.6,0.89,0.8,0.22), upper=
c(0,0.1,0.15,0.7,0.7,0.037,90,125,1.2,0.98,0.95,0.47)
,method="L-BFGS-B")

max[1]<--0.2334785
max[2]<--0.01748831
max[3]<-0.03435184
max[4]<-0.4356599
max[5]<-0.4335856
max[6]<-0.03190299
max[7]<-66.12879
max[8]<-120
max[9]<-0.9000099
max[10]<-0.9557929
max[11]<-0.8988601
max[12]<-0.3082733


maximo<--aux_max_15$value






log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],x,max[9],max[10],max[11],max[12],index),-10^7)
}


aux_max_16<-opm(c(max[8]),fn=log_likelihood_prac, lower=9, upper=140,method="L-BFGS-B")

max<-rep(0,12)
max[1]<--0.2334785
max[2]<--0.01748831
max[3]<-0.03435184
max[4]<-0.4356599
max[5]<-0.4335856
max[6]<-0.03190299
max[7]<-66.12879
max[8]<-120.1132
max[9]<-0.9000099
max[10]<-0.9557929
max[11]<-0.8988601
max[12]<-0.3082733


log_likelihood_prac(max)

maximo<-log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)





#-----Esto no va en la maximizacion

#------------Change to lon data---------------
log_likelihood_prac<-function(x){
-max(log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],index),-10^7)
}




maximos_f<-rep(0,12)
maximos_f<-opm(max,fn=log_likelihood_prac, lower=c(-0.5,-0.5,-0.5,0.005,0.005,0.005,40,100,0.1,0.6,0.6,0.1), upper=c(0.3,0.3,0.3,0.7,0.7,0.4,80,160,2,0.97,0.92,0.4),method="L-BFGS-B")

#-----------------






#--------------------Profile likelihood---------

#------------Change to lon data---------------
log_likelihood_1<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),x,max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}



log_likelihood_2<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],x,max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_3<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],x,max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_4<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],x,max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_5<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],x,max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_6<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],x,max[7],max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_7<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],x,max[8],max[9],max[10],max[11],max[12],index)
}
log_likelihood_8<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],x,max[9],max[10],max[11],max[12],index)
}
log_likelihood_9<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],x,max[10],max[11],max[12],index)
}
log_likelihood_10<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],x,max[11],max[12],index)
}
log_likelihood_11<-function(x){
log_likelihood_3d_prac(T,n,c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],x,max[12],index)
}
log_likelihood_12<-function(x){
log_likelihood_3d_prac(T=T,n=n,telemetry=c(mu_1,mu_2,mu_3_r),max[1],max[2],max[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],x,index)
}




par(mfrow=c(4,3))

max




val_right<-rep(0,12)
val_left<-rep(0,12)

#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.5,0.1,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_1(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_1,max[1]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[1],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[1],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[1],sqrt(norm_1))
val1<-dnorm(valc1_1,max[1],sqrt(norm_1))

val_left[1]<-valc1_1
val_right[1]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[12])),xlab=expression(z[12]),ylab=expression(paste("Profile Likelihood ",z[12])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[1],max[1]),c(0,6),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,6),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,6),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------

#-----------rho_13 #------(-0.1,0.1) 
m<-50
reg_aux<-seq(-0.20,0.18,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_2(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_2,max[2]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[2],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[2],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[2],sqrt(norm_1))
val1<-dnorm(valc1_1,max[2],sqrt(norm_1))

val_left[2]<-valc1_1
val_right[2]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[13])),xlab=expression(z[13]),ylab=expression(paste("Profile Likelihood ",z[13])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[2],max[2]),c(0,10),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,10),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,10),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)



#----------



#-----------rho_23 #----(-0.05,0.15)
m<-50
reg_aux<-seq(-0.15,0.25,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_3(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_3,max[3]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[3],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[3],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[3],sqrt(norm_1))
val1<-dnorm(valc1_1,max[3],sqrt(norm_1))

val_left[3]<-valc1_1
val_right[3]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",z[23])),xlab=expression(z[23]),ylab=expression(paste("Profile Likelihood ",z[23])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[3],max[3]),c(0,10),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,10),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,10),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------



#-----------s_1  #----(0.05,0.065)
m<-50
reg_aux<-seq(0.35,0.55,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_4(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_4,max[4]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[4],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[4],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[4],sqrt(norm_1))
val1<-dnorm(valc1_1,max[4],sqrt(norm_1))

val_left[4]<-valc1_1
val_right[4]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[1])),xlab=expression(sigma[1]),ylab=expression(paste("Profile Likelihood ",sigma[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[4],max[4]),c(0,20),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,20),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,20),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------




#-----------s_2  #------(0.029,0.038)
m<-50
reg_aux<-seq(0.35,0.55,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_5(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_5,max[5]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[5],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[5],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[5],sqrt(norm_1))
val1<-dnorm(valc1_1,max[5],sqrt(norm_1))

val_left[5]<-valc1_1
val_right[5]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[2])),xlab=expression(sigma[2]),ylab=expression(paste("Profile Likelihood ",sigma[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[5],max[5]),c(0,20),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,20),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,20),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------s_3   #------(0.028,0.037)
m<-50
reg_aux<-seq(0.025,0.04,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_6(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_6,max[6]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[6],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[6],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[6],sqrt(norm_1))
val1<-dnorm(valc1_1,max[6],sqrt(norm_1))

val_left[6]<-valc1_1
val_right[6]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",sigma[3])),xlab=expression(sigma[3]),ylab=expression(paste("Profile Likelihood ",sigma[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[6],max[6]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----------


#-----------b_1  #----------(9,12)
m<-50
reg_aux<-seq(52,80,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_7(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_7,max[7]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[7],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[7],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[7],sqrt(norm_1))
val1<-dnorm(valc1_1,max[7],sqrt(norm_1))

val_left[7]<-valc1_1
val_right[7]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[1])),xlab=expression(beta[1]),ylab=expression(paste("Profile Likelihood ",beta[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[7],max[7]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------



#-----------b_2  #-----(8.5,11)    
m<-50
reg_aux<-seq(95,145,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_8(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_8,max[8]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[8],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[8],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[8],sqrt(norm_1))
val1<-dnorm(valc1_1,max[8],sqrt(norm_1))

val_left[8]<-valc1_1
val_right[8]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[2])),xlab=expression(beta[2]),ylab=expression(paste("Profile Likelihood ",beta[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[8],max[8]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------



#-----------b_3  #-----(0.6,1.2)
m<-50
reg_aux<-seq(0.3,1.45,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_9(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_9,max[9]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[9],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[9],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[9],sqrt(norm_1))
val1<-dnorm(valc1_1,max[9],sqrt(norm_1))

val_left[9]<-valc1_1
val_right[9]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",beta[3])),xlab=expression(beta[3]),ylab=expression(paste("Profile Likelihood ",beta[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[9],max[9]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_1   #--(0.89,0.94)
m<-50
reg_aux<-seq(0.93,0.98,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_10(reg_aux[i])-maximo
}


hessian(log_likelihood_10,max[10])


#norm_1<-1/(-hessian(log_likelihood_10,max[10]))
#norm_1<-norm_1[1,1]

norm_1<-0.00003
normal_approx_1<-dnorm(reg_aux,max[10],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[10],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[10],sqrt(norm_1))
val1<-dnorm(valc1_1,max[10],sqrt(norm_1))

val_left[10]<-valc1_1
val_right[10]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[1])),xlab=expression(H[1]),ylab=expression(paste("Profile Likelihood ",H[1])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[10],max[10]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_2   #---(0.8,0.88)
m<-50
reg_aux<-seq(0.84,0.95,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_11(reg_aux[i])-maximo
}

norm_1<-1/(-hessian(log_likelihood_11,max[11]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[11],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[11],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[11],sqrt(norm_1))
val1<-dnorm(valc1_1,max[11],sqrt(norm_1))

val_left[11]<-valc1_1
val_right[11]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[2])),xlab=expression(H[2]),ylab=expression(paste("Profile Likelihood ",H[2])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[11],max[11]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#----------


#-----------H_3  #------(0.22,0.47)
m<-50
reg_aux<-seq(0.15,0.54,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_12(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_12,max[12]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,max[12],sqrt(norm_1))

valc1_1<-qnorm(0.005,max[12],sqrt(norm_1))
valc2_1<-qnorm(0.995,max[12],sqrt(norm_1))
val1<-dnorm(valc1_1,max[12],sqrt(norm_1))

val_left[12]<-valc1_1
val_right[12]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",H[3])),xlab=expression(H[3]),ylab=expression(paste("Profile Likelihood ",H[3])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(max[12],max[12]),c(0,250),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,250),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,250),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)

#---------





rho<-coeff(max[1],max[2],max[3],max[10],max[11],max[12])



log_likelihood_r_12<-function(x){
log_likelihood_3d_prac_r(T=T,n=n,telemetry=c(mu_1,mu_2,mu_3_r),x,rho[2],rho[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}

log_likelihood_r_13<-function(x){
log_likelihood_3d_prac_r(T=T,n=n,telemetry=c(mu_1,mu_2,mu_3_r),rho[1],x,rho[3],max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}


log_likelihood_r_23<-function(x){
log_likelihood_3d_prac_r(T=T,n=n,telemetry=c(mu_1,mu_2,mu_3_r),rho[1],rho[2],x,max[4],max[5],max[6],max[7],max[8],max[9],max[10],max[11],max[12],index)
}


maximo



par(mfrow=c(1,3))
val_left_r<-rep(0,3)
val_right_r<-rep(0,3)



#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.5,0.1,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_12(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_r_12,rho[1]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,rho[1],sqrt(norm_1))

valc1_1<-qnorm(0.005,rho[1],sqrt(norm_1))
valc2_1<-qnorm(0.995,rho[1],sqrt(norm_1))
val1<-dnorm(valc1_1,rho[1],sqrt(norm_1))



val_left_r[1]<-valc1_1
val_right_r[1]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",rho[12])),xlab=expression(rho[12]),ylab=expression(paste("Profile Likelihood ",rho[12])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(rho[1],rho[1]),c(0,6),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,6),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,6),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)





#-----------rho_13 #------(-0.1,0.1) 
m<-50
reg_aux<-seq(-0.08,0.08,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_13(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_r_13,rho[2]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,rho[2],sqrt(norm_1))

valc1_1<-qnorm(0.005,rho[2],sqrt(norm_1))
valc2_1<-qnorm(0.995,rho[2],sqrt(norm_1))
val1<-dnorm(valc1_1,rho[2],sqrt(norm_1))

val_left_r[2]<-valc1_1
val_right_r[2]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",rho[13])),xlab=expression(rho[13]),ylab=expression(paste("Profile Likelihood ",rho[13])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(rho[2],rho[2]),c(0,25),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,25),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,25),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)



#----------



#-----------rho_23 #----(-0.05,0.15)
m<-50
reg_aux<-seq(-0.1,0.12,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_23(reg_aux[i])-maximo
}


norm_1<-1/(-hessian(log_likelihood_r_23,rho[3]))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,rho[3],sqrt(norm_1))

valc1_1<-qnorm(0.005,rho[3],sqrt(norm_1))
valc2_1<-qnorm(0.995,rho[3],sqrt(norm_1))
val1<-dnorm(valc1_1,rho[3],sqrt(norm_1))

val_left_r[3]<-valc1_1
val_right_r[3]<-valc2_1


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood ",rho[23])),xlab=expression(rho[23]),ylab=expression(paste("Profile Likelihood ",rho[23])),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(rho[3],rho[3]),c(0,25),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,25),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,25),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


val_left_r
val_right_r
rho


#----------


#-----Chi aproximation----------------#




par(mfrow=c(1,3))

f<-qchisq(0.99,1)

root_1_l_r<-rep(0,12)
root_1_r_r<-rep(0,12)



#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.45,0.05,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_12(reg_aux[i])-maximo
}





like_ratio_r_12<-function(x){
	-2*(log_likelihood_r_12(x)-maximo)-f
}


root_1_l_r[1]<-uniroot(like_ratio_r_12,lower=-0.5,upper=rho[1],tol=1e-9)$root
root_1_r_r[1]<-uniroot(like_ratio_r_12,lower=rho[1],upper=0.1,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",rho[12])),xlab=expression(rho[12]),ylab=expression(paste("Likelihood Ratio Statistic ",rho[12])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.5,0.1),c(f,f),col="blue",lwd=2)
lines(c(rho[1],rho[1]),c(0,15),col="red",lwd=2)
lines(c(root_1_l_r[1],root_1_l_r[1]),c(0,15),col="green",lwd=2)
lines(c(root_1_r_r[1],root_1_r_r[1]),c(0,15),col="green",lwd=2)

#-------------------------------------------------








#-----------rho_13 #------(-0.1,0.1) 
m<-50
reg_aux<-seq(-0.065,0.05,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_13(reg_aux[i])-maximo
}




like_ratio_r_13<-function(x){
	-2*(log_likelihood_r_13(x)-maximo)-f
}


root_1_l_r[2]<-uniroot(like_ratio_r_13,lower=-0.1,upper=rho[2],tol=1e-9)$root
root_1_r_r[2]<-uniroot(like_ratio_r_13,lower=rho[2],upper=0.1,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",rho[13])),xlab=expression(rho[13]),ylab=expression(paste("Likelihood Ratio Statistic ",rho[13])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.1,0.1),c(f,f),col="blue",lwd=2)
lines(c(rho[2],rho[2]),c(0,15),col="red",lwd=2)
lines(c(root_1_l_r[2],root_1_l_r[2]),c(0,15),col="green",lwd=2)
lines(c(root_1_r_r[2],root_1_r_r[2]),c(0,15),col="green",lwd=2)


#----------



#-----------rho_23 #----(-0.05,0.15)
m<-50
reg_aux<-seq(-0.06,0.1,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_r_23(reg_aux[i])-maximo
}




like_ratio_r_23<-function(x){
	-2*(log_likelihood_r_23(x)-maximo)-f
}


root_1_l_r[3]<-uniroot(like_ratio_r_23,lower=-0.1,upper=rho[3],tol=1e-9)$root
root_1_r_r[3]<-uniroot(like_ratio_r_23,lower=rho[3],upper=0.1,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",rho[23])),xlab=expression(rho[23]),ylab=expression(paste("Likelihood Ratio Statistic ",rho[23])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.1,0.12),c(f,f),col="blue",lwd=2)
lines(c(rho[3],rho[3]),c(0,15),col="red",lwd=2)
lines(c(root_1_l_r[3],root_1_l_r[3]),c(0,15),col="green",lwd=2)
lines(c(root_1_r_r[3],root_1_r_r[3]),c(0,15),col="green",lwd=2)

root_1_l_r
root_1_r_r
rho


#----------










#-------------Chi general




par(mfrow=c(4,3))

f<-qchisq(0.99,1)

root_1_l<-rep(0,12)
root_1_r<-rep(0,12)


#-----------rho_12 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.45,0.05,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_1(reg_aux[i])-maximo
}


like_ratio_1<-function(x){
	-2*(log_likelihood_1(x)-maximo)-f
}


root_1_l[1]<-uniroot(like_ratio_1,lower=-0.5,upper=-0.3,tol=1e-9)$root
root_1_r[1]<-uniroot(like_ratio_1,lower=-0.1,upper=0.1,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[12])),xlab=expression(z[12]),ylab=expression(paste("Likelihood Ratio Statistic ",z[12])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.47,0.1),c(f,f),col="blue",lwd=2)
lines(c(max[1],max[1]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[1],root_1_l[1]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[1],root_1_r[1]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------rho_13 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.15,0.14,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_2(reg_aux[i])-maximo
}


like_ratio_2<-function(x){
	-2*(log_likelihood_2(x)-maximo)-f
}


root_1_l[2]<-uniroot(like_ratio_2,lower=-0.2,upper=-0.1,tol=1e-9)$root
root_1_r[2]<-uniroot(like_ratio_2,lower=0,upper=0.18,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[13])),xlab=expression(z[13]),ylab=expression(paste("Likelihood Ratio Statistic ",z[13])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.2,0.18),c(f,f),col="blue",lwd=2)
lines(c(max[2],max[2]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[2],root_1_l[2]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[2],root_1_r[2]),c(0,15),col="green",lwd=2)

#-------------------------------------------------



#-----------rho_23 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(-0.12,0.2,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_3(reg_aux[i])-maximo
}


like_ratio_3<-function(x){
	-2*(log_likelihood_3(x)-maximo)-f
}


root_1_l[3]<-uniroot(like_ratio_3,lower=-0.12,upper=-0.05,tol=1e-9)$root
root_1_r[3]<-uniroot(like_ratio_3,lower=0.1,upper=0.2,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",z[23])),xlab=expression(z[23]),ylab=expression(paste("Likelihood Ratio Statistic ",z[23])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(-0.2,0.18),c(f,f),col="blue",lwd=2)
lines(c(max[3],max[3]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[3],root_1_l[3]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[3],root_1_r[3]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------s_1 #-------------(-0.30,0) 
m<-20
reg_aux<-seq(0.37,0.53,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_4(reg_aux[i])-maximo
}


like_ratio_4<-function(x){
	-2*(log_likelihood_4(x)-maximo)-f
}


root_1_l[4]<-uniroot(like_ratio_4,lower=0.35,upper=0.4,tol=1e-9)$root
root_1_r[4]<-uniroot(like_ratio_4,lower=0.45,upper=0.55,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[1])),xlab=expression(sigma[1]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.35,0.55),c(f,f),col="blue",lwd=2)
lines(c(max[4],max[4]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[4],root_1_l[4]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[4],root_1_r[4]),c(0,15),col="green",lwd=2)

#-------------------------------------------------





#-----------s_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.365,0.52,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_5(reg_aux[i])-maximo
}


like_ratio_5<-function(x){
	-2*(log_likelihood_5(x)-maximo)-f
}


root_1_l[5]<-uniroot(like_ratio_5,lower=0.35,upper=0.4,tol=1e-9)$root
root_1_r[5]<-uniroot(like_ratio_5,lower=0.45,upper=0.55,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[2])),xlab=expression(sigma[2]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.35,0.55),c(f,f),col="blue",lwd=2)
lines(c(max[5],max[5]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[5],root_1_l[5]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[5],root_1_r[5]),c(0,15),col="green",lwd=2)

#-------------------------------------------------



#-----------s_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.027,0.038,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_6(reg_aux[i])-maximo
}


like_ratio_6<-function(x){
	-2*(log_likelihood_6(x)-maximo)-f
}


root_1_l[6]<-uniroot(like_ratio_6,lower=0.025,upper=max[6],tol=1e-9)$root
root_1_r[6]<-uniroot(like_ratio_6,lower=max[6],upper=0.04,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",sigma[3])),xlab=expression(sigma[3]),ylab=expression(paste("Likelihood Ratio Statistic ",sigma[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.027,0.038),c(f,f),col="blue",lwd=2)
lines(c(max[6],max[6]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[6],root_1_l[6]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[6],root_1_r[6]),c(0,15),col="green",lwd=2)

#-------------------------------------------------


#-----------b_1 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(55,77,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_7(reg_aux[i])-maximo
}


like_ratio_7<-function(x){
	-2*(log_likelihood_7(x)-maximo)-f
}


root_1_l[7]<-uniroot(like_ratio_7,lower=55,upper=max[7],tol=1e-9)$root
root_1_r[7]<-uniroot(like_ratio_7,lower=max[7],upper=80,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[1])),xlab=expression(beta[1]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(52,80),c(f,f),col="blue",lwd=2)
lines(c(max[7],max[7]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[7],root_1_l[7]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[7],root_1_r[7]),c(0,15),col="green",lwd=2)

#-------------------------------------------------


#-----------b_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(100,140,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_8(reg_aux[i])-maximo
}


like_ratio_8<-function(x){
	-2*(log_likelihood_8(x)-maximo)-f
}


root_1_l[8]<-uniroot(like_ratio_8,lower=95,upper=max[8],tol=1e-9)$root
root_1_r[8]<-uniroot(like_ratio_8,lower=max[8],upper=145,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[2])),xlab=expression(beta[2]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(95,145),c(f,f),col="blue",lwd=2)
lines(c(max[8],max[8]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[8],root_1_l[8]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[8],root_1_r[8]),c(0,15),col="green",lwd=2)

#-------------------------------------------------

#-----------b_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.52,1.25,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_9(reg_aux[i])-maximo
}


like_ratio_9<-function(x){
	-2*(log_likelihood_9(x)-maximo)-f
}


root_1_l[9]<-uniroot(like_ratio_9,lower=0.3,upper=max[9],tol=1e-9)$root
root_1_r[9]<-uniroot(like_ratio_9,lower=max[9],upper=1.45,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",beta[3])),xlab=expression(beta[3]),ylab=expression(paste("Likelihood Ratio Statistic ",beta[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.3,1.45),c(f,f),col="blue",lwd=2)
lines(c(max[9],max[9]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[9],root_1_l[9]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[9],root_1_r[9]),c(0,15),col="green",lwd=2)

#-------------------------------------------------





#-----------H_1 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.935,0.97,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_10(reg_aux[i])-maximo
}


like_ratio_10<-function(x){
	-2*(log_likelihood_10(x)-maximo)-f
}


root_1_l[10]<-uniroot(like_ratio_10,lower=0.93,upper=max[10],tol=1e-9)$root
root_1_r[10]<-uniroot(like_ratio_10,lower=max[10],upper=0.98,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[1])),xlab=expression(H[1]),ylab=expression(paste("Likelihood Ratio Statistic ",H[1])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.93,0.98),c(f,f),col="blue",lwd=2)
lines(c(max[10],max[10]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[10],root_1_l[10]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[10],root_1_r[10]),c(0,15),col="green",lwd=2)

#-------------------------------------------------



#-----------H_2 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.86,0.93,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_11(reg_aux[i])-maximo
}


like_ratio_11<-function(x){
	-2*(log_likelihood_11(x)-maximo)-f
}


root_1_l[11]<-uniroot(like_ratio_11,lower=0.84,upper=max[11],tol=1e-9)$root
root_1_r[11]<-uniroot(like_ratio_11,lower=max[11],upper=0.95,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[2])),xlab=expression(H[2]),ylab=expression(paste("Likelihood Ratio Statistic ",H[2])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.84,0.95),c(f,f),col="blue",lwd=2)
lines(c(max[11],max[11]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[11],root_1_l[11]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[11],root_1_r[11]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----------H_3 #-------------(-0.30,0) 
m<-50
reg_aux<-seq(0.19,0.48,length=m)
val_perfil_aux<-rep(0,m)

for(i in 1:m){
	val_perfil_aux[i]<-log_likelihood_12(reg_aux[i])-maximo
}


like_ratio_12<-function(x){
	-2*(log_likelihood_12(x)-maximo)-f
}


root_1_l[12]<-uniroot(like_ratio_12,lower=0.15,upper=max[12],tol=1e-9)$root
root_1_r[12]<-uniroot(like_ratio_12,lower=max[12],upper=0.54,tol=1e-9)$root

plot(reg_aux,-2*val_perfil_aux,type="l",main=expression(paste("Likelihood Ratio Statistic ",H[3])),xlab=expression(H[3]),ylab=expression(paste("Likelihood Ratio Statistic ",H[3])),lwd=2,cex.main=2,cex.lab=1.2)
lines(c(0.15,0.54),c(f,f),col="blue",lwd=2)
lines(c(max[12],max[12]),c(0,15),col="red",lwd=2)
lines(c(root_1_l[12],root_1_l[12]),c(0,15),col="green",lwd=2)
lines(c(root_1_r[12],root_1_r[12]),c(0,15),col="green",lwd=2)

#-------------------------------------------------




#-----Interval rho's Normal Aproximation----- 

#-----MLE rho's-----

coeff(max[1],max[2],max[3],max[10],max[11],max[12])



l_rho_12_n<-coeff(-0.4405,max[2],max[3],max[10],max[11],max[12])[1]
r_rho_12_n<-coeff(-0.0263,max[2],max[3],max[10],max[11],max[12])[1]


l_rho_13_n<-coeff(max[1],-0.1392,max[3],max[10],max[11],max[12])[2]
r_rho_13_n<-coeff(max[1],0.1043,max[3],max[10],max[11],max[12])[2]

l_rho_23_n<-coeff(max[1],max[2],-0.0895,max[10],max[11],max[12])[3]
r_rho_23_n<-coeff(max[1],max[2],0.1582,max[10],max[11],max[12])[3]


#----Intervalr rho_12 Normal Aproximation---

l_rho_12_n
r_rho_12_n


#----Intervalr rho_13 Normal Aproximation---


l_rho_13_n
r_rho_13_n


#----Intervalr rho_23 Normal Aproximation---


l_rho_23_n
r_rho_23_n





#-----Interval rho's Normal Aproximation----- 



l_rho_12_n<-coeff(-0.4405,max[2],max[3],max[10],max[11],max[12])[1]
r_rho_12_n<-coeff(-0.0263,max[2],max[3],max[10],max[11],max[12])[1]


l_rho_13_n<-coeff(max[1],-0.1392,max[3],max[10],max[11],max[12])[2]
r_rho_13_n<-coeff(max[1],0.1043,max[3],max[10],max[11],max[12])[2]

l_rho_23_n<-coeff(max[1],max[2],-0.0895,max[10],max[11],max[12])[3]
r_rho_23_n<-coeff(max[1],max[2],0.1582,max[10],max[11],max[12])[3]


#----Intervalr rho_12 Normal Aproximation---

l_rho_12_n
r_rho_12_n


#----Intervalr rho_13 Normal Aproximation---


l_rho_13_n
r_rho_13_n


#----Intervalr rho_23 Normal Aproximation---


l_rho_23_n
r_rho_23_n






mu_dat<-c(mu_1,mu_2,mu_3_r)


maximos_f

z<-c(maximos_f$p1,maximos_f$p2,maximos_f$p3)
s<-c(maximos_f$p4,maximos_f$p5,maximos_f$p6)
b<-c(maximos_f$p7,maximos_f$p8,maximos_f$p9)
h<-c(maximos_f$p10,maximos_f$p11,maximos_f$p12)

r<-coeff(z[1],z[2],z[3],h[1],h[2],h[3])

sim<-100



vel<-simulation_vel_3D(T,r,s,b,h,mu_dat,index,sim)





vel_mean<-rep(0,3*n)
for(i in 1:(3*n)){
	vel_mean[i]<-mean(vel[,i])
}


par(mfrow=c(3,2))



plot(time,mu_1,type="l",xlab="t",ylab="degree longitude",main="Trajectory longitude",lwd=2,cex.main=2,cex.lab=1.5,col="Blue")

plot(t,c(0,vel[1,1:n]),type="l",main="Velocity prediction longitude",xlab="t",ylab="degree longitude",lwd=2,cex.main=2,cex.lab=1.5,col="Gray",ylim=c(-0.035,0.02))
for(i in 2:sim){
	lines(t,c(0,vel[i,1:n]),col="Gray")
}


lines(t,c(0,vel_mean[1:n]),col="Green")


legend("bottomleft", legend=c("Predictive trajectory samples","Trajectory predicted mean"),
       col=c("Gray","Green"), lty=1, cex=1.1)


plot(time,mu_2,type="l",xlab="t",ylab="degree latitude",main="Trajectory latitude",lwd=2,cex.main=2,cex.lab=1.5,col="Blue")

plot(t,c(0,vel[1,(n+1):(2*n)]),type="l",main="Velocity prediction latitude",xlab="t",ylab="degree latitude",lwd=2,cex.main=2,cex.lab=1.5,col="Gray",ylim=c(-0.03,0.02))
for(i in 2:sim){
	lines(t,c(0,vel[i,(n+1):(2*n)]),lwd=2,cex.main=2,cex.lab=1.5,col="Gray")
}


lines(t,c(0,vel_mean[(n+1):(2*n)]),lwd=2,cex.main=2,cex.lab=1.5,col="Green")



legend("bottomleft", legend=c("Predictive trajectory samples","Trajectory predicted mean"),
       col=c("Gray","Green"), lty=1, cex=1.1)


plot(time,mu_3,type="l",xlab="t",ylab="metres above mean sea level",main="Trajectory altitude",lwd=2,cex.main=2,cex.lab=1.5,col="Blue")

plot(t,1000*c(0,vel[1,(2*n+1):(3*n)]),type="l",main="Velocity prediction altitude",xlab="t",ylab="meters above mean sea level",lwd=2,cex.main=2,cex.lab=1.5,col="Gray",ylim=c(-100,100))
for(i in 2:sim){
	lines(t,1000*c(0,vel[i,(2*n+1):(3*n)]),lwd=2,cex.main=2,cex.lab=1.5,col="Gray")
}


lines(t,1000*c(0,vel_mean[(2*n+1):(3*n)]),lwd=2,cex.main=2,cex.lab=1.5,col="Green")


#a<-(k/2)*c(1:(n/k))
 
legend("topleft", legend=c("Predictive trajectory samples","Trajectory predicted mean"),
       col=c("Gray","Green"), lty=1, cex=1.1)







#-----------------------Simulaciones 3D------------------------------------------------#

#-----Parametros movimiento 1
s_1<-1
b_1<-2
h_1<-0.9

mu_i_1<-mu_1[1]
v_i_1<-0
#------Parametros movimiento 2
s_2<-10
b_2<-2
h_2<-0.7

mu_i_2<-mu_2[1]
v_i_2<-0


#------Parametros movimiento 3
s_3<-max[6]
b_3<-max[9]
h_3<-max[12]
mu_i_3<-mu_3_r[1]
v_i_3<-0


z_12<-0.95
z_13<-max[2]
z_23<-max[3]



#-----Varibles globales
T<-120
n<-240
t<-c(0,(1:n)/2)
D<-T/n



telemetry<-simulation_3d(mu_i_1,mu_i_2,mu_i_3,z_12,z_13,z_23,s_1,s_2,s_3,b_1,b_2,b_3,h_1,h_2,h_3)

mu_1_s<-telemetry[1:(n+1)]
mu_2_s<-telemetry[(n+2):(2*n+2)]
mu_3_s<-telemetry[(2*n+3):(3*n+3)]



plot(mu_1_s,mu_2_s,type="l")
points(mu_1_s[1],mu_2_s[1],col="red")

par(mfrow=c(1,3))
plot(t,mu_1_s,type="l")
plot(t,mu_2_s,type="l")
plot(t,mu_3_s,type="l")
plot(time,mu_1,type="l",xlab="t",ylab="longitude",main="Trajectory Longitude",lwd=2,cex.main=2,cex.lab=1.5)
plot(time,mu_2,type="l",xlab="t",ylab="latitude",main="Trajectory Latitude",lwd=2,cex.main=2,cex.lab=1.5)
plot(time,mu_3_r,type="l",xlab="t",ylab="altitude",main="Trajectory Altitude",lwd=2,cex.main=2,cex.lab=1.5)




mfpca=function(Y_sparse,time_sparse,group,Q,nknots,nloop,burnin,thin,orth,sim){

	library(splines)
	library(MASS)
	library(Matrix)
	#library(MCMCpack)
	source("basis_setup_sparse.R")


	if(sim==TRUE){
		#rm(list=ls())
		load("sim_N100_Orth.RData")
		nknots_true=params[[1]]
		Q_true=params[[2]]
		R_true=params[[3]]
		D_true=params[[4]]
		Sigma_omega_true=params[[5]]
		Theta_true=params[[6]]
		Mu_true=params[[7]]
		Fpc_true=params[[7]]
		#time_sparse=TIME_SPARSE[[1]]
		#Y_sparse=Y_SPARSE[[1]]
		Alpha_true=ALPHA[[1]]
		F_true=F_SPARSE[[1]]
		Omega_true=OMEGA_SPARSE[[1]]
		#group=rep(1,N)
		#Q=Q_true
		#orth=FALSE
		#nloop=1110
		#burnin=10
		#thin=5
		est_mu=1
		est_alpha=1
		est_d=1
		est_r=1
		est_sigma_omega=1
    	est_theta=1
	}		
	if(sim==FALSE){
		est_mu=1
		est_alpha=1
		est_d=1
		est_r=1
		est_sigma_omega=1
    	est_theta=1
	}
	N=length(Y_sparse)
	P=length(Y_sparse[[1]])

	basis_stuff=basis_setup_sparse(time_sparse,nknots,plotit=TRUE,orth=orth)
	knots=basis_stuff[[1]]
	phi_t=basis_stuff[[3]]
	time_cont=basis_stuff[[4]]
	phi_t_cont=basis_stuff[[5]]
	K=nknots+4

	P=length(Y_sparse[[1]]) # repetition: not necessary
	N=length(Y_sparse) # repetition: not necessary
	id=1:N
	G=length(unique(group))
	N_g=sum(group==unique(group)[1])
	if(G>1){
		for(g in 2:G){
			N_g=c(N_g,sum(group==unique(group)[g]))
		}
	}
                        
	## initialize parameters   

	Alpha=array(0,dim=c(sum(Q),N)) # FPCA scores for all N subjects
	Mu=cbind(array(0,dim=c(sum(K),1))) # Basis coefficients for mean curves
	Theta=array(0,dim=c(sum(K),sum(Q))) # FPCA loadings
	for(p in 1:P){
		Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]=mvrnorm(n=1,rep(0,K[p]),diag(rep(1,K[p]))) # multivariate normal
	    Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]=
	    		Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]/
	    		sqrt(sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),sum(Q[min(1,p-1):(p-1)])+1]^2))
	    if(Q[p]>1){
			for(q in (sum(Q[min(1,p-1):(p-1)])+2):sum(Q[1:p])){
				Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=mvrnorm(n=1,rep(0,K[p]),diag(rep(1,K[p])))
				for(r in 1:(q-sum(Q[min(1,p-1):(p-1)])-1)){
					Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]-
						sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]*Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r])/
						sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r]^2)*Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q-r]
	    	            Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]=Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]/
	    	            		sqrt(sum(Theta[(sum(K[min(1,p-1):(p-1)])+1):(sum(K[1:p])),q]^2))
	    	    }
			} 
		}
	}		
	if(G==1){ # Correlation of FPCA scores
		R=diag(rep(1,sum(Q)))
	}	
	if(G>1){
		R=array(0,dim=c(sum(Q),sum(Q),G)) 
		for(g in 1:G){
			R[,,g]=diag(rep(1,sum(Q)))
		}
	}
	D=diag(rep(1,sum(Q))) #Variance of FPCA scores
	Sigma_alpha=D%*%R%*%D #Variance-Covariance of FPCA scores
	#phi=array(0,dim=c(sum(Q)))
	#delta=array(1,dim=c(sum(Q)))
	#nu=0.1
	#mu_i=array(0,dim=c(N))
	#sigma_i=0.1
	#mu_g=array(0,dim=c(G))
	#sigma=0.1
	Sigma_omega=0.1*rep(1,P) # Variances of errors
	#eta_eps=1
	#S_eps=.1*diag(rep(1,P))
	c_omega=1000
	psi=0
	nu=1000^.5
	delta=0.01^.5
	eta1_eps=.01
	eta2_eps=.01
	array_ind=0

    Alpha_mean=Alpha
    Mu_mean=Mu
    Theta_mean=Theta
    R_mean=R
    D_mean=D
    Sigma_omega_mean=Sigma_omega

			
	# set up arrays
	ALPHA_array=NULL
	MU_array=NULL
	D_array=NULL
	R_array=NULL
	THETA_array=NULL
	SIGMA_OMEGA_array=NULL
	
	# run MCMC 
   	for(iter in 1:nloop){

		print(iter)

		#draw Mu
		if(est_mu==0){
			Mu=Mu_true
		}
		if(est_mu==1){
			Sigma_mu=1/c_omega*diag(rep(1,sum(K)))
			mu_mu=rep(0,sum(K))
			#for(i in 1:N){
			#	Y_i=Y_sparse[[i]]
			#	times_i=time_sparse[[i]]
			#	Y_i_vec=Y_i[[1]]
			#	Phi_D_t_i=t(phi_t[[1]][[i]])
			#	Sigma_omega_i=Sigma_omega[1]*diag(length(Y_i[[1]]))
			#	for(p in 2:P){
			#		Y_i_vec=c(Y_i_vec,Y_i[[p]])
			#		Phi_D_t_i=bdiag(Phi_D_t_i,t(phi_t[[p]][[i]]))
			#		Sigma_omega_i=bdiag(Sigma_omega_i,Sigma_omega[p]*diag(length(Y_i[[p]])))
			#	}
			#	Y_i_vec_resid=Y_i_vec-Phi_D_t_i%*%Theta%*%cbind(Alpha[,i])
			#	Sigma_mu=Sigma_mu+t(Phi_D_t_i)%*%solve(Sigma_omega_i)%*%Phi_D_t_i
			#	mu_mu=mu_mu+t(Phi_D_t_i)%*%solve(Sigma_omega_i)%*%Y_i_vec_resid
			#}
			for(i in 1:N){
				Y_i=Y_sparse[[i]]
				times_i=time_sparse[[i]]
				Y_i_vec=Y_i[[1]]
				Phi_D_t_i=t(phi_t[[1]][[i]])
				Sigma_omega_i=Sigma_omega[1]*diag(length(Y_i[[1]]))
				for(p in 2:P){
					Y_i_vec=c(Y_i_vec,Y_i[[p]])
					Phi_D_t_i=bdiag(Phi_D_t_i,t(phi_t[[p]][[i]]))
					Sigma_omega_i=bdiag(Sigma_omega_i,Sigma_omega[p]*diag(length(Y_i[[p]])))
				}
				Sigma_mu=Sigma_mu+t(Phi_D_t_i)%*%solve((Phi_D_t_i%*%Theta%*%Sigma_alpha%*%t(Theta)%*%t(Phi_D_t_i)+Sigma_omega_i))%*%Phi_D_t_i
				mu_mu=mu_mu+t(Phi_D_t_i)%*%solve((Phi_D_t_i%*%Theta%*%Sigma_alpha%*%t(Theta)%*%t(Phi_D_t_i)+Sigma_omega_i))%*%cbind(Y_i_vec)
			}
			Sigma_mu=.5*(Sigma_mu+t(Sigma_mu))
			Sigma_mu=solve(Sigma_mu)
			mu_mu=Sigma_mu%*%mu_mu
			Mu=as.vector(mvrnorm(n=1,mu_mu,Sigma_mu))
		}


		#draw ALPHA
		if(est_alpha==0){
			Alpha=Alpha_true
		}
		if(est_alpha==1){
			for(i in 1:N){
				Y_i=Y_sparse[[i]]
				times_i=time_sparse[[i]]
				Y_i_vec=Y_i[[1]]
				Phi_D_t_i=t(phi_t[[1]][[i]])
				Sigma_omega_i=Sigma_omega[1]*diag(rep(1,length(times_i[[1]])))
				for(p in 2:P){
					Y_i_vec=c(Y_i_vec,Y_i[[p]])
					Phi_D_t_i=bdiag(Phi_D_t_i,t(phi_t[[p]][[i]]))
					Sigma_omega_i=bdiag(Sigma_omega_i,Sigma_omega[p]*diag(rep(1,length(times_i[[p]]))))
				}
				Y_i_vec_resid=Y_i_vec-Phi_D_t_i%*%Mu
				Sigma_alpha_i=t(Theta)%*%t(Phi_D_t_i)%*%solve(Sigma_omega_i)%*%Phi_D_t_i%*%Theta+2*solve(Sigma_alpha)
				Sigma_alpha_i=.5*(Sigma_alpha_i+t(Sigma_alpha_i))
				Sigma_alpha_i=solve(Sigma_alpha_i)
				mu_alpha_i=t(Theta)%*%t(Phi_D_t_i)%*%solve(Sigma_omega_i)%*%(Y_i_vec_resid)
				mu_alpha_i=Sigma_alpha_i%*%mu_alpha_i
	    		    Alpha[,i]=as.vector(mvrnorm(n=1,mu_alpha_i,Sigma_alpha_i))
			}
		}
		
				
		#draw D
		if(est_d==0){
			D=D_true
		}	
		if(est_d==1){
			C=solve(R)
			s_y=Alpha%*%t(Alpha)
			for(q in 1:sum(Q)){
				a=0
				if(q==1){
					for(r in (q+1):sum(Q)){
						a=a+s_y[q,r]*C[q,r]/D[r,r]
					}	
				}
				if(q>1 & q<sum(Q)){
					for(r in c(1:(q-1),(q+1):sum(Q))){
						a=a+s_y[q,r]*C[q,r]/D[r,r]
					}
				}
				if(q==sum(Q)){
					for(r in 1:(q-1)){
						a=a+s_y[q,r]*C[q,r]/D[r,r]
					}
				}
				ran=1
				a_d_p=seq(0.01,ran,0.01)
				c_temp=0
				w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
				if(max(w_d_p)<1e-100){
					c_temp=max(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p)
					w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
				}
				if (max(w_d_p)>1e+100){
					c_temp=max(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p)
					w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
				}
	 			while(w_d_p[length(w_d_p)]>max(w_d_p)*0.05){
					ran=ran+1
					a_d_p=seq(0.01,ran,0.01)
					c_temp=0
					w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
					if(max(w_d_p)<1e-100){
						c_temp=max(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p)
						w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*
							C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
					}
					if (max(w_d_p)>1e+100){
						c_temp=max(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*C[q,q]/(a_d_p^2)/2-a/a_d_p)
						w_d_p=exp(-log(a_d_p)*(N+1)-(log(a_d_p)-psi)^2/(delta^2)/2-s_y[q,q]*
						 C[q,q]/(a_d_p^2)/2-a/a_d_p-c_temp)
					}
				}
				f_d_p=w_d_p/sum(w_d_p*0.01)
				sum_f_d_p=0*f_d_p
				for(l in 1:length(f_d_p)){
					sum_f_d_p[l]=sum(f_d_p[1:l])*0.01
				}
				alpha_d_p=runif(1,0,1)
				pos=(1:length(sum_f_d_p))[sum_f_d_p>=alpha_d_p]
				pos=pos[1]
				if(pos>1){
					D[q,q]=a_d_p[pos]+(alpha_d_p-sum_f_d_p[pos-1])/(sum_f_d_p[pos]-sum_f_d_p[pos-1])*0.01
				}
				if(pos==1){
					D[q,q]=a_d_p[pos]
				}
					
			}
		}

		
		#draw R
		if(est_r==0){
			R=R_true
		}
		if(est_r==1){
			R_N=R
			B=array(0,dim=c(sum(Q),sum(Q)))
			for(i in 1:N){
				B=B+solve(D)%*%Alpha[,i]%*%t(Alpha[,i])%*%t(solve(D))
			}
			for(p in 1:(P-1)){
				for(q1 in (sum(Q[min(1,p-1):(p-1)])+1):sum(Q[1:p])){
					for(q2 in (sum(Q[1:p])+1):sum(Q)){
						R_temp=R_N
						R_temp[q1,q2]=1
						R_temp[q2,q1]=1
						f_1=det(R_temp)
						R_temp[q1,q2]=-1
						R_temp[q2,q1]=-1
						f_minus1=det(R_temp)
						R_temp[q1,q2]=0
						R_temp[q2,q1]=0
						f_0=det(R_temp)
						a=(f_1+f_minus1-2*f_0)/2
						b=(f_1-f_minus1)/2
						c=f_0
						if(a<0){
							root1=(-b+sqrt(b^2-4*a*c))/2/a
							root2=(-b-sqrt(b^2-4*a*c))/2/a
							l_q1q2=max(root1,-1)
							l_q1q2=ceiling(l_q1q2*100)/100
							u_q1q2=min(root2,1)
							u_q1q2=floor(u_q1q2*100)/100
							r_q1q2_N=l_q1q2+(u_q1q2-l_q1q2)*runif(1,0,1)
						}
						if(a>0){
							if ((b^2-4*a*c)<0){
								r_q1q2_N=-1+2*runif(1,0,1)
							}
							if ((b^2-4*a*c)>=0){
								root1=(-b-sqrt(b^2-4*a*c))/2/a
								root2=(-b+sqrt(b^2-4*a*c))/2/a
								l_q1q21=-1
								u_q1q21=max(root1,-1)
								u_q1q21=floor(u_q1q21*100)/100
								l_q1q22=min(root2,1)
								l_q1q22=ceiling(l_q1q22*100)/100
								u_q1q22=1
								alpha_temp=runif(1,0,1)
								if (alpha_temp <= (u_q1q21-l_q1q21)/(u_q1q21-l_q1q21+u_q1q22-l_q1q22)){
									r_q1q2_N=l_q1q21+(u_q1q21-l_q1q21+u_q1q22-l_q1q22)*alpha_temp
	                              }            
	                              if (alpha_temp > (u_q1q21-l_q1q21)/(u_q1q21-l_q1q21+u_q1q22-l_q1q22)){
									r_q1q2_N=l_q1q22+(u_q1q21-l_q1q21+u_q1q22-l_q1q22)*alpha_temp-(u_q1q21-l_q1q21)
								}
							}
						}
						R_temp_C=R_N
						R_temp_N=R_N
						R_temp_N[q1,q2]=r_q1q2_N
						R_temp_N[q2,q1]=R_temp_N[q1,q2]
						r_q1q2_C=R_temp_C[q1,q2]
						alpha_q1q2=runif(1,0,1)
						ratio_N_C=exp((log(det(R_temp_N))-log(det(R_temp_C)))*(-N/2)-1/2*sum(diag((solve(R_temp_N)%*%B)))+
							1/2*sum(diag(solve(R_temp_C)%*%B))-(r_q1q2_N)^2/(nu^2)/2+(r_q1q2_C)^2/(nu^2)/2);
						ratio_q1q2=min(1,ratio_N_C)
						if(ratio_q1q2>alpha_q1q2){
							R_N[q1,q2]=r_q1q2_N
							R_N[q2,q1]=R_N[q1,q2]
						}
					}                          
				}
			}
			R=R_N
			#print(R)
		}


		# Compute Sigma_alpha
		Sigma_alpha=D%*%R%*%D

	    #draw SIGMA_OMEGA
		if(est_sigma_omega==0){
			Sigma_omega=Sigma_omega_true
		}
		if(est_sigma_omega==1){
			indk=1
			indq=1
			for(p in 1:P){
				eta1_omega_p=eta1_eps
				eta2_omega_p=eta2_eps
				for(i in 1:N){
					Y_pi=Y_sparse[[i]][[p]]
					eta1_omega_p=eta1_omega_p+.5*length(Y_pi)
					Y_pi_resid=Y_pi-t(phi_t[[p]][[i]])%*%cbind(Mu[indk:sum(K[1:p])])-
						t(phi_t[[p]][[i]])%*%Theta[indk:sum(K[1:p]),indq:sum(Q[1:p])]%*%cbind(Alpha[indq:sum(Q[1:p]),i])
					eta2_omega_p=eta2_omega_p+.5*sum(Y_pi_resid^2)
				}
				Sigma_omega[p]=1/rgamma(1,eta1_omega_p,eta2_omega_p)	
				indk=indk+K[p]
				indq=indq+Q[p]
			}
		}


	    #draw THETA
		if(est_theta==0){
			#Theta=Theta_true
		}
	    if(est_theta==1){
	    	indk=1
			indq=1
			for(p in 1:P){
				mu_theta_p=array(0,dim=c(Q[p]*K[p],1))
				Sigma_theta_p=diag(rep(1,Q[p]*K[p]))/c_omega
				for(i in 1:N){
					Y_pi=Y_sparse[[i]][[p]]
					Alpha_pi=cbind(Alpha[indq:sum(Q[1:p]),i])
					for(t in 1:length(Y_pi)){
						Phi_D_pit=cbind(phi_t[[p]][[i]][,t])
						Y_pit_resid=Y_pi[t]-t(Phi_D_pit)%*%cbind(Mu[indk:sum(K[1:p])])
						if(Q[p]>1){
							for(q in 2:Q[p]){
								Phi_D_pit=bdiag(Phi_D_pit,cbind(phi_t[[p]][[i]][,t]))
							}
						}
						Sigma_theta_p=Sigma_theta_p+Phi_D_pit%*%Alpha_pi%*%t(Alpha_pi)%*%t(Phi_D_pit)/Sigma_omega[p]
						mu_theta_p=mu_theta_p+Phi_D_pit%*%Alpha_pi%*%cbind(Y_pit_resid)/Sigma_omega[p]
					}
				}
				Sigma_theta_p=.5*(Sigma_theta_p+t(Sigma_theta_p))
				Sigma_theta_p=solve(Sigma_theta_p)
				Sigma_theta_p=.5*(Sigma_theta_p+t(Sigma_theta_p))
				mu_theta_p=Sigma_theta_p%*%mu_theta_p
				Theta_p=mvrnorm(1,mu_theta_p,Sigma_theta_p)
				Theta[indk:sum(K[1:p]),indq:sum(Q[1:p])]=array(Theta_p,dim=c(K[p],Q[p]))

				poolvar=D[indq:sum(Q[1:p]),indq:sum(Q[1:p])]%*%D[indq:sum(Q[1:p]),indq:sum(Q[1:p])]
				temp_sigma=Theta[indk:sum(K[1:p]),indq:sum(Q[1:p])]%*%poolvar%*%t(Theta[indk:sum(K[1:p]),indq:sum(Q[1:p])])
				eigen_temp_sigma=eigen(temp_sigma)
				v_temp=eigen_temp_sigma$vectors
				d_temp=eigen_temp_sigma$values
				for(com in 1:length(d_temp)){
					if(!(d_temp[com]-Re(d_temp[com])==0)){
						d_temp[com]=-1*10^5
					}
				}
				pos_temp=array(0,dim=c(Q[p],1))
				for(pos in 1:Q[p]){
					pos_temp[pos]=(1:length(d_temp))[max(d_temp)==d_temp]
					d_temp[pos_temp[pos]]=-1e+5
				}
				Theta[indk:sum(K[1:p]),indq:sum(Q[1:p])]=v_temp[,pos_temp]
				for(q in 1:Q[p]){
					Theta[,sum(Q[min(1,p-1):((p-1))])+q]=sign(Theta[indk,sum(Q[min(1,p-1):((p-1))])+q])*
						Theta[,sum(Q[min(1,p-1):((p-1))])+q]
				}
				indk=indk+K[p]
				indq=indq+Q[p]
			}
			print(Theta)
		}


		if(iter%%thin==0 & iter>=burnin){
			array_ind=array_ind+1
			ALPHA_array[[array_ind]]=Alpha
			MU_array[[array_ind]]=Mu
			D_array[[array_ind]]=D
			R_array[[array_ind]]=R
			SIGMA_OMEGA_array[[array_ind]]=Sigma_omega
			THETA_array[[array_ind]]=Theta
			
			Alpha_mean=((array_ind-1)*Alpha_mean+Alpha)/array_ind
			Mu_mean=((array_ind-1)*Mu_mean+Mu)/array_ind
			D_mean=((array_ind-1)*D_mean+D)/array_ind
			R_mean=((array_ind-1)*R_mean+R)/array_ind
			Sigma_omega_mean=((array_ind-1)*Sigma_omega_mean+Sigma_omega)/array_ind
			Theta_mean=((array_ind-1)*Theta_mean+Theta)/array_ind
						
			if(sim==TRUE){
				print(cbind(Alpha_mean[,1],Alpha_true[,1]))
				print(cbind(Mu_mean,Mu_true))
				print(cbind(D_mean,D_true))
				print(cbind(R_mean,R_true))
				print(cbind(Sigma_omega_mean,Sigma_omega_true))
				print(cbind(Theta_mean,Theta_true))
			}	
			if(sim==FALSE){
				print(cbind(Mu_mean))
				print(cbind(Sigma_omega_mean))
				print(cbind(D_mean))
				print(cbind(R_mean))
				print(cbind(Theta_mean))
			}	
		}
	}	

	mcmc_results=list()
	mcmc_results[[1]]=ALPHA_array
	mcmc_results[[2]]=MU_array
	mcmc_results[[3]]=D_array
	mcmc_results[[4]]=R_array
	mcmc_results[[5]]=SIGMA_OMEGA_array
	mcmc_results[[6]]=THETA_array
	mcmc_results[[7]]=basis_stuff
	if(sim==TRUE){
		mcmc_results[[8]]=params
	}	

	return(mcmc_results)
}

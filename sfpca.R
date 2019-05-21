### prepare data for sfpca model
prepare_data = function(data, unique_subject_id, time_var, response, standardize.y=FALSE, scale.time=FALSE){
	# data: target longitudinal data for analysis (must be a data frame)
	# unique_subject_id: the column name corresponding to unique subject id in the data (string)
	# time_var: the column name corresponding to the time variable in the data (string)
	# response: the column name of the intersted response variable
	# standardize.y: the option of whether or not to standardize response variable (True/False) with mean 0 and sd 1
	# scale.time: the option of whether or not to scale the time variable to be within [0, 1] (True/False)

	# create new ID
	data$ID = as.numeric(as.numeric(data[, unique_subject_id]))
	N = length(unique(data$ID)) # total number of unique subjects

	# create time 
	if (scale.time == TRUE){
		data$time = (data[, time_var] - min(data[, time_var])) / (max(data[, time_var]) - min(data[, time_var]))
	} else{
		data$time = data[, time_var]
	}

	T = length(unique(data$time)) # total number of sampling time points

	# create response 
	if (standardize.y == TRUE){
		data$response = (data[, response] - mean(data[, response], na.rm=T)) / sd(data[, response], na.rm=T)
	} else {
		data$response = data[, response]

	}
	
	# re-order the data by ID and time
	data = data[order(data$ID, data$time), ]

	# create visits vector, response and time matrix
	ID.list = unique(data$ID)
	visits.vector = vector(mode = "numeric", length = N)
	response.list = NULL
	time.matrix = matrix(rep(0, N*T), nrow=N)

	# visits index for each individual when stacking the data
	visits.stop = vector(mode = "numeric", length = N)

	# size index for each individual in covariance matrix
	cov.start = cov.stop = vector(mode = "numeric", length = N)

	for(i in 1:N){ 
		# visits vector
		subject_i = data[data$ID==ID.list[i],]
		subject_i$n_visits = dim(subject_i)[1]	
		visits.vector[i] = unique(subject_i$n_visits)

		# visits index
		visits.stop[i] = sum(visits.vector)

		# covariance size index
		cov.stop[i] = sum(visits.vector^2)

	    # response matrix
	    # response.matrix[i, ] = c(subject_i$response, rep(0, T - unique(subject_i$n_visits)))
	    response.list = c(response.list, subject_i$response)

		# time matrix
		time.matrix[i, ] = c(subject_i$time, rep(0, T - unique(subject_i$n_visits)))

		rm(subject_i)
	}	
	visits.start = c(1, visits.stop[-N] + 1)
	cov.start = c(1, cov.stop[-N] + 1)
	cov.size = sum(visits.vector^2)

	prepared_data = list(data=data, num_subjects=N, num_times=T, response.list=response.list, time.matrix=time.matrix,
		                 visits.vector=visits.vector, visits.start=visits.start, visits.stop=visits.stop,
		                 cov.start=cov.start, cov.stop=cov.stop, cov.size=cov.size)
	return(prepared_data)
}


### set up spline basis for sparse data
basis_setup_sparse = function(prepared_data, nknots, orth=TRUE, delta=1/10000){
	# prepared_data: longitudinal data after applying prepared_data() function
	# knots: user-defined number of knots
	# orth: default setting for orth should be TRUE (after discussed with Wes on 02/13/2019)

	time_var = prepared_data$data$time
	num_subjects = prepared_data$num_subjects
	num_times = prepared_data$num_times
	S = prepared_data$time.matrix
	V = prepared_data$visits.vector

	# continuous time interval
	time_unique = sort(unique(time_var))
	time_min = min(time_unique)
	time_max = max(time_unique)
	time_cont = seq(time_min, time_max / delta) * delta # chop the entire time interval into many small subintervals
	time_cont = round(time_cont / delta)*delta # to avoid rounding error?

	# specify placement of knots
	qs = 1/(nknots + 1)
	knots = quantile(time_unique, qs)
	if(nknots > 1){
		for(q in 2:nknots){
			knots = c(knots, q*quantile(time_unique,qs))
		}
	}

	knots = as.vector(knots)


	# obtain cubic spline basis
	library('splines')

	#### stack subjects and visits
	## 1. for densely sampled time points
	phi_t_cont=list()
	phi_t_cont = bs(time_cont, knots=knots, degree=3,intercept=TRUE) # cubic spline, degree=spline_degree

	### the same as in setup_basis_sparse so far

	# Gram-Schmidt Orthonormalization
	temp = phi_t_cont
	K = nknots + 4 # num of spline basis 

	for(k in 1:K){
		if(orth==TRUE){
			if(k > 1){
				for(q in 1:(k-1)){
					temp[,k]=temp[,k]-(sum(temp[,k]*temp[,k-q])/
						sum(temp[,k-q]^2))*temp[,k-q];
				}
			}
		}		
	    temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
	}

	phi_t_cont=t(sqrt(1/delta)*temp)

	## 2. for sparsely sampled time points
	phi_t_stacked=NULL
	phi_t=list()
	for(i in 1:num_subjects){
		phi_t[[i]] = array(0,dim=c(K, V[i])) # phi_t: K (number of basis function) * number of total visit for each subject

		for(k in 1:K){
			for(t in 1:V[i]){
				phi_t[[i]][k, t] = phi_t_cont[k, abs(time_cont - S[i, t]) == min(abs(time_cont - S[i, t]))]
			}
		}

		# stack subjects and visits: number of visits as rows, and number of basis as columns
		phi_t_stacked = rbind(phi_t_stacked, t(phi_t[[i]]))
	}

	results_basis = list(knot_place=knots, time_cont=time_cont, orth_spline_basis_sparse=phi_t, 
						 orth_spline_basis_sparse_stacked=phi_t_stacked, orth_spline_basis_cont=phi_t_cont)
	return(results_basis)
}

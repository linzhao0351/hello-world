#@author: Lin Zhao
#created on: Oct 31, 2018

#####
# pseudo code for MCMC estimation
#####

using DataFrames

##### Initialization

### data required

waitlist = Array{Bool, 2}(undef, nStudents, nPrograms)
Distance = Matrix{Float64}(undef, nStudents, nPrograms)


onPlatformPrograms = Dataframe()
# should be like
# ID      Name
#  1	  on program 1
#  2	  on program 2
# ...
#  J	  on program J
onPlatformProgramsDict = Dict{String, Int64}()
# should be like
# on program 1 => 1
# on program 2 => 2
# ...
# on program 3 => J

offPlatformPrograms = Dataframe()
# should be like
# ID      Name
#  1	  on program 1
#  2	  on program 2
# ...
#  K	  on program K
offPlatformProgramsDict = Dict{String, Int64}()
# should be like (Name => NUM)
# off program 1 => 1
# off program 2 => 2
# ...
# off program 3 => J

### constant params

nStudents
nPrograms
nOnPlatform
nOffPlatform

nIteration

lStudentX
lProgramX

### define types

mutable struct student
	ID::String
	X::Array{Float64}
	application::DataFrame(order = Int64[],
						   programID = String[],
						   status = String[])
	finalPlacement::Int64 # should be enrolled program's id
	finalPlacementType::Bool # 0: On; 1: Off
	waitlistedPrograms::DataFrame(order = Int64[],
								  programID = String[],
								  status = String[])
	## Following not sure
	rankInPrograms::Dict{String, Int} # should be (program id, rank)
	initialPlacement::Int64
	application_offPlatform::DataFrame(order = Int64[],
						   programID = String[],
						   status = String[])
end

Students = Array{student, 1}(undef, nStudents)

mutable struct program
	NUM::Int64
	ID::String
	X::Array{Float64}
	onPlatform::Bool
	## Following not sure
	pi::Int64
	piWaitlist::Int64
	capacity::Int64
end

Programs = Array{program, 1}(undef, nPrograms)

function sIndexRule(i::Int64, j::Int64)
	# Given student i and program j, determine sIndex[i, j]

	return sIndex_ij
end

sIndex = zeros(Int64, (nStudents, nOnPlatform)) # This should be determined ex-ante
sIndexOff = zeros(Int64, (nStudents, nOffPlatform))


### data augmentation: These values will be updated in each iteration, but are temporal

utilityOn = zeros(Float64, (nStudents, nOnPlatform))
utilityOff = zeros(Float64, (nStudents, nOffPlatform))

feas = zeros(Bool, (nStudents, nOffPlatform))
accept = zeros(Bool, (nStudents, nOffPlatform))
offer = zeros(Bool, (nStudents, nOffPlatform))

probWaitlistOffer = Array{DataFrame, 1}(undef, nStudents)
# DataFrame here should be
# order		programID 		prob
#   1		  ID1			p_i1
#   2		  ID2			p_i2

### parameters: These values will be updated ub each iteration and will be stored

### On Platform: u_ij = delta_j + lambda * D_ij + eta_ij + epsilon_ij

## (1) delta_j = x_j * beta_bar + ksi_j
# beta_bar ~ MN(0, Sigma_beta_bar)
beta_bar = Array{Array{Float64, 1}, 1}(undef, nIteration)
prior_Sigma_beta_bar = # prior

# ksi_j ~ N(0, sigma_ksi^2), sigma_ksi^2 ~ IGamma(\bar{nu_ksi}, \bar{sigma_ksi}^2)
ksi = zeros(Float64, (nIteration, nPrograms))

sigma2_ksi = zeros(Float64, nIteration)
prior_nu_ksi = # prior
prior_sigma2_ksi = # prior

## (2) lambda ~ N(0, \bar{sigma_lambda}^2)
lambda = zeros(Float64, nIteration)
prior_sigma2_lambda = # prior

## (3) eta_ij = \sum \sum eta_o_lk * studentX_k * programX_l + \sum eta_u_l * programX_l
# eta_o ~ MN(0, Sigma_eta_o)
# eta_u ~ MN(0, Sigma_eta_u), Sigma_eta_u ~ IW(nu_eta_u, Sigma_eta_u)
# element would be coefficient
eta_o = Array{Array{Float64, 1}, 1}(undef, nIteration)
prior_Sigma_eta_o = # prior

eta_u = Matrix{Array{Float64, 1}}(undef, nIteration, nStudents)

Sigma_eta_u = Array{Matrix{Float64}, 1}(undef, nIteration)
prior_nu_eta_u = # prior
prior_Sigma_eta_u =  #prior

## (4) epsilon ~ N(0, sigma2_epsilon), sigma2_epsilon ~ IGamma(nu_epsilon_bar, sigma2_epsilon_bar)
epsilon = zeros(Float64, nIteration)

sigma2_epsilon = zeros(Float64, nIteration)
prior_nu_epsilon = # prior
priro_sigma2_epsilon = # prior

### alpha
alpha = zeros(Float64, nIteration)
a_bar = # prior
b_bar = # prior


### costRenege
# c_renege ~ log Normal(mu_regene, sigma2_renege)
# - mu_renege ~ N(mu_renege_bar, sigma2_mu_renege_bar)
# - sigma2_renege ~ IGamma(nu_sigma2_renege_bar, sigma2_sigma2_renege_bar)
costRenege = zeros(Float64, (nIteration, nStudents))

mu_cRenege = zeros(Float64, (nIteration, nStudents))
prior_mu0_cRenege = # prior
prior_tau2_cRenege = # prior

sigma2_cRenege = zeros(Float64, (nIteration, nStudents))
prior_nu0_cRenege = # prior
prior_sigma2_cRenege = # prior


### shocks
# w ~ N(-pi_k, sigma2_v)
# denote -pi_k as mu_w
# denote sigma2_v as sigma2_w
# mu_w ~ N(0, sigma2_mu_w_bar)
# sigma2_w ~ IGamma(nu_sigam_w_bar, s_sigma_w_bar)

shocks = Array{Matrix{Float64}, 1}(undef, nIteration)
# Matrix should be (nStudents * nPrograms)

mu_w = zeros(Float64, nIteration)
prior_mu_w = # prior
prior_sigma2_w

sigma2_w = zeros(Float64, nIteration)
prior_nu_w = # prior
prior_s2_w = # prior













### Starting values

# All parameters require proper starting values





##### Estimation Steps
nIter = N

### Estimating Waitlist Chances: get p[i, k]_hat

function update_probWaitlistOffer(i::Int64, alpha::Float64)
	# update prob_ik for student i
	# Method 1:
	# application = Students[i].application
	# NOTE: performance can be improved later #
	# for j = 1:nOnPlatform
	#	if isWaitlisted(i, j) == true # student i is waitlisted at program j
	#		probWaitlistOffer[i, j] = prob_receive_waitlist_offer(alpha, sIndex[i, j])
	#	else
	#		probWaitlistOffer[i, j] = 0
	#	end
	# end

	# Method 2: Initialize probWaitlistOffer as zeros
	waitlistedPrograms = Student[i].waitlistedPrograms
	for jj = 1:length(waitlistedPrograms)
		programID = waitlistedPrograms[:programID][jj]
		programNUM = onPlatformProgramsDict[programID]
		probWaitlistOffer[i, programNUM] = prob_receive_waitlist_offer(alpha, sIndex[i, programNUM])
	end
end

function F_piWaitlist(sIndex_ik::Int64)
	# return the probability of sIndex[i, j] > piWaitlist[j], i.e., F_piWaitlist

	return prob
end

function prob_receive_waitlist_offer(alpha::Float64, sIndex_ik::Int64)
	# return p_ik
	p_ik = (1 - alpha) * (1 - F_piWaitlist(sIndex_ik))
	return p_ik
end

function isWaitlisted(i::Int64, j::Int64)
	application = Students[i].application

	waitlisted = false
	order_j = findall(application[:program] .== onPlatformPrograms[j])
	if length(order_j) > 0
		status = application[:status][order_j]
		waitlisted = (application[:status] == "waitlisted")
	end
	return waitlisted
end


### Update Preferences

function setLowerBound(oldbound, newbound)
	if newbound > oldbound
		lowerBound = newbound
	else
		lowerBound = oldbound
	end
	return lowerBound
end

function setUpperBound(oldbound, newbound)
	if newbound < oldbound
		lowerBound = newbound
	else
		lowerBound = oldbound
	end
	return lowerBound
end


# Step 1(a)

function update_onPlatform_utility(i::Int64, j::Int64)
	application = Students[i].application
	programID_j = onPlatformPrograms[j]

	truncL = -Inf
	truncU = Inf

	# If j is on the waitlist (then must be feasible)
	if isWaitlisted(i, j) == true
		# some constraint 1:
		# - constraint from other wailisted program
		# 	- programs on the waitlist but have higher order give larger utility
		# 	- programs on the waitlist but have lower order give less utility
		# - Initial placement must give lower utility
		# - larger than the utility of all feasible programs
		# - off platform constraint

		# constraint from other wailisted program
		waitlistedPrograms = Students[i].waitlistedPrograms
		num_waitlistedPrograms = length(waitlistedPrograms[:order])
		if num_waitlistedPrograms > 1
			order_j = waitlistedPrograms[:order][findall(waitlistedPrograms[:programID] .== programID)]

			for jj = 1:num_waitlistedPrograms

				programID_jj = waitlistedPrograms[:programID][jj]
				programNUM_jj = onPlatformProgramsDict[jj]

				if programID_jj != programID_j
					if waitlistedPrograms[:order][jj] > order_j
						truncU = setUpperBound(truncU, utilityOn[i, programNUM_jj])
					end

					if waitlistedPrograms[:order][jj] < order_j
						truncL = setLowerBound(truncL, utilityOn[i, programNUM_jj])
					end
				end
			end
		end

		# constraint from initial placement
		initialPlacement = Students[i].initialPlacement
		initialPlacementNUM = onPlatformProgramsDict[initialPlacement]
		truncL = setLowerBound(truncL, utilityOn[i, initialPlacementNUM])

		# constraint from all other feasible programs
		for jj = 1:nOnPatform
			programID_jj = onPlatformPrograms[jj]
			if programID_jj != programID_j
				if isFeasible(i, jj) == true
					truncL = setLowerBound(truncL, utilityOn[i, jj])
				end
			end
		end

		# constraint from off platform programs
		# NOTE: seems that deltaV only imposes constraint on u_ik and u_ij
	end


	# If j is the initial placement (then must be feasible)
	if programID_j == Students[i].initialPlacement
		# some constraint 2:
		# - programs on the waitlist give larger utility
		# - all programs that are feasible give less utility
		# - off platform constraint

		# constraint from waitlist programs
		waitlistedPrograms = Students[i].waitlistedPrograms
		num_waitlistedPrograms = length(waitlistedPrograms[:order])
		if num_waitlistedPrograms > 0
			for jj = 1:num_waitlistedPrograms
				programNUM_jj = onPlatformProgramsDict[jj]
				truncU = setUpperBound(truncU, utilityOn[i, programNUM_jj])
			end
		end

		# constraint from all other feasible programs
		for jj = 1:nOnPlatform
			programID_jj = onPlatformPrograms[jj]
			if programID_jj != programID_j
				if isFeasible(i, jj) == true
					truncL = setLowerBound(truncL, utilityOn[i, jj])
				end
			end
		end

		# constraint from off platform programs
		for kk = 1:nOffPlatform
			if feas[i, k] == 1 && accept[i, k] == 0
				deltaVBound = calcDeltaVBound(i, kk, boundfor = "u_ij")
				truncL = setLowerBound(truncL, deltaVBound)
			end

			if accept_ik == 1
				deltaVBound = calcDeltaVBound(i, kk, boundfor = "u_ij")
				truncU = setUpperBound(truncU, deltaVBound)
			end

		end
	end

	# j is other programs (may not be feasible)
	if isWaitlisted(i, j) == false && programID_j != Students[i].initialPlacement
		if isFeasible(i, j) == true
			# some constraint 3
			# - smaller than initial placemant and waitlist
			# - smaller than waitlist feasible and chosen program ???

			# constraint from waitlist programs
			waitlistedPrograms = Students[i].waitlistedPrograms
			num_waitlistedPrograms = length(waitlistedPrograms[:order])
			if num_waitlistedPrograms > 0
				for jj = 1:num_waitlistedPrograms
					programNUM_jj = onPlatformProgramsDict[jj]
					truncU = setUpperBound(truncU, utilityOn[i, programNUM_jj])
				end
			end

			# constraint from initialPlacement
			initialPlacement = Students[i].initialPlacement
			initialPlacementNUM = onPlatformProgramsDict[initialPlacement]
			truncU = setUpperBound(truncU, utilityOn[i, initialPlacementNUM])

			# constraint from selected programs
			for jj = 1:length(application[:order])
				if isWaitlistFeasible(i, jj) == true
					programNUM_jj = onPlatformProgramsDict[application[:program][jj]]
					truncU = setUpperBound(truncU, utilityOn[i, programNUM_jj])
				end
			end
		else
			# Don't care
		end
	end

	# calculate Xb
	# sample u_ij from truncated normal N(Xb, sigma2_epsilon_bar)
	# update utilityOn[i, j]
	uij_mean = calcUtilityMean(i, j)
	new_uij = rand(TruncatedNormal(uijMean, epsilon[nIter], truncL, truncU), 1)[1]
	utilityOn[i, j] = new_uij

end

function isFeasible(i::Int64, j::Int64)
	# check if student i if feasible at program j
	return sIndex[i, j] > Programs[j].pi
end

function isWaitlistFeasible(i::Int64, j::Int64)
	return  sIndex[i, j] > Programs[j].piWaitlist
end

function calcDeltaVBound(i::Int64, k::Int64, boundfor::String)
	# Given utilityOn and utilityOff
	probs = probWaitlistOffer[i][:prob]
	productPI = 1
	for ii = 1:length(probs)
		productPI = productPI * (1 - probs[ii])
	end

	if boundfor == "u_ij"
		deltaVBound = utilityOff[i, k] - ((1 - productPI) / productPI) * costRenege[i]
	end

	if boundfor == "u_ik"
		j = onPlatformProgramsDict[Students[i].initialPlacement]
		deltaVBound = utilityOn[i, j] + ((1 - productPI) / productPI) * costRenege[i]
	end

	if boundfor == "c"
		j = onPlatformProgramsDict[Students[i].initialPlacement]
		deltaVBound = (utilityOff[i, k] - utilityOn[i, j]) * productPI / (1 - productPI)
	end

	return deltaVBound
end

function calcUtilityMean(i::Int64, j::Int64)
	# U = Xb + eps, this function calculates the Xb part
	studentX = Students[i].X
	programX = Programs[j].X

	# t for temporary

	# delta_j
	t_beta_bar = beta_bar[nIter]
	term1 = studentX' * t_beta_bar + ksi[nIter]

	# lambda * D_ij
	t_lambda = lambda[nIter]
	term2 = t_lambda * Distance[i, j]

	# eta_ij
	interactX = reshape(studentX * programX', (lStudentX * lProgramX, 1))
	t_eta_o = eta_o[nIter]
	term3_1 = interactX' * t_eta_o

	t_eta_u = eta_u[nIter, i]
	term3_2 = programX' * t_eta_u

	return term1 + term2 + term3_1 + term3_2
end





# Step 1(b): update w_ik

function update_shocks(i::Int64, k::Int64)
	truncL = -Inf
	truncU = Inf

	if Students[i].finalPlacement == onPlatformPrograms[k]
		truncL = setLowerBound(truncL, -sIndexOff[i, k])
	end

	# if restrictions on utility would be violated when k is feasible:
	# - if accept_ik = 0, then must have delta V_kj < 0, otherwise need w_ik s.t. k is not feasible
	# - if accept_ik = 1, then must have delta V_kj > 0, otherwise need w_ik s.t. k is feasible
	deltaVi_kj = calcDeltaV(i, k)

	if accept[i, k] == 0
		if deltaVi_kj >= 0
			truncU = setUpperBound(truncU, -sIndexOff[i, k])
		end
	end

	if accept[i, k] == 1
		if deltaVi_kj <= 0
			truncL = setLowerBound(truncL, -sIndexOn[i, k])
		end
	end

	# sample w_ik from truncated Normal N(-pi_k, \sigma_v^2)
	new_w_ik = rand(TruncatedNormal(mu_w[nIter], sigma2_w[nIter], truncL, truncU), 1)[1]
	shocks[nIter]][i, k] = new_w_ik

end

function calcDeltaV(i::Int64, k::Int64)
	probs = probWaitlistOffer[i][:prob]
	productPI = 1
	for ii = 1:length(probs)
		productPI = productPI * (1 - probs[ii])
	end

	deltaVi_kj = -(1 - productPI) * costRenege[nIter, i] + productPI * (utilityOff[i, k] - utilityOn[i, j])
	return deltaVi_kj
end


function update_feas(i::Int64, k::Int64)
	feas[i, k] = (sIndex[i, k] + shocks[i, k] > 0)
end




# Step 1(c): update accept

function update_accept(i::Int64, k::Int64)
	j = Students[i].initialPlacement
	deltaVi_kj = calcDeltaV(i, k)
	if feas[i, k] == 1 && deltaVi_kj > 0
		accept[i, k] = 1
	end
end


# Step 1(d): update offer

function update_offer(i::Int64, j::Int64)
	application = Students[i].application
	programID = onPlatformPrograms[j]

	if programID in application[:programID]
		rank_j = findall(application[:programID] .== programID)

		if programID == Students[i].finalPlacement
			offer[i, j] = 1
		end

		if application[:status][rank_j] == "Waitlist" # data type to be altered
			if Students[i].finalPlacementType == "OffPlatform" # data type to be altered
				offer[i, j] = 0
			end

			if Students[i].finalPlacementType == "OnPlatform"
				rankFinalPlacement = findall(application[:programID] .== Students[i].finalPlacement)

				if rankFinalPlacement > rank_j
					offer[i, j] = rand(Binomial(1, alpha[nIter]), 1)[1]
				else
					offer[i, j] = 0
				end
			end
		end
	else
		offer[i, j] = 0
	end
end


# Step 1(e): update u_ik

function update_offPlatform_utility(i::Int64, k::Int64)

	truncL = -Inf
	truncU = Inf

	if Students[i].finalPlacement == offPlatformPrograms[k]
		# constraint from initialPlacement: deltaVi_kj > 0
		deltaVBound = calcDeltaVBound(i, k, boundfor = "u_ik")
		truncL = setLowerBound(truncL, deltaVBound)

		# constraint from other offPlatform programs
		application_offPlatform = Students[i].application_offPlatform
		for kk = 1:length(application_offPlatform[:order])
			if kk != k
				if application_offPlatform[:status] == "admitted"
					# contraint: ui_k > ui_kk
					programID_kk = application_offPlatform[kk]
					programNUM_kk = offPlatformProgramsDict[programID_kk]
					truncL = setLowerBound(truncL, utilityOff[i, programNUM_kk])
				end
			end
		end
	end

	# If Students[i].finalPlacement is off from waitlit
	if Students[i].finalPlacement in Students[i].waitlistedPrograms[:programID]
		# constraint ui_j' - c_i > ui_k
		finalPlacementNUM = onPlatformProgramsDict[Students[i].finalPlacement]
		truncU = setUpperBound(truncU, utilityOn[finalPlacementNUM] - costRenege[nIter, i])
	end

	if # Students[i].finalPlacement is other offPlatform
		# some constraint

	end

	if # Students[i].finalPlacement is initial onPlatform placement
		# some constraint
	end

	# update off platform utility
	uik_mean = calcUtilityMean(i, k) - costSignUp
	new_uik = rand(TruncatedNormal(uik_mean, ksi[nIter], truncL, truncU), 1)[1]
	utilityOff[i, k] = new

end


# Step 1(f): update ci_renege

function update_costRenege(i::Int64)

	truncL = -Inf
	truncU = Inf

	# If Student[i].finalPlacement is off from waitlist
	if Students[i].finalPlacement in Students[i].waitlistedPrograms[:programID]
		# constraint: ui_j' - ci_renege > ui_k
		finalPlacementNUM = onPlatformProgramsDict[Students[i].finalPlacement]
		truncU = setUpperBound(truncU, utilityOn[finalPlacementNUM] - utilityOff[i, k])
	end

	# If Student[i].finalPlacement is some offPlatform program k
	if Student[i].finalPlacementType == 1
		# constraint: deltaVi_kj > 0
		k = offPlatformProgramsDict[Students[i].finalPlacement]
		deltaVBound = calcDeltaVBound(i, k, boundfor = "c")
		truncU = setUpperBound(truncU, deltaVBound)
	end

	# update c_renege_i
	new_c = exp(rand(TruncatedNormal(mu_renege, sigma2_renege, truncL, truncU), 1)[1])
	costRenege[nIter, i] = new_c
end





### Update parameters

function update_costRenege_mean()
	# Normal
	# number of observation used = nStudents
	t_cRenege = costRenege[nIter, :]
	t_cRenege_mean = mean(t_costRenege)

	post_mu_cRenege = (prior_mu0_cRenege / prior_tau2_cRenege + nStudents * t_cRenege_mean / sigma2_cRenege) /
					  (1 / prior_tau2_cRenege + nStudents / sigma2_cRenege)
	post_sigma2 = (1 / prior_tau2_cRenege + nStudents / sigma2_cRenege)^(-1)

	mu_cRenege[nIter] = rand(Normal(post_mu_cRenege, post_sigma2), 1)[1]
end

function update_costRenege_variance()
	# IGamma
	# number of observation used = nStudents
	t_cRenege = costRenege[nIter, :]
	s2 = sum((t_cRenege .- mu_cRenege[nIter]).^2) / n

	post_nu_cRenege = prior_nu0_cRenege + nStudents
	post_sigma2_cRenege = (prior_nu0_cRenege * prior_sigma2_cRenege + nStudents * s2) / post_nu_cRenege

	sigma2_cRenege[nIter] = rand(InverseGamma(post_nu_cRenege / 2, post_nu_cRenege * post_sigma2_cRenege / 2), 1)[1]
end

function update_beta_bar()
	# Normal
	# number of observation used = nPrograms * nStudents

	# construct Data
	Y = zeros(Float64, (nStudents, nPrograms))
	for i = 1:nStudents
		for j = 1:nPrograms
			Y[i, j] = deductUtility(i, j, deduct = "beta_bar")
		end
	end
	Y = reshape(Y, (nStudents * nPrograms), 1)

	X =  zeros(Float64, nStudents * nPrograms)
	for j = 1:nPrograms
		programX = Programs[j].X
		for i = 1:nStudents
			X[((j - 1) * nStudents + i), :] = programX'
		end
	end

	post_Sigma_beta_bar = (prior_Sigma_beta_bar^(-1) + X' * X ./ sigma2_epsilon)^(-1)
	beta_bar[nIter] = rand(MultivariateNormal(zeros(Float64, lProgramX), post_Sigma_beta_bar), 1)[1]
end

function update_ksi(j::Int64)
	# Normal
	# number of observation used = nPrograms
end

function update_sigma2_ksi(j::Int64)
	 # Normal
	 # number of observations used = nPrograms
end

function update_lambda()
	# Normal
	# number of observations used = nPrograms * nStudents
end

function update_random_coefficient(i::Int64)
	# Normal
	# number of observation used = nPrograms (of the same student)
end

function update_random_coefficient_variance()
	# IW
	# number of observation used = nStudents
end

function update_shock_mean()
	# Normal
	# number of observation used = nStudents * nOffPlatform
end

function update_shock_variance()
	# IGamma
	# number of observation used = nStudents * nOffPlatform
end

function update_alpha()
	# Beta
	# number of observation used = nStudents * nOnPatform
	# Is it reasonable ???
end

function deductUtility(i::Int64, j::Int64, deduct::String)
	studentX = Students[i].X
	programX = Programs[j].X

	# delta_j
	t_beta_bar = beta_bar[nIter]
	term1_1 = studentX' * t_beta_bar

	term1_2 = ksi[nIter]

	# lambda * D_ij
	t_lambda = lambda[nIter]
	term2 = t_lambda * Distance[i, j]

	# eta_ij
	interactX = reshape(studentX * programX', (lStudentX * lProgramX, 1))
	t_eta_o = eta_o[nIter]
	term3_1 = interactX' * t_eta_o

	t_eta_u = eta_u[nIter, i]
	term3_2 = programX' * t_eta_u

	utilityMean = term1_1 + term1_2 + term2 + term3_1 + term3_2

	if deduct == "beta_bar"
		return utilityOn[i, j] - (utilityMean - term1_1)
	end

	if deduct == "ksi"
		return utilityOn[i, j] - (utilityMean - term1_2)
	end

	if deduct == "lambda"
		return utilityOn[i, j] - (utilityMean - term2)
	end

	if deduct == "eta_o"
		return utilityOn[i, j] - (utilityMean - term3_1)
	end

	if deduct == "eta_u"
		return utilityOn[i, j] - (utilityMean - term3_2)
	end
end

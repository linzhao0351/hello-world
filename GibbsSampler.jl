#@author: Lin Zhao
#created on: Oct 31, 2018

#####
# pseudo code for MCMC estimation
#####

using DataFrames

##### Initialization

### data required

waitlist = Array{Bool, 2}(undef, nStudents, nPrograms)

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

offPlarformPrograms = Dataframe()
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

nStudentX
nProgramX

### define types

mutable struct student
	ID::String
	X::Array{Float64}
	application::DataFrame(order = Int64[],
						   programID = String[],
						   status = String[])
	finalPlacement::Int64 # should be enrolled program's id
	finalPlacementType::Bool
	waitlistedPrograms::DataFrame(order = Int64[],
								  programID = String[])
	## Following not sure
	rankInPrograms::Dict{String, Int} # should be (program id, rank)
	initialPlacement::Int64
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

sIndex = zeros(Int64, (nStudents, nPrograms)) # This should be determined ex-ante

### data augmentation: These values will be updated in each iteration, but are temporal

utilityOn = zeros(Float64, (nStudent, nOnPlatform))
utilityOff = zeros(Float64, (nStudent, nOffPlatform))

feas = zeros(Bool, (nStudent, nOffPlatform))
accept = zeros(Bool, (nStudent, nOffPlatform))
offer = zeros(Bool, (nStudent, nOffPlatform))

probWaitlistOffer = zeros(Float64, (nStudents, nOnPlatform))


### parameters: These values will be updated ub each iteration and will be stored

### On Platform: u_ij = delta_j + lambda * D_ij + eta_ij + epsilon_ij

## (1) delta_j = x_j * beta_bar + ksi_j
# beta_bar ~ MN(0, Sigma_beta_bar)
beta_bar = Array{Array{Float64, 1}, 1}(undef, nIteration)
Sigma_beta_bar = # prior

# ksi_j ~ N(0, sigma_ksi^2), sigma_ksi^2 ~ IGamma(\bar{nu_ksi}, \bar{sigma_ksi}^2)
ksi = zeros(Float64, (nIteration, nPrograms))

sigma2_ksi = zeros(Float64, nIteration)
nu_ksi_bar = # prior
sigam2_ksi_bar = # prior

## (2) lambda ~ N(0, \bar{sigma_lambda}^2)
lambda = zeros(Float64, nIteration)
sigma2_lambda_bar = # prior

## (3) eta_ij = \sum \sum eta_o_lk * studentX_k * programX_l + \sum eta_u_l * programX_l
# eta_o ~ MN(0, Sigma_eta_o), Sigma_eta_o ~ IW(nu_eta_o_bar, S_eta_o_bar)
eta_o = Matrix{Array{Float64, 1}}(undef, nIteration, nStudent) # element would be coefficient

Sigma_eta_o = Array{Matrix{Float64}, 1}(undef, nIteration)
nu_eta_o_bar = # prior
S_eta_o_bar = # prior

# eta_u ~ MN(0, Sigma_eta_u_bar)
eta_u = Array{Array{Float64, 1}, 1}(undef, nIteration)
Sigma_eta_u_bar = # prior

## (4) epsilon ~ N(0, sigma2_epsilon), sigma2_epsilon ~ IGamma(nu_epsilon_bar, sigma2_epsilon_bar)
epsilon = zeros(Float64, nIteration)
nu_epsilon_bar = # prior
sigma2_epsilon_bar = # prior

### alpha
alpha = zeros(Float64, nIteration)
a_bar = # prior
b_bar = # prior


### costRenege
# c_renege ~ log Normal(mu_regene, sigma2_renege)
# - mu_renege ~ N(mu_renege_bar, sigma2_mu_renege_bar)
# - sigma2_renege ~ IGamma(nu_sigma2_renege_bar, sigma2_sigma2_renege_bar)
costRenege = zeros(Float64, (nIteration, nStudents))

mu_renege = zeros(Float64, (nIteration, nStudents))
mu_renege_bar = # prior
sigma2_mu_renege_bar = # prior

sigma2_renege = zeros(Float64m (nIteration, nStudents))
nu_sigma2_renege_bar = # prior
sigma2_sigma2_renege_bar = # prior


### shocks
# w ~ N(-pi_k, sigma2_v)
# denote -pi_k as mu_w
# denote sigma2_v as sigma2_w
# mu_w ~ N(0, sigma2_mu_w_bar)
# sigma2_w ~ IGamma(nu_sigam_w_bar, s_sigma_w_bar)

Array{Matrix{Float64}, 1}(undef, nIteration)
shocks = Array{Matrix{Float64}, 1}(undef, nIteration)

mu_w = zeros(Float64, nIteration)
sigma2_mu_w_bar = # prior

sigma2_w = zeros(Float64, nIteration)
nu_sigam_w_bar = # prior
s_sigma_w_bar = # prior













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
		status = application[:status]
		if status == "waitlisted"
			waitlisted = true
		else
			waitlisted = false
		end
	else
		waitlisted = false
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

	if isWaitlisted(i, j) == true # j is on the waitlist (then must be feasible)
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

	if # j is the initial placement (then must be feasible)
		# some constraint 2:
		# - programs on the waitlist give larger utility
		# - all programs that are feasible give less utility
		# - off platform constraint
	end

	if # j is other programs (may not be feasible)
		if isFeasible(i, j)
			# some constraint 3
			# - smaller than initial placemant and waitlist
			# - off platform constraint
		else
			# Don't care
		end
	end

	# calculate Xb
	# sample u_ij from truncated normal N(Xb, sigma2_epsilon_bar)
	# update utilityOn[i, j]
end

function isFeasible(i::Int64, j::Int64)
	# check if student i if feasible at program j
	return true
end

function isWaitlistFeasible(i::Int64, j::Int64)
	if sIndex[i, j] > Programs[j].piWaitlist
		return true
	else
		return false
	end
end

function calc_deltaVBound(i::Int64, j::Int64, k::Int64)
	# Given utilityOn and utilityOff



	return deltaVBound
end

# Step 1(b): update w_ik

function update_shocks(i::Int64, k::Int64)
	truncL = -Inf
	truncU = Inf

	if Students[i].finalPlacement == k
		if -sIndex[i, k] > truncL
			truncL = -sIndex[i, k]
		end
	end

	# if restrictions on utility would be violated when k is feasible:
	# - violate for all k such that feas_ik = 0 and accept_ik = 0, delta V_kj < 0
	if accept[i, k] == 0
		for j = 1:nOnPlatform
			if isFeasible(i, j) == true
				deltaVi_kj = calc_deltaV(i, j, k)
				if deltaVi_kj >= 0
					if -sIndex[i, k] < truncU
						truncU = -sIndex[i, k]
					end
				end
			end
		end
	end

	# sample w_ik from truncated Normal N(-pi_k, \sigma_v^2)

end

function calc_deltaV(i::Int64, j::Int64, k::Int64)


	return deltaVi_kj
end


function update_feas(i::Int64, k::Int64)
	if sIndex[i, k] > shocks[i, k]
		feas[i, k] = true
	else
		feas[i, k] = false
	end
end


# Step 1(c): update accept

function update_accept(i::Int64, k::Int64)
	j = Students[i].initialPlacement
	deltaVi_kj = calc_deltaV(i, j, k)
	if feas[i, k] == true & deltaVi_kj > 0
		accept[i, k] = 1
	end
end


# Step 1(d): update offer

function update_offer(i::Int64, j::Int64)
	application = Students[i].application
	if j in application[2,:]
		rank_j = findall(application[2,:] .== j)

		if application[3, rank_j] == Students[i].finalPlacement
			offer[i, j] = 0
		end

		if application[3, rank_j] == "Waitlist" # data type to be altered
			if Students[i].finalPlacementType == "OffPlatform" # data type to be altered
				offer[i, j] = 0
			end

			if Students[i].finalPlacementType == "OnPlatform"
				rankFinalPlacement = findall(application[2,:] . == Students[i].finalPlacement)

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

	if Students[i].finalPlacement == k
		# constraint deltaVi_kj > 0

		for kk = 1:nOffPlatform
			if kk != k
				# contraint: ui_k > ui_kk
			end
		end

	end

	if # Students[i].finalPlacement is off from waitlit
		# constraint ui_j' - c_i > ui_k
	end

	if # Students[i].finalPlacement is other offPlatform
		# some constraint
	end

	if # Students[i].finalPlacement is initial onPlatform placement
		# some constraint
	end

end


# Step 1(f): update ci_renege

function update_costRenege(i::Int64)

	truncL = -Inf
	truncU = Inf

	if # Student[i].finalPlacement is off from waitlist
		# constraint: ui_j' - ci_renege > ui_k
	end

	if # Student[i].finalPlacement is some off Platform program k
		# constraint: deltaVi_kj > 0
	end
end





### Update parameters

function update_costRenege_mean()
	# Normal
	# number of observation used = nStudents
end

function update_costRenege_variance()
	# IGamma
	# number of observation used = nStudents
end

function update_beta_bar()
	# Normal
	# number of observation used = nPrograms * nStudents
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


data <- read.table("gencounts.txt",header=TRUE)

# get number of positions and number of individuals
npos <- dim(data)[1]
nind <- dim(data)[2]/4

# redimension data as a 3-dim array
countmatrix <- as.vector(as.matrix(data))
dim(countmatrix) <- c(npos,4,nind)

# m_i(a)
basecount <- array(1,dim=c(npos,4))

# n_ij^tot
indcount <- array(0,dim=c(npos,nind))

for (i in 1:npos)	{
	for (j in 1:nind)	{
		for (a in 1:4)	{
			basecount[i,a] <- basecount[i,a] + countmatrix[i,a,j]
			indcount[i,j] <- indcount[i,j] + countmatrix[i,a,j]
		}
	}
}

pi <- array(0,dim=c(npos,4))
for (i in 1:npos)	{
	pi[i,] <- basecount[i,] / sum(basecount[i,])	
}


# define array for total log likelihood 
loglarray <- array(0,dim=30)
epsilonarray <- array(0,dim=30)

# a loop over all values of epsilon, from 0.01 to 0.30
for (e in 1:30)	{

	# define current error rate 
	eps <- 0.01 * e

	# define arrays for:

	# best genotype ...
	genotype <- array(0,dim=c(npos,nind))

	# ... and corrresponding posterior probability score
	score <- array(0,dim=c(npos,nind))

	lnL <- 0

	for (i in 1:npos)	{
		for (j in 1:nind)	{

			# store prior and likelihood for the 10 genotypes
			prior <- array(0,dim=4)
			likelihood <- array(0,dim=4)

			# loop over all possible genotypes
			for (a in 1:4)	{

                prior[a] <- pi[i,a];
                likelihood[a] <- (1-eps)^countmatrix[i,a,j] * (eps/3)^(indcount[i,j]-countmatrix[i,a,j]);
			}

			# joint probability is product of prior and likelihood
			# thus, sum in log
			jointprob <- prior * likelihood

			# now make the sum
			marginalprob <- sum(jointprob)

			# and normalize, to get posterior probabilities
			posterior <- jointprob / marginalprob

			# log likelihood
			lnL <- lnL + log(marginalprob)

			# get best genotype
			genotype[i,j] <- which.max(posterior)
			score[i,j] <- max(posterior)
		}
	}

	# store value of lnL and epsilon
	loglarray[e] <- lnL
	epsilonarray[e] <- eps
}

# plot lnL as a function of epsilon
plot(loglarray ~ epsilonarray)



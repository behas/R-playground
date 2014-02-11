#!/usr/bin/env Rscript


#############################################
# R-file for simulating PageRank Algorithm	#
#											#
# (c) 2009 University of Vienna				#
# Author: bernhard.haslhofer@univie.ac.at	#
#############################################

# Note: the algorithm has been implemented following the instructions in
# "Introduction to Information Retrieval" by C. Manning et. al, Chapter 21
# The book is available at: http://www-csli.stanford.edu/~hinrich/information-retrieval-book.html

# Step 1: takes a matrix with dimension N x N and
# replaces all elements in all rows that have no 1's by 1/N
step1 <- function(m) {
	
	if(class(m) == "matrix") {
		
		n <- 1 / dim(m)[1]
		
		for(i in 1:nrow(m)) {
			
			if(sum(m[i,]) == 0) {
				
				if(debug)
					cat("Replacing each 1 in row",i,"by", n ,"\n\n")

				m[i,] <- m[i,] + n
			}
		
		}
		
		return (m)
		
	} else {
		cat("cannot handle class", class(m), "\n\n")
	}
	
}

# Step 2: takes a matrix with dimension N x N and
# divides each 1 in A by the number of 1's in its row
step2 <- function(m) {
	
	if(class(m) == "matrix") {
	
		for(i in 1:nrow(m)) {
		
			row <- m[i,]
			
			sum1 <- sum ( row[is.numeric(row) & row == 1] )
			
			if(sum1 > 0) {

				if(debug)
					cat("Dividing each 1 in row",i,"by",sum1,"\n")
				
				for(j in 1:length(row)) {
					
					if(row[j] == 1) {
						row[j] <- row[j] / sum1
					}
					
				}
				
				m[i,] <- row
			}
			
		
		}
		
		return (m)
		
	} else {
			cat("cannot handle class", class(m), "\n\n")
	}
	
	
}


# Step 1-4: calulates the transition probability matrix from a given adjacency matrix
# and a given teleport factor
calculateTransitionProbabilty <- function(adjacencyMatrix, teleportFactor) {
	
	m1 <- step1(adjacencyMatrix);
	if(debug)
		cat("Computing P - intermediate step 1:",m1,"\n")

	
	m2 <- step2(m1);
	if(debug)
		cat("Computing P - intermediate step 2:",m2,"\n")
	
	m3 <- m2 * (1-teleportFactor);
	if(debug)
		cat("Computing P - intermediate step 3:",m3,"\n")
	
	m4 <- m3 + teleportFactor / dim(adjacencyMatrix)[1];
	
	return (m4);
	
	
}


# Performs a random walk with a given transition probability matrix,
# a given probability distribution vector, and a given number of iterations.
# Records the sequence of proability vectors in a result matrix
randomWalk <- function(probMatrix,probDistrVector,noIterations) {
	
	if(length(probDistrVector != ncol(probMatrix))) {
		cat("The length of the probability vector must be equal to the dimension of the transition probabilty matrix")
	}
	
	# create result matrix
	result <- matrix(nrow=noIterations, ncol=length(probDistrVector))
	
	# copy the columnames
	colnames(result) <- colnames(probMatrix)	

	rownames(result) <- paste(c("x"), 0:(noIterations-1), sep="")
		
	for(i in 1:noIterations) {
				
		if( i == 1) {

			result[i,] <- probDistrVector
			
		} else {
			
			result[i,] <- result[i-1,] %*% probMatrix		
		}
		
		
	}
	
	return (result)
	
}

# plots the result to a graph
drawResult <- function(m) {
	
	# pdf(file="result.pdf", height=3.5, width=3.5)
	
	matplot(result, type="l", lty = 1:ncol(m), col = 1:ncol(m), lwd = 2, xlab = "Steps", ylab = "Rank")

	title(main="Page Rank Result", col.main="red")

	legend("topright", colnames(result), lty = 1:ncol(m), col = 1:ncol(m), lwd = 2)	
	
	box()
	
	# dev.off()
	
}

# set the debug option
debug <- TRUE

# load data
aTable <- read.table("matrix1.txt", header=TRUE)

# convert data from data-frame to numerical matrix
adjacencyMatrix <- data.matrix(aTable)

# print out the loaded Matrix
adjacencyMatrix

# calculate transition probability matrix with a given teleport factor
probabilityMatrix <- calculateTransitionProbabilty(adjacencyMatrix,0.10)

# print out the transition probability matrix
format(probabilityMatrix, digits=2)

# define a probability distribution vector
probDistrVector <- c(1,0,0,0,0,0,0)

# calculate the result with a given number of iterations
result <- randomWalk(probabilityMatrix,probDistrVector,10)

# print out the result
format(result, digits=2)

# draw the result
drawResult(result)
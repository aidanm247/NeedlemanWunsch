#install.packages("stringr") #This is the way to install from CRAN
library(stringr) 

#Function that takes in two strings and gives the score and alignment based on the Needleman-Wunsch algorithm
    needlemanWunsch <- function(sequence1 , sequence2, m, mm, gp){  
        
        match <- m
        mismatch <- mm
        gapPenalty <- gp

        #create a new matrix set to all zeros
        mat <- matrix(0, nrow = str_length(sequence2)  + 1 , ncol = str_length(sequence1)+ 1)
        mat <- initializeMatrix(mat, sequence1, sequence2, match , mismatch, gapPenalty)

        #fill in the matrix
        for(i in 2:(str_length(sequence2) + 1)){
            for(j in 2:(str_length(sequence1) + 1)){
               mat[i,j] <- findCellValue(sequence1, sequence2, mat, i, j, match, mismatch, gapPenalty)
            }
        }
        #output
        print(mat)
        traceScore(mat, sequence1, sequence2 , m, mm, gp)
    }

#Function to initalize the matrix
    initializeMatrix <- function(mat, sequence1, sequence2, m, mm, gp){
        #set the first row based on the gap penalty
        for(i in 1:str_length(sequence1) + 1 ){
            mat[1,i] <- (i-1) * gp
        }
         #set the first column based on the gap penalty
        for(i in 1:str_length(sequence2) + 1 ){
            mat[i,1] <- (i-1) * gp
        }
        return(mat)
    }

#Function to find the value of a given cell in the matrix
    findCellValue <- function(seq1, seq2,  matrix,  row, col, m, mm , gp){
        match <- m
        mismatch <- mm
        gapPenalty <- gp

        #get the values of the cells to the left, above, and above-left
        left <- matrix[row, col-1]
        above <- matrix[row-1, col]
        aboveLeft <- matrix[row-1, col-1]

        #get the characters of the sequences
        n1 <- substr(seq1, col-1, col-1)
        n2 <- substr(seq2, row-1, row-1)

        #calculate the score for a match, mismatch, and gap
        if(n1 == n2){
            matchScore <- aboveLeft + match
            
        }else{
            matchScore <- aboveLeft + mismatch
        }
        gapScore1 <- above + gapPenalty
        gapScore2 <- left + gapPenalty

        #return the max of the three scores
        return(max(matchScore, gapScore1, gapScore2))
    }

    #Function to traceback through the matrix and return the score
    traceScore <- function(matrix, sequence1, sequence2, m, mm, gp){
        currRow <- str_length(sequence2) + 1
        currCol <- str_length(sequence1) + 1 
        currentCell <- matrix[currRow,currCol]
        match <- m
        mismatch <- mm
        gapPenalty <- gp

        topAlignment <- ""
        bottomAlignment <- ""

        #find where in the matrix the current cell got it's value from
        while(currRow > 0 && currCol > 0){
            left <- matrix[currRow, currCol-1]
            above <- matrix[currRow-1, currCol]
            aboveLeft <- matrix[currRow-1, currCol-1]

            n1 <- substr(sequence1, currCol-1, currCol-1)
            n2 <- substr(sequence2, currRow-1, currRow-1)

            if(n1 == n2 ){
                matchScore <- aboveLeft + match
            }else{
                matchScore <- aboveLeft + mismatch
            }
            aboveGap <- above + gapPenalty
            gapScore2 <- left + gapPenalty

            #add guards to prevent out of bounds errors
            if(currRow != 1 && currCol != 1 && currentCell == matchScore){
                topAlignment <- paste(n2, topAlignment, sep = "")
                bottomAlignment <- paste(n1, bottomAlignment, sep = "")
                currRow <- currRow - 1
                currCol <- currCol - 1
            }else if(currRow != 1 && currentCell == aboveGap){
                currRow <- currRow - 1
                topAlignment <- paste(n2, topAlignment, sep = "")
                bottomAlignment <- paste("-", bottomAlignment, sep = "")
            }else if(currCol != 1 && currentCell == gapScore2){
                currCol <- currCol - 1
                topAlignment <- paste("-", topAlignment, sep = "")
                bottomAlignment <- paste(n1, bottomAlignment,   sep = "")
            }else{
                break
            }
            currentCell <- matrix[currRow, currCol]
        }
       
        score <- 0
        middle <- ""

        for(i in 1:str_length(topAlignment)){
            if(substr(topAlignment, i, i) == substr(bottomAlignment, i, i)){
                score <- score + match
                middle <- paste(middle, "|", sep = "")
            }else if(substr(topAlignment, i, i) == "-" || substr(bottomAlignment, i, i) == "-"){
                score <- score + gapPenalty
                middle <- paste(middle, " ", sep = "")
            }else{  
                score <- score + mismatch
                middle <- paste(middle, ":", sep = "")
            }
            
        }
        print(bottomAlignment)
        print(middle)
        print(topAlignment)
        print(paste("Score: ", score))
     }


#Use this line to test the code
#needlemanWunsch("seq1", "seq2", match, mismatch, gap)
#example: needlemanWunsch("GCGTAC", "GATACG", 3, -1, -2)

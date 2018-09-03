read_seq2geno.tab<-function(f, na.value= NA){
    mat<- <- data.matrix(read.table(f, 
        sep= '\t', header= T, check.names= F, 
        stringsAsFactors= F, row.names= 1))
    if (any(is.na(mat))){
	quit()
    }else{
      return(mat)
    }
}


BEGIN {
    FPAT = "([^,]+)|(\"[^\"]+\")"
}

{
    len = length($4)
    i = substr($4, 2, len - 2) 
    if (i + 0 > threshold){
        len = length($1)
        print substr($1, 2, len - 2)
    } 

}

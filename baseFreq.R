nline <- function(fastq) {
    nline <- 0
    while(length(line <- scan(fastq,skip=nline,n=1,what='character',quiet=T)) != 0)
        nline <- nline + 1
    nline
}
baseFreq0 <- function(fastq, n = 10,start =1, nread = -1,plot=F) {
    if(nread == -1)
        nread = nline(fastq) / 4; ##total number of reads in this file
    skip <- 1
     freq <- matrix(0, ncol = n, nrow = 5)
    rownames(freq) <- c('A','C','G','T','N')
    colnames(freq) <- as.character(1:n)
    read.count <- 0

    while(length(line <- scan(fastq,skip=skip,n=1,what='character',quiet=T)) != 0) {

        sequence <- strsplit(line,'')[[1]]

        for(i in 1:(start + n -1) ) {
            switch(sequence[i],
                   A = {freq[1,i] <- freq[1,i] + 1},
                   C = {freq[2,i] <- freq[2,i] + 1},
                   G = {freq[3,i] <- freq[3,i] + 1},
                   T = {freq[4,i] <- freq[4,i] + 1},
                   N = {freq[5,i] <- freq[5,i] + 1},
                   stop(paste0('non ATCGN characters found: ',sequence[i],'line: ',paste(sequence,sep='',collaspe=''),'nline: ',nline,'\n')))
        }
        read.count <- read.count + 1
        if(read.count == nread)
            break
        skip <- skip + 4
    }
    freq <- freq / nread
    if(plot) {
        if(!require(seqLogo))
            stop('try to install seqLogo package for sequence logo plot\n')
        ##we have to remove N before feed makePWM with freq;assign weight of N equally to ATCG
        pwm.mtr <- freq[1:4,]
        for(i in 1:ncol(pwm.mtr) ) {
            if((N.w <- freq[5,i]) != 0 )
                pwm.mtr[,i] <- pwm.mtr[,i] + (N.w/4)
        }
        plot(makePWM(pwm.mtr)) #N is not supported by seqLogo?
    }
    round(freq,2)
}

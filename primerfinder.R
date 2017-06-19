#' Find Primers for RNA staining
findPrimers <- function(in.file, tm.min=45, tm.max=65, ending=c("TA","AG"), uid=NULL) {
  
  
  if (is.null(uid)) uid <- paste("uid", gsub("[^[:digit:]]", "", as.character(as.numeric(Sys.time()))), sep ="")

  output.dir <- "output"
  
  uid.dir <- paste(getwd(), "/", output.dir, "/", uid, "/", sep="")
  if (dir.exists(uid.dir)) {
    
  } else {
    dir.create(uid.dir)
  }
  
  cat("\nStart FindPrimers analysis\n")
  cat("\nAll results will be saved to the folder", uid.dir, "\n")
  
  ##################################
  # read input sequence from file #
  ##################################
  
  cat("\nReading RNA bases from", in.file, "\n")
  input <- read.table(file=paste(in.file, sep=""), sep="\t")
  seq.in <- ""
  for (i in 1:dim(input)[1]) {
    seq.i <- as.character(gsub("[^[:alpha:]]", "", as.character(input[i,1])))
    seq.in <- c(seq.in, seq.i)
  }
  seq.in <- paste(seq.in, collapse = "")
  seq.in <- toupper(as.character(seq.in))
  
  total_bp <- nchar(seq.in)
  
  
  if (total_bp<=40) {
    stop("Less than 40 bases detected.")
  } else {
    cat("Loaded", total_bp, "bases\n")
  }
  
  if (length(grep("[^(A|G|C|T|U)]", seq.in))>0) {
    stop("The input contains invalid bases. Only A, T, C, G, U are accepted at this stage.\n")
  }
  
  
  ##################################
  # generate input file to primer3 #
  ##################################
  
  
  primer3.in.filename <- paste(uid.dir, "primer3.in.txt", sep="")
  cat("\nLoading primer3 settings from primer3.in.template")
  primer3.in <- read.table(file="primer3.in.template", sep=",", stringsAsFactors = FALSE)
  primer3.in[1,1] <- paste("SEQUENCE_ID=", uid, sep="")
  primer3.in[2,1] <- paste("SEQUENCE_TEMPLATE=", seq.in, sep="")
  write.table(primer3.in, file=primer3.in.filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  ##################################
  #########   run primer3   ########
  ##################################
  
  cat("\nCalling primer3 to detect candidate primers\n")
  #print(primer3.in[3:(dim(primer3.in)[1]-1), 1])

  wd <- getwd()
  setwd(paste(wd, "/primer3", sep=""))
  primer3.command <- paste(getwd(), "/primer3_core.exe ", primer3.in.filename, sep="")
  system(primer3.command, ignore.stdout = TRUE, show.output.on.console = FALSE)
  primer3.out <- read.table(file=paste(uid, ".int", sep=""), skip=2, comment.char = "", head=T, row.names=1)
  # The index provided by primer3 should add 1
  primer3.out$start <- primer3.out$start+1
  rownames(primer3.out) <- NULL
  file.remove(file=paste(uid, ".int", sep=""))
  setwd(paste(wd, sep=""))
  
  if (dim(primer3.out)[1]<=0) {
    cat("\nNo candidates detected\n")
    return(1)
  }
  
  write.table(data.frame(PID=rownames(primer3.out), primer3.out, check.names=F), file=paste(uid.dir,"primers.1.all.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=FALSE)
  cat("\nDetected", dim(primer3.out)[1], "candidates\n")

  ##########################################
  ##### pick primers with valid ending #####
  ##########################################
  
  primer3.out$valid_end <- 0
  for (i in 1:dim(primer3.out)[1]) {
    i.end <- substr(as.character(primer3.out$sequence[i]), 20, 21)
    if (i.end %in% ending) primer3.out$valid_end[i] <- 1
  }
  primer3.out <- primer3.out[which(primer3.out$valid_end==1), ]
  
  
  if (dim(primer3.out)[1]<=0) {
    cat("\nNo candidates with valid ending {", paste(ending, collapse =","), "}\n")
    return(1)
  }
  
  write.table(data.frame(PID=rownames(primer3.out), primer3.out, check.names=F), file=paste(uid.dir,"primers.2.validending.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=FALSE)
  cat("\nKeep", dim(primer3.out)[1], "candidates ending with {", paste(ending, collapse =","), "}\n")
  
  #######################################################################################
  ##### pick primer pairs and validate melting temperature and self-complementarity #####
  #######################################################################################
  
  cat("\nPicking primer pairs and checking melting temperature and self-complementarity\n")
  primers <- data.frame(primer3.out[, c("sequence","start")], stringsAsFactors = FALSE)
  colnames(primers) <- c("left.seq","left.start")
  
  primers <- primers[which((primers$left.start+39)<=total_bp), ]
  
  if (dim(primers)[1] <=0) {
    cat("\nNo candidates start from at least 40 bases from the end of the input sequence\n")
    return(1)
  }
  
  
  primers$right.seq <- NA
  primers$left.seq.rc <- NA
  primers$right.seq.rc <-NA
  primers$left.bTM <- 0
  primers$left.aTM <- 0
  primers$left.nnTM <- 0
  primers$right.bTM <- 0
  primers$right.aTM <- 0
  primers$right.nnTM <- 0
  primers$left.hairpin <- 1
  primers$left.dimers <- 1
  primers$right.hairpin <- 1
  primers$right.dimers <- 1
  primers$valid <- 0
  primers$left.seq <- NA
  
  for (i in 1:dim(primers)[1]) {
    
  #  if ((nchar(seq.in)-primers$left.start[i])>=20) {
      primers$left.seq[i] <- substr(seq.in, primers$left.start[i], primers$left.start[i]+19)
      primers$right.seq[i] <- substr(seq.in, primers$left.start[i]+20, primers$left.start[i]+39)
      primers$left.seq.rc[i] <- as.character(getReverseComplement(as.character(primers$left.seq[i])))
      primers$right.seq.rc[i] <- as.character(getReverseComplement(as.character(primers$right.seq[i])))
      
      primers$left.bTM[i] <- getBasicTemperature(primers$left.seq.rc[i])$min
      primers$left.aTM[i] <- getAdjustedTemperature(primers$left.seq.rc[i])$min
      primers$left.nnTM[i] <- getNNTemperature(primers$left.seq.rc[i])$min
      primers$left.hairpin[i] <- checkHairpin(primers$left.seq.rc[i])
      primers$left.dimers[i] <- checkDimers(primers$left.seq.rc[i])
      
      primers$right.bTM[i] <- getBasicTemperature(primers$right.seq.rc[i])$min
      primers$right.aTM[i] <- getAdjustedTemperature(primers$right.seq.rc[i])$min
      primers$right.nnTM[i] <- getNNTemperature(primers$right.seq.rc[i])$min
      primers$right.hairpin[i] <- checkHairpin(primers$right.seq.rc[i])
      primers$right.dimers[i] <- checkDimers(primers$right.seq.rc[i])
  
      i.tms <- as.numeric(primers[i, c("left.bTM","left.aTM","left.nnTM","right.bTM","right.aTM", "right.nnTM")])
      i.sc <- sum(as.numeric(primers[i, c("left.hairpin","left.dimers","right.hairpin","right.dimers")]))
      if (all(i.tms<=tm.max) & all(i.tms>=tm.min) & i.sc==0) primers$valid[i] <- 1
  #  }
  }
  
  write.table(primers, file=paste(uid.dir, "primers.3.pairs.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=FALSE)
  
  primers <- primers[which(primers$valid==1), ]
  
  cat("\nDetected", dim(primers)[1], "valid primer pairs\n")
  
  if (dim(primers)[1]>0) {
    write.table(primers[which(primers$valid==1), ], file=paste(uid.dir, "primers.4.final.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=FALSE)
  }
  
  cat("\nDone!\n")
}

#' Check whether there are any other IUPAC bases in the sequence
#' 
#' @param seq The input sequence
#' @return If the sequence contains IUPAC bases, return 1; otherwise, return 0
#' 
anyIUPACBases <- function(seq) {
  if (length(grep("(M|R|W|S|Y|K|V|H|D|B|N)", seq))>0) {
    return(1)
  }
  return(0)
}

#' check the validity of the input sequence
#' 
#' @param seq The input sequence
#' @return If the sequence contains non-valid bases, return 0; otherwise, return 1
#' 
checkBases <- function(seq) {
  if (length(grep("(A|G|C|U|T|M|R|W|S|Y|K|V|H|D|B|N)", seq))>0) {
    return(0)
  }
  return(1)
}

#' Get Reverse Complement of Nucleotide sequences
#' 
#' @param seq The input Nucleotide sequence
#' @param isDAN TRUE if the sequence is a DNA sequence; otherwise, FALSE.
#' @param do.complement If TRUE, get the complement of the input sequence.
#' @param do.reverse If TRUE, reverse the sequence.
#' @return The converted sequence
#'  
getReverseComplement <- function(seq, isDNA=TRUE, do.complement=TRUE, do.reverse=TRUE) {
  
  ret <- c()
  for (i in 1:nchar(seq)) {
    i.char <- substr(seq, i, i)
    i.ret <- ""
    if (i.char=="A") {
      if (isDNA==TRUE) {i.ret <- "T"}
      else {i.ret <- "U"}
    } else if (i.char=="T" | i.char=="U") {
      i.ret <- "A"
    } else if (i.char=="G") {
      i.ret <- "C"
    } else if (i.char=="C") {
      i.ret <- "G"
    } else if (i.char=="M") {
      i.ret <- "K"
    } else if (i.char=="K") {
      i.ret <- "M"
    } else if (i.char=="R") {
      i.ret <- "Y"
    } else if (i.char=="Y") {
      i.ret <- "R"
    } else if (i.char=="W") {
      i.ret <- "W"
    } else if (i.char=="S") {
      i.ret <- "S"
    } else if (i.char=="V") {
      i.ret <- "B"
    } else if (i.char=="B") {
      i.ret <- "V"
    } else if (i.char=="H") {
      i.ret <- "D"
    } else if (i.char=="D") {
      i.ret <- "H"
    }
    
    if (do.complement) {
      ret <- c(ret, i.ret)
    } else {
      ret <- c(ret, i.char)
    }
  }
  
  if (do.reverse) ret <- rev(ret)
  ret <- paste(ret, collapse="")
  
  return(ret)
}


countOligo <- function(seq) {
  
  seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- length(seq.vector)
  
  aCount <- length(which(seq.vector == "A"))
  gCount <- length(which(seq.vector == "G"))
  cCount <- length(which(seq.vector == "C"))
  tCount <- length(which(seq.vector == "T"))
  sCount <- length(which(seq.vector == "S"))
  uCount <- length(which(seq.vector == "U"))
  wCount <- length(which(seq.vector == "W"))
  iCount <- length(which(seq.vector == "I"))
  rCount <- length(which(seq.vector == "R"))
  yCount <- length(which(seq.vector == "Y"))
  vCount <- length(which(seq.vector == "V"))
  hCount <- length(which(seq.vector == "H"))
  dCount <- length(which(seq.vector == "D"))
  bCount <- length(which(seq.vector == "B"))
  nCount <- length(which(seq.vector == "N"))
  
  aaCount <- countNeighbors("AA", seq) + countNeighbors("TT", seq) + countNeighbors("UU", seq)
  atCount <- countNeighbors("AT", seq) + countNeighbors("AU", seq)
  taCount <- countNeighbors("TA", seq) + countNeighbors("UA", seq)
  caCount <- countNeighbors("CA", seq) + countNeighbors("TG", seq) + countNeighbors("UG", seq)
  gtCount <- countNeighbors("GT", seq) + countNeighbors("AC", seq) + countNeighbors("GU", seq)
  ctCount <- countNeighbors("CT", seq) + countNeighbors("AG", seq) + countNeighbors("CU", seq)
  gaCount <- countNeighbors("GA", seq) + countNeighbors("TC", seq) + countNeighbors("UC", seq)
  cgCount <- countNeighbors("CG", seq) 
  gcCount <- countNeighbors("GC", seq) 
  ggCount <- countNeighbors("GG", seq) + countNeighbors("CC", seq) 
  
  IUpairVals_min <- IUpairVals_max <- c(0,0,0)
  for (i in 2:seq.n) {
    base0 <- substr(seq, i-1, i-1)
    base <- substr(seq, i, i)
    temp <- c(0,0,0)
    temp=calcIUpair(base0, base, i, seq, "min")
    for (j in 1:3) {
      IUpairVals_min[j] <- IUpairVals_min[j] + temp[j]
    }
    temp=calcIUpair(base0, base, i, seq, "max")
    for (j in 1:3) {
      IUpairVals_max[j] <- IUpairVals_max[j] + temp[j]
    }
  }
  
  return(list(aCount=aCount, gCount=gCount, cCount=cCount, tCount=tCount,
              sCount=sCount, uCount=uCount, wCount=wCount, iCount=iCount,
              rCount=rCount, yCount=yCount, vCount=vCount, hCount=hCount,
              dCount=dCount, bCount=bCount, nCount=nCount, aaCount=aaCount,
              atCount=atCount, taCount=taCount, caCount=caCount, gtCount=gtCount,
              ctCount=ctCount, gaCount=gaCount, cgCount=cgCount, gcCount=gcCount,
              ggCount=ggCount, IUpairVals_min=IUpairVals_min, IUpairVals_max=IUpairVals_max))
}

calcIUpair <- function(base0, base, i, seq, choice) {
  reValue <- c(0, 0, 0)
  IUpacBase <- ""
  pair1 <- ""
  pair2 <- ""
  temp1 <- c(0,0,0)
  temp2 <- c(0,0,0)
  base2 <- substr(seq, i+1, i+1)
  
  if (isIUpackBase(base0)==1) return(ret)
  
  if (isIUpackBase(base)==1) {
    if(base=="M"){IUpacBase="AC"}
    else if(base=="R"){IUpacBase="AG"}
    else if(base=="W"){IUpacBase="AT"}
    else if(base=="S"){IUpacBase="CG"}
    else if(base=="Y"){IUpacBase="CT"}
    else if(base=="K"){IUpacBase="GT"}
    else if(base=="V"){IUpacBase="ACG"}
    else if(base=="H"){IUpacBase="ACT"}
    else if(base=="D"){IUpacBase="AGT"}
    else if(base=="B"){IUpacBase="CGT"}
    else if(base=="N"){IUpacBase="ACGT"}
    
    j<-1
    while(substr(IUpacBase, j, j) != "") {
      base <- substr(IUpacBase, j, j) 
      pair1 <- paste(base0, base, sep="")
      
      if(pair1=="AA") {temp1[1]= 1.2; temp1[2]=8.0; temp1[3]=21.9 ;}
      else if(pair1=="AT"){temp1[1]= 0.9; temp1[2]=5.6; temp1[3]=15.2;}
      else if(pair1=="TA"){temp1[1]=0.9; temp1[2]=6.6; temp1[3]= 18.4;}
      else if(pair1=="CA"){temp1[1]=1.7; temp1[2]=8.2; temp1[3]=21.0;}
      else if(pair1=="GT"){temp1[1]= 1.5; temp1[2]=9.4; temp1[3]=25.5;}
      else if(pair1=="CT"){temp1[1]= 1.5; temp1[2]=6.6; temp1[3]=16.4;}
      else if(pair1=="GA"){temp1[1]=1.5; temp1[2]=8.8; temp1[3]=23.5;}
      else if(pair1=="CG"){temp1[1]= 2.8; temp1[2]=11.8; temp1[3]=29.0;}
      else if(pair1=="GC"){temp1[1]=2.3 ; temp1[2]=10.5; temp1[3]=26.4;}
      else if(pair1=="GG"){temp1[1]=2.1 ; temp1[2]=10.9; temp1[3]=28.4;}
      
      if(base2==""){
        for(k in 1:2)
        {	temp2[k]=0.0;	}
        
      }else if(isIUpackBase(base2)==0){
        pair2=paste(base, base2, sep="")
        if(pair2=="AA"){temp2[1]= 1.2; temp2[2]=8.0; temp2[3]=21.9;}
        else if(pair2=="AT"){temp2[1]= 0.9; temp2[2]=5.6; temp2[3]=15.2;}
        else if(pair2=="TA"){temp2[1]=0.9; temp2[2]=6.6; temp2[3]= 18.4;}
        else if(pair2=="CA"){temp2[1]=1.7; temp2[2]=8.2; temp2[3]=21.0;}
        else if(pair2=="GT"){temp2[1]= 1.5; temp2[2]=9.4; temp2[3]=25.5;}
        else if(pair2=="CT"){temp2[1]= 1.5; temp2[2]=6.6; temp2[3]=16.4;}
        else if(pair2=="GA"){temp2[1]=1.5; temp2[2]=8.8; temp2[3]=23.5;}
        else if(pair2=="CG"){temp2[1]= 2.8; temp2[2]=11.8; temp2[3]=29.0;}
        else if(pair2=="GC"){temp2[1]=2.3; temp2[2]=10.5; temp2[3]=26.4;}
        else if(pair2=="GG"){temp2[1]=2.1; temp2[2]=10.9; temp2[3]=28.4;}
      }else if(isIUpackBase(base2)==1){
        base0=base; base=base2; i<- i+1; 
        temp2=CalcIUpair(base0,base,i,seq,choice);
        i <- i-1;
      }
      
      for(k in 1:3) {
        if(j==0){
          reValue[k]=temp1[k]+temp2[k];
        }else{
          if ((choice=="max")&&(reValue[k]<temp1[k]+temp2[k]))
          {	reValue[k]=temp1[k]+temp2[k];	
          }else if((choice=="min")&&(reValue[k]>temp1[k]+temp2[k]))
          {	reValue[k]=temp1[k]+temp2[k]; 	
          }
        }
      }
      j <- j + 1;
    }
    
  }
  return(reValue)
}

isIUpackBase <- function(base) {
  if (base %in% c("M","R","W","S","Y","K","V","H","D","B","N")) {
    return(1)
  } else {
    return(0)
  }
}

countNeighbors <- function(pattern, seq) {
  seq.n <- nchar(seq)
  ret <- 0
  i <- 0
  while(i>=0 & i<seq.n) {
    seq <- substr(seq, i, seq.n)
    i=regexpr(pattern, seq, fixed=TRUE)[[1]]
    if (i>=0) {
      ret <- ret + 1
      i <- i + 1
    }
  }
  return(ret)
}

#' Calculate the basic melting temperature of a sequence
#' 
#' @param seq The input sequence
#' @return The min and max of basic melting temperature
#' 
getBasicTemperature <- function(seq) {
  
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  
  eGC_min <- counts$gCount+counts$cCount+counts$sCount
  eGC_max <- seq.n-counts$aCount-counts$tCount-counts$wCount-counts$uCount
  
  bTM.min <- bTM.max <- 0
  
  if (seq.n>0) {
    if (seq.n<14) {
      bTM.min <- round(2*(seq.n-eGC_min)+4*eGC_min, 1)
      bTM.max <- round(2*(seq.n-eGC_max)+4*eGC_max, 1)
    } else {
      bTM.min <- round(64.9+41*((eGC_min-16.4)/seq.n), 1)
      bTM.max <- round(64.9+41*((eGC_max-16.4)/seq.n), 1)
    }
  }
  return(list(min=bTM.min, max=bTM.max))
}

getGC <- function(seq) {
  gcValmin <- gcValmax <- 0
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  eGC_min <- counts$gCount+counts$cCount+counts$sCount
  eGC_max <- seq.n-counts$aCount-counts$tCount-counts$wCount-counts$uCount
  gcValmin <- round(100*eGC_min/seq.n,1)
  gcValmax <- round(100*eGC_max/seq.n, 1)
  return(list(min=gcValmin, max=gcValmax))
}

#' Calculate the salt adjusted melting temperature of a sequence
#' 
#' @param seq The input sequence
#' @return The min and max of adjusted melting temperature
#' 
getAdjustedTemperature <- function(seq, isDNA=T, saltConcentration=50) {
  
  ret <- 0
  #seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  
  gc <- getGC(seq)
  
  eGC_min <- counts$gCount+counts$cCount+counts$sCount
  eGC_max <- seq.n-counts$aCount-counts$tCount-counts$wCount-counts$uCount
  
  aTM.min <- aTM.max <- 0

  if (seq.n>0) {
    if (isDNA) {
      if (seq.n<14) {
        aTM.min <- round(2*(seq.n-eGC_min)+4*eGC_min+21.6+(7.21*log(saltConcentration/1000)), 1)
        aTM.max <- round(2*(seq.n-eGC_max)+4*eGC_max+21.6+(7.21*log(saltConcentration/1000)), 1)
      } else {
        aTM.min <- round(100.5+(0.41*gc$min)-(820/seq.n)+(7.21*log(saltConcentration/1000)), 1)
        aTM.max <- round(100.5+(0.41*gc$max)-(820/seq.n)+(7.21*log(saltConcentration/1000)), 1)
      }
    } else {
      if (seq.n>20) {
        aTM.min <- round(79.8+(0.584*gc$min)-(11.8*(gc$min/100)*(gc$min/100))-(820/seq.n)+(8.03*log(saltConcentration/1000)), 1)
        aTM.max <- round(79.8+(0.584*gc$max)-(820/seq.n)+(8.03*log(saltConcentration/1000)), 1)
      } else {
        warning("aTM: The input sequence is too short")
      }
    }
  }
  return(list(min=aTM.min, max=aTM.max))
}


DeltaG <- function(seq) {
  dg.min <- dg.max <- 0
  ret <- 0
  #seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  
  if (seq.n>7) {
    val <- -5.0
    val <- val + 1.2*counts$aaCount + 0.9*counts$atCount + 0.9*counts$taCount
    val <- val + 1.7*counts$caCount + 1.5*counts$gtCount + 1.5*counts$ctCount
    val <- val + 1.5*counts$gaCount + 2.8*counts$cgCount + 2.3*counts$gcCount
    val <- val + 2.1*counts$ggCount
    
    dg.min <- val + counts$IUpairVals_min[1]
    dg.max <- val + counts$IUpairVals_max[1]
    
    dg.min <- round(1000*dg.min, 1)/1000
    dg.max <- round(1000*dg.max, 1)/1000
  }
  return(list(min=dg.min, max=dg.max))
}

DeltaH <- function(seq, isDNA=T) {
  dh.min <- dh.max <- 0
  #seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  
  if (seq.n>7) {
    val <- 0
    if (isDNA) {
      val <- val + 8.0*counts$aaCount + 5.6*counts$atCount + 6.6*counts$taCount
      val <- val + 8.2*counts$caCount + 9.4*counts$gtCount + 6.6*counts$ctCount
      val <- val + 8.8*counts$gaCount + 11.8*counts$cgCount + 10.5*counts$gcCount
      val <- val + 10.9*counts$ggCount
    } else {
      val <- val + 6.8*counts$aaCount + 9.38*counts$atCount + 7.69*counts$taCount
      val <- val + 10.44*counts$caCount + 11.4*counts$gtCount + 10.48*counts$ctCount
      val <- val + 12.44*counts$gaCount + 10.64*counts$cgCount + 14.88*counts$gcCount
      val <- val + 13.39*counts$ggCount
    }
    dh.min <- val + counts$IUpairVals_min[2]
    dh.max <- val + counts$IUpairVals_max[2]
    
    dh.min <- round(1000*dh.min, 1)/1000
    dh.max <- round(1000*dh.max, 1)/1000
  }
  return(list(min=dh.min, max=dh.max))
}

DeltaS <- function(seq, isDNA=T) {
  ds.min <- ds.max <- 0
  #seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- nchar(seq)
  counts <- countOligo(seq)
  
  if (seq.n>7) {
    val <- 0
    if (isDNA) {
      val <- val + 21.9*counts$aaCount + 15.2*counts$atCount + 18.4*counts$taCount
      val <- val + 21.0*counts$caCount + 25.5*counts$gtCount + 16.4*counts$ctCount
      val <- val + 23.5*counts$gaCount + 29.0*counts$cgCount + 26.4*counts$gcCount
      val <- val + 28.4*counts$ggCount
    } else {
      val <- val + 19.0*counts$aaCount + 26.7*counts$atCount + 20.5*counts$taCount
      val <- val + 26.9*counts$caCount + 29.5*counts$gtCount + 27.1*counts$ctCount
      val <- val + 32.5*counts$gaCount + 26.7*counts$cgCount + 36.9*counts$gcCount
      val <- val + 32.7*counts$ggCount
    }
    ds.min <- val + counts$IUpairVals_min[3]
    ds.max <- val + counts$IUpairVals_max[3]
    
    ds.min <- round(1000*ds.min, 1)/1000
    ds.max <- round(1000*ds.max, 1)/1000
  }
  return(list(min=ds.min, max=ds.max))
}

#' Calculate the nearest neighbor melting temperature of a sequence
#' 
#' @param seq The input sequence
#' @return The min and max of nn melting temperature
#' 
getNNTemperature <- function(seq, isDNA=T, primerConcentration=50, saltConcentration=50) {
  
  nnTM.min <- nnTM.max <- 0
  #seq.vector <- unlist(strsplit(seq, split=""))
  seq.n <- nchar(seq)
  
  if (seq.n>7) {
    K <- 1/(primerConcentration*1e-9)
    R <- 1.987
    RlnK <- R*log(K)
    RlogK <- round(1000*RlnK,1)/1000
    
    nnTM.min <- 1000*((DeltaH(seq, isDNA=isDNA)$min-3.4)/(DeltaS(seq, isDNA=isDNA)$min+RlnK))
    nnTM.min <- nnTM.min-272.9+7.21*log(saltConcentration/1000)
    
    nnTM.max <- 1000*((DeltaH(seq, isDNA=isDNA)$max-3.4)/(DeltaS(seq, isDNA=isDNA)$max+RlnK))
    nnTM.max <- nnTM.max-272.9+7.21*log(saltConcentration/1000)
    
    nnTM.min <- round(nnTM.min, 1)
    nnTM.max <- round(nnTM.max, 1)
  }
  
  return(list(min=nnTM.min, max=nnTM.max))
}


#' Check self-complementarity of a primer sequence
checkSelfComplementarity <- function(seq) {
  return(1)
}

getIndexOf <- function(seq, subSeq, startIndex, minMatch) {
  ret <- c(-1, -1)
  for (k in minMatch:nchar(subSeq)) {
    theMatch <- regexpr(substr(subSeq, 1, k), substr(seq, startIndex, nchar(seq)), fixed=TRUE)[[1]]
    if (theMatch>0) {
      ret[1] <- theMatch
      ret[2] <- k
    }
  }
  return(ret)
}

#' debugging
checkHairpin <- function(seq, min.hairpin.len=4, bubbleSize=3) {
  seq.rc <- getReverseComplement(seq)
  seq.n <- nchar(seq)
  
  theResult <- c()
  count <- 0
  
  ret <- 0
  theResults <- list()
  
  maxSeqLength <- abs(seq.n/2-bubbleSize)
  for (compPos in 1:(seq.n-2*min.hairpin.len)) {
    maxMatch <- 0
    for (seqPos in 1:(seq.n-maxSeqLength)) {
      theResult <- getIndexOf(substr(seq, 1, seqPos+maxSeqLength-1),
                              substr(seq.rc, compPos, seq.n), 
                              seqPos, min.hairpin.len)
      if (theResult[1]>-1) {
        #print(theResult)
        ret <- 1
        # skip hairpin insert
        a <- theResult[1]
        b <- theResult[1]+theResult[2]-1
        c <- seq.n-compPos-theResult[2]
        d <- seq.n-compPos-1
        theResults <- insertHairpinArray(a, b, c, d, theResults)
        if (theResult[2] > maxMatch) maxMatch <- theResult[2]
        seqPos <- theResult[1] + theResult[2] - min.hairpin.len
        if (seqPos + min.hairpin.len >= maxSeqLength) {
          compPos <- compPos + maxMatch - min.hairpin.len
          break;
        }
      } else {
        if (maxMatch > min.hairpin.len) compPos <- compPos + maxMatch-min.hairpin.len; 
        break;  
      }
      #cat(compPos,", ",seqPos, "\n")
    }
  }
  return(ret)
}

insertHairpinArray <- function(a, b, c, d, results) {
  arrayCount <- length(results)
  if (a>=c | a>=b | c>=d | b>=c) {
    return(results)
  }
  if (arrayCount>0) {
    for (i in 1:arrayCount) {
      if (results[[i]][1]<=a & results[[i]][2]>=b & results[[i]][3]<=c & results[[i]][4]>=d) return(results)
      if (results[[i]][1]>=a & results[[i]][2]<=b & results[[i]][3]>=c & results[[i]][4]<=d) {
        results[[i]][1]=a
        results[[i]][2]=b
        results[[i]][3]=c
        results[[i]][4]=d
        return(results)
      }
    }
  }
  results[[arrayCount+1]] <- rep(0,0,0,0)
  results[[arrayCount+1]][1] <- a
  results[[arrayCount+1]][2] <- b
  results[[arrayCount+1]][3] <- c
  results[[arrayCount+1]][4] <- d
  return(results)
}


checkDimers <- function(seq, min.align.len=5, maxMismatchNum=1) {
  seq.n <- nchar(seq)
  seq.rc <- getReverseComplement(seq)
  
  mat <- matrix(NA, seq.n, seq.n)
  
  rows <- unlist(strsplit(seq.rc, split=""))
  cols <- unlist(strsplit(seq, split=""))
  
  rownames(mat) <- cols
  colnames(mat) <- rows
  
  # assumption, only ACTGU in the sequence
  # fill alignment matrix
  for (i in 1:2) {
    for (j in 1:2) {
      if (cols[i] == rows[j]) {
       # mat[i,j] <- 1
        mat[j,i]=1
        if (i>1 & j>1) mat[j, i] <- mat[j,i]+mat[j-1, i-1] # mat[i, j] <- mat[i,j]+mat[i-1, j-1]
      } else {
        # mat[i,j] <- 0
        mat[j,i] <- 0
      }
    }
  }

  for (i in 3:(seq.n-1)) {
    for (j in 3:(seq.n-1)) {
      if (cols[i]==rows[j]) {
        # mat[i,j] <- mat[i-1, j-1]+1
        mat[j,i] <- mat[j-1, i-1]+1
      } else {
        
        mat[j,i] <- 0 # mat[i,j] <- 0
        # if (!is.na(mat[i-1,j-1])) {
        if (!is.na(mat[j-1,i-1])) {
          # if ((mat[i-1,j-1]>1) & (cols[i+1]==rows[j+1])) {
          if ((mat[j-1,i-1]>1) & (cols[i+1]==rows[j+1])) {
            # mat[i,j] <- mat[i-1, j-1]
            mat[j, i] <- mat[j-1, i-1]
          }
        }
      }
    }
  }
  
  
  i=seq.n
  j=i
  if (cols[i]==rows[j]) {
    #mat[i,j]=1
    #mat[i,j] <- mat[i,j] + mat[i-1,j-1]
    mat[j, i]=1
    mat[j, i] <- mat[j, i] + mat[j-1, i-1]
  } else {
    #mat[i,j] <- 0
    mat[j, i] <- 0
  }
  
  mat[is.na(mat)] <- 0
  
  # Make aligned array
  count <- 1
  theResults <- list()
  for (i in 1:seq.n) {
    for (j in 1:seq.n) {
      #if (mat[i,j]==1) {
      if (mat[j,i]==1) {
        mismatches <- 0
        hasMatch <- 1
        lastMatch <- 1
        maxInc <- seq.n - ifelse(i<=j, j, i)
        k=1
        if (maxInc>0) {
          for (k in 1:(maxInc)) {
            # hasMatch <- mat[i+k, j+k]
            hasMatch <- mat[j+k, i+k]
            if (!hasMatch) break
            if (hasMatch<=lastMatch) {
              if (mismatches>=maxMismatchNum) break
              mismatches <- mismatches + 1
            }
            lastMatch <- hasMatch
          }
        } 
        if ((k-mismatches) >= min.align.len) {
          theResults[[count]] <- c(i, j, i+k-1, j+k-1, mismatches)
          count <- count + 1
        }
      }
    }
  }
  
 # return(theResults)
  
  return(length(theResults))
}

args<-commandArgs(trailingOnly=TRUE)
indata=args[1]
if (is.na(indata)) indata <- "./input/input.txt"
findPrimers(indata)




gpc.power <- function(Dprime,p,q,prevalence,R1,R2,n.cases,n.controls,alpha,power,
                      show.html=TRUE){
# get power for given alpha using Genetic Power Calculator
# wget http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi?fA=0.1&k=0.13&rAa=1.0&rAA=1.2&dprime=0.5&m1=0.3&n=10000000&ccratio=1.5&alpha=5.3526e-10&power=0.8 -O 1.html

#browser()

gpc.url <- "http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi"

cat(paste("connecting to Genetic Power Calculator at\n\t", gpc.url))

query <- paste("?fA=", q, "&k=", prevalence, "&rAa=", R1, 
                 "&rAA=", R2, "&dprime=", Dprime, "&m1=", p, "&n=", format(n.cases, scientific=FALSE),
                 "&ccratio=", n.controls/n.cases, "&alpha=", alpha, "&power=", 
                 power, sep = "")

tfile <- tempfile("GeneticPowerCalculator")
tfile1 <- paste(tfile,"_1.html",sep="")
tfile2 <- paste(tfile,"_2.html",sep="")
wget.cmd <- paste("wget \"",gpc.url,query,"\" -O ", tfile1, " ;", sep="")
wget.cmd
system(wget.cmd)

cat("after wget\n")
system("echo done");


if(show.html)
  system(paste("firefox", tfile1))

system(paste("nohup ./gpc-get-power-from-html.sh ", tfile1," |grep -v \"^cat\" > ", tfile2))
cat("after ./gpc-get-power-from-html.sh \n")

power.res <- read.table(tfile2)
dimnames(power.res) <- list(c("dominant","recessive","general (2df)", "allelic"),
                            c("alphac","power",paste("ncases(pwr=", round(100*power,1),"%)",sep="")))

power.res
}


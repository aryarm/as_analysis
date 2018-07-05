allele.imbalance<- function(rna, rna.err){
   
   #Step1
   # estimate sequencing error, dispersion DNA counts and dispersion RNA counts
   # optimization function 1
   mtmp <- function(par, rna.x, rna.n, rna.err, ge){
     
      # define parameters
      disp.rna <- par
   
      # likelihood given genotype err (not actually het)
      term1<- 0.5 * dbetabinom(rna.x, rna.n, m= rna.err/3, s= disp.rna)
      term2<- 0.5 * dbetabinom(rna.n-rna.x, rna.n, m= rna.err/3, s= disp.rna)
      # genotyping error likelihood for a single site, multiplied by its genotyping error
      err.like <- (ge)*(term1 + term2)

      # likelihood given het sites (actually het)
      het.like = (1-ge)*dbetabinom(rna.x, rna.n, m= 0.5, s= disp.rna)
     
      # log likelihood
      ll<- -sum(log(err.like + het.like))
     
      return(ll)
    }
   
   # optimize
   m0<- optim(par= 1, mtmp,
              rna.x= rna$ref.matches, rna.n = rna$N,
              ge= rna$genotype.error, rna.err= rna.err,
              method="Brent",lower = 1e-5,
              upper = 100)
   
   # get parameter
   disp.rna = m0$par
   
   # get genes
   genes= as.character(unique(rna$gene))
   
   # allelic.imbalance  is how far proportion is away from 0.5
   # first site : L( x | allelic.imbalance)=het.like* P( x | 0.5 + allelic.imbalance) + err.like
   # subsequent sites : L( x | allelic.imbalance) = het.like(0.5 * Pr(x | 0.5 - d) + 0.5 * Pr( x | 0.5 + d))) + err.like
   # chi-square test against null hypothesis of allelic.imbalance = 0
   
   myfunc<- function(i){
      #Step 1 subset by gene
      # rna
      rna= rna[rna$gene %in% genes[i],]
      rna= rna[order(rna$start),]
      
      # step 2 likelihood estimate
      # likelihood function estimate for ase
      # ll<- function(allelic.imbalance,  dna.x, dna.n, rna.x, rna.n, ge, dna.err, rna.err, disp.dna, disp.rna){
      ll<- function(allelic.imbalance,  rna.x, rna.n, ge, rna.err, disp.rna){
         
         # for 1st site
         # likelihood of allelic imbalance (d1)
         d1= dbetabinom(rna.x[1], rna.n[1], m= 0.5+ allelic.imbalance, s= disp.rna) 
         
         # likelihood given genotype err (not actually het)
         term1<- 0.5 * dbetabinom(rna.x[1], rna.n[1], m= rna.err/3, s= disp.rna)
         term2<- 0.5 * dbetabinom(rna.n[1]-rna.x[1], rna.n[1], m= rna.err/3, s= disp.rna)
         err.like <- (ge[1])*(term1 + term2)
         
         # likelihood of 1st sites
         p1= d1+ err.like
         
         # for subsequent sites
         len= nrow(rna)
         
         if(len >1){
            
            # likelihood of allelic imbalance of n+1, n+2....n+k sites
            d2= 0.5 * dbetabinom(rna.x[2:len], rna.n[2:len], 0.5+ allelic.imbalance, disp.rna) 
            d3= 0.5 * dbetabinom(rna.x[2:len], rna.n[2:len], 0.5- allelic.imbalance, disp.rna) 
            
            # likelihood given genotype err (not actually het)
            term1<- 0.5 * dbetabinom(rna.x[2:len], rna.n[2:len], m= rna.err/3, s= disp.rna)
            
            term2<- 0.5 * dbetabinom(rna.n[2:len]-rna.x[2:len], rna.n[2:len], m= rna.err/3, s= disp.rna)
            
            err.like <- (ge[2:len])*(term1 + term2)
            
            p2= (d2 + d3) + err.like
            
            return(-sum(log(c(p1,p2))))
          }else{
               return(-sum(log(p1)))
          }
      }
      
      # optimize allelic.imbalance
      m1= optimize(ll, c(-0.5, 0.5),
                   rna.x= rna$ref.matches, rna.n= rna$N,
                   ge= rna$genotype.error,
                   rna.err= rna.err, disp.rna= disp.rna)
      
      # estimates of allelic.imbalance
      alt.ll <- m1$objective
      estimate <- m1$minimum
      
      # NULL hypothesis
      null.ll= ll(allelic.imbalance = 0,
                  rna.x= rna$ref.matches, rna.n= rna$N,
                  ge= rna$genotype.error, rna.err= rna.err,
                  disp.rna= disp.rna)
      
      # Likelihood ratio test
      lrt.stat <- 2 * (null.ll - alt.ll)
      pval <- pchisq(lrt.stat, df=1, lower.tail=F)
      result= data.frame(gene= genes[i], pval =pval, d= estimate)
      return(result)
   }
   rlist= llply(1:length(genes), myfunc)
   res= do.call(rbind,rlist)
   res$fdr= p.adjust(res$pval, "fdr")
   return(res)
}
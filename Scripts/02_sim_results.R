### Bayesian analysis with longitudinal studies by indication
### Reagan Mozer and Mark E. Glickman


### This script aggregates the results of the simulation study and creates Table 2 in the manuscript.


# Read in results and make table
library(tidyverse)
library(knitr)
library(kableExtra)

add_res = function(m1){
  tmp=data.frame(m1$means[1:5,],n=n,pX=pX,seed=this.seed)
  names(tmp)=c("Target","Av.Median","Av.Mean","Av.Lower95%CI","Av.Upper95%CI","Av.Range95%CI",
               "Coverage","Av.AutoCorr(Lag10)","nsim","n","pX","seed")
  tmp$param = c("q1","q2","p1","p2","p3")
  
  tmp1 = tmp %>% select(n, pX, seed, param, nsim, Av.Median, Coverage) %>% 
    pivot_wider(names_from=param, values_from=c("Av.Median","Coverage"))
  return(tmp1)
}

f=list.files("Results/")
res = data.frame()

for (j in 1:length(f)){
  load(f[j])
  curr = add_res(m1)
  res=rbind(res,curr)
  gdata::keep(res, f, add_res,sure=T)
}



out = res %>% select(n, pX, seed, starts_with("Av.Median"), starts_with("Coverage"))
names(out)=gsub("Av.Median","est",names(out))
names(out)=gsub("Coverage","cov",names(out))


out2=out %>%select(n, pX, ends_with("q1"),ends_with("q2"),ends_with("p1"),ends_with("p2"),ends_with("p3"))
tab2 = knitr::kable(out2, format="latex",digits=2,row.names=F, col.names=c("","",rep(c("Est.","CR"),5)), 
                    linesep=c(rep("",3),"\\hline"),escape=F,
                    caption="Simulation results")%>% add_header_above(c("$n$", "$\\pi_X$", "$q_1$"=2, "$q_2$"=2,
                                                                        "$p_1$"=2,"$p_2$"=2,"$p_3$"=2),escape=F)



sink("Tables/sim_results.tex")
tab2
sink()
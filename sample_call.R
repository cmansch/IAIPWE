for (file in list.files("./value_estimator", pattern=".R")){
  source(paste0("./value_estimator/", file)) # these include the data settings
}
q2 <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + I(1-a1):a2 +
                                I(1-a1):x21 + I(1-a1):a2:x21 + # for a1==0 arms
                                a1 + a1:x11 + a1:x12 + a1:a2 + # for a2==0 arms
                                a1:x21 + a1:a2:x21 - 1,
                              solver.method = 'lm',
                              predict.method = 'predict.lm')
 
q2r <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + I(1-a1):x21 +
                                  a1 + x11:a1 + x12:a1 + x21:a1 - 1,
                                solver.method = 'lm',
                                predict.method = 'predict.lm')
q1 <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + 
                                 a1 + x11:a1 + x12:a1 - 1,
                               solver.method = 'lm',
                               predict.method = 'predict.lm')
 
q_list <- list(q1, list("r0" = q2,
                        "r1" = q2r) )
 
p1 <- modelObj::buildModelObj(model = ~ 1,
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))

p2 <- modelObj::buildModelObj(model = ~ I(a1==0):I(r2==0) -1, # + I(a1==0):I(r2==1) + I(a1==1):I(r2==0),
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))
pi_list <- list(p1, p2)

b<-1
df <- design1(vp=1, n=500, seed=b*10000, s2=100, a1p=0.5, a2p=0.5, r2p = 0.4)
regimes <- list(c(0,0), c(0,1), c(1,0), c(1,1))

regime_dfs <- regimelist(embRegimes = regimes,
                           dat = df)

IAIPW(df, pi_list, q_list, regime_dfs, feasibleSetsIndicator=TRUE, t_s=500)

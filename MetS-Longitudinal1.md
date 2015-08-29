``` r
require(coin)
require(compareGroups)
require(ggplot2)
require(epitools)
require(tidyr)
require(dplyr)


load("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/R_analysis/TABLAS/tE.Rdata")



tE=subset(tE,!is.na(gen)&!is.na(num.ms2))
tE=subset(tE,!HIQHEART%in%'Yes') ##quita a los que tienen CVD
########quita a los que estan tomando tratamiento
########de diabetes, HTA o lipidos
##tE=subset(tE,!MEQDIAB%in%'Yes')
##tE=subset(tE,!MEQLIPID%in%'Yes')
##tE=subset(tE,!MEQBLPR%in%'Yes')

####define la variable 'dosis de alelo E4'
tE$e4[tE$gen%in%c('E2E2','E2E3','E3E3')]='00'
####tE$e4[tE$gen%in%c('E3E4')]='01'
tE$e4[tE$gen%in%c('E2E4','E3E4')]='01'
tE$e4[tE$gen%in%c('E4E4')]='11'
tE$e4=as.factor(tE$e4)


tE$num.ms2=ifelse(tE$num.ms2%in%c('3','4','5'),'3+',tE$num.ms2)
tE$num.ms2=as.factor(tE$num.ms2)
tE$HOMA=tE$LBXGLU*tE$LBXIN/405
label(tE$HOMA)='HOMA'

tE$waist.ms=factor(tE$waist.ms,labels=c('No','Yes'))
label(tE$LBXAGE)='Age'

tE[,64:78]=lapply(tE[,64:78],as.factor)

####define la variable E4
tE$e[tE$gen%in%c('E2E2','E2E3','E3E3')]='00'
tE$e[tE$gen%in%c('E2E4','E4E4','E3E4')]='01'
tE$e=as.factor(tE$e)
```

TABLAS BIVARIADAS
=================

MetS
----

``` r
##############################################Grupos MetS
t=compareGroups(metS2~LBXAGE+bmi3cat+BMXWAIST+SMQSTAT+ALDQWK+WKTYPE+LBXTC+LBXTR+LBXHDD+LBXGLU+LBXIN+LBXGH+HOMA+e4+e,data=tE)
```

    ## Warning in compare.i(X[, i], y = y, selec.i = selec[i], method.i =
    ## method[i], : Some levels of 'e' are removed since no observation in that/
    ## those levels

``` r
c=createTable(t,show.ratio =T,hide.no='No')
c
```

    ## 
    ## --------Summary descriptives table by 'metS2'---------
    ## 
    ## ___________________________________________________________________________ 
    ##                     0            1              OR        p.ratio p.overall 
    ##                   N=3147       N=1261                                       
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
    ## Age            47.9 (9.42)  52.5 (5.20)  1.09 [1.07;1.10] <0.001   <0.001   
    ## bmi3cat:                                                           <0.001   
    ##     norm       912 (29.1%)   65 (5.17%)        Ref.        Ref.             
    ##     obese      360 (11.5%)  660 (52.5%)  25.6 [19.5;34.3]  0.000            
    ##     overweight 1861 (59.4%) 533 (42.4%)  4.01 [3.08;5.29]  0.000            
    ## BMXWAIST       94.0 (8.76)   105 (8.74)  1.16 [1.15;1.17] <0.001   <0.001   
    ## SMQSTAT:                                                           <0.001   
    ##     Non-smoker 1118 (35.7%) 376 (30.0%)        Ref.        Ref.             
    ##     Ex-smoker  841 (26.9%)  433 (34.6%)  1.53 [1.30;1.81] <0.001            
    ##     Current    1170 (37.4%) 443 (35.4%)  1.13 [0.96;1.32]  0.147            
    ## ALDQWK         61.9 (57.5)  72.9 (65.6)  1.00 [1.00;1.00] <0.001   <0.001   
    ## WKTYPE:                                                            <0.001   
    ##     H          2732 (86.8%) 1034 (82.1%)       Ref.        Ref.             
    ##     S          415 (13.2%)  226 (17.9%)  1.44 [1.20;1.72] <0.001            
    ## LBXTC           211 (37.0)   220 (39.3)  1.01 [1.00;1.01] <0.001   <0.001   
    ## LBXTR           124 (84.0)   212 (133)   1.01 [1.01;1.01] <0.001   <0.001   
    ## LBXHDD         54.4 (10.8)  47.5 (10.2)  0.93 [0.93;0.94] <0.001   <0.001   
    ## LBXGLU         94.3 (14.1)   110 (24.9)  1.06 [1.06;1.07] <0.001   <0.001   
    ## LBXIN          5.84 (3.77)  10.1 (7.05)  1.22 [1.20;1.25] <0.001   <0.001   
    ## LBXGH          5.39 (0.42)  5.79 (0.81)  4.82 [3.97;5.86] <0.001   <0.001   
    ## HOMA           1.38 (1.00)  2.78 (2.25)  2.40 [2.21;2.60] <0.001   <0.001   
    ## e4:                                                                 0.047   
    ##     00         2583 (82.1%) 1009 (80.0%)       Ref.        Ref.             
    ##     01         541 (17.2%)  234 (18.6%)  1.11 [0.93;1.31]  0.240            
    ##     11          23 (0.73%)   18 (1.43%)  2.01 [1.06;3.74]  0.033            
    ## e:                                                                  0.121   
    ##     00         2583 (82.1%) 1009 (80.0%)       Ref.        Ref.             
    ##     01         564 (17.9%)  252 (20.0%)  1.14 [0.97;1.35]  0.113            
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

MetS~APOE4 presence
-------------------

### Normal weight

``` r
####################################BMI
t1=compareGroups(e4~LBXAGE+bmi+SMQSTAT+ALDQWK+WKTYPE,data=tE,subset=bmi<25)
c1=createTable(t1,show.p.trend=T,hide.no='No')

t3=compareGroups(e4~LBXTC+LBXTR+LBXHDD+LBXGLU+LBXIN+LBXGH+HOMA,
                 data=tE,subset=bmi<25)
c3=createTable(t3,show.p.trend=T,hide.no='No')

t2=compareGroups(e4~waist.ms+hta.ms+glu.ms+tg.ms2+hdl.ms2+num.ms2,data=tE,subset=bmi<25)
c2=createTable(t2,show.p.trend=T,hide.no='0')

tot1=rbind('Anthropometric and Lifestyle'=c1,'Plasma Biochemistry'=c3,'metS2 Criteria'=c2)###une las tablas
tot1
```

    ## 
    ## --------Summary descriptives table by 'e4'---------
    ## 
    ## ________________________________________________________________________ 
    ##                        00          01          11      p.overall p.trend 
    ##                       N=793       N=176        N=8                       
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
    ## Anthropometric and Lifestyle:
    ##     Age            44.9 (10.8) 44.3 (11.4) 49.1 (8.82)   0.436    0.841  
    ##     bmi            23.3 (1.49) 23.2 (1.51) 23.3 (1.50)   0.994    0.966  
    ##     SMQSTAT:                                             0.256    0.579  
    ##         Non-smoker 273 (34.7%) 61 (35.1%)   0 (0.00%)                    
    ##         Ex-smoker  157 (20.0%) 35 (20.1%)   3 (37.5%)                    
    ##         Current    356 (45.3%) 78 (44.8%)   5 (62.5%)                    
    ##     ALDQWK         58.0 (53.9) 60.6 (65.9) 76.2 (44.1)   0.577    0.391  
    ##     WKTYPE:                                              0.789    0.695  
    ##         H          686 (86.5%) 152 (86.4%)  8 (100%)                     
    ##         S          107 (13.5%) 24 (13.6%)   0 (0.00%)                    
    ## Plasma Biochemistry:
    ##     LBXTC          201 (38.6)  209 (37.2)  216 (50.2)    0.032    0.009  
    ##     LBXTR          113 (77.1)  112 (63.9)  123 (44.0)    0.893    0.907  
    ##     LBXHDD         55.8 (11.5) 55.7 (11.9) 54.1 (10.7)   0.906    0.746  
    ##     LBXGLU         93.9 (18.5) 94.4 (17.4) 97.6 (29.9)   0.806    0.589  
    ##     LBXIN          4.66 (2.63) 4.48 (2.30) 5.31 (4.97)   0.573    0.649  
    ##     LBXGH          5.37 (0.59) 5.39 (0.57) 6.00 (1.67)   0.035    0.183  
    ##     HOMA           1.10 (0.76) 1.05 (0.57) 1.57 (2.12)   0.170    0.877  
    ## metS2 Criteria:
    ##     waist.ms:                                            1.000    0.638  
    ##         No         792 (99.9%) 176 (100%)   8 (100%)                     
    ##         Yes         1 (0.13%)   0 (0.00%)   0 (0.00%)                    
    ##     hta.ms         270 (34.0%) 64 (36.4%)   1 (12.5%)    0.389    0.979  
    ##     glu.ms         188 (23.7%) 47 (26.7%)   2 (25.0%)    0.685    0.429  
    ##     tg.ms2         167 (21.1%) 38 (21.6%)   5 (62.5%)    0.030    0.209  
    ##     hdl.ms2        48 (6.05%)   7 (3.98%)   1 (12.5%)    0.250    0.509  
    ##     num.ms2:                                             0.159    0.457  
    ##         0          355 (44.8%) 71 (40.3%)   1 (12.5%)                    
    ##         1          258 (32.5%) 66 (37.5%)   6 (75.0%)                    
    ##         2          127 (16.0%) 28 (15.9%)   0 (0.00%)                    
    ##         3+         53 (6.68%)  11 (6.25%)   1 (12.5%)                    
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

### Overweight

``` r
t1=compareGroups(e4~LBXAGE+bmi+SMQSTAT+ALDQWK+WKTYPE,data=tE,subset=bmi>=25&bmi<30)
c1=createTable(t1,show.p.trend=T,hide.no=0)

t3=compareGroups(e4~LBXTC+LBXTR+LBXHDD+LBXGLU+LBXIN+LBXGH+HOMA,
                 data=tE,subset=bmi>=25&bmi<30)
c3=createTable(t3,show.p.trend=T,hide.no='No')

t2=compareGroups(e4~waist.ms+hta.ms+glu.ms+tg.ms2+hdl.ms2+num.ms2,data=tE,subset=bmi>=25&bmi<30)
c2=createTable(t2,show.p.trend=T,hide.no='0')

tot2=rbind('Anthropometric and Lifestyle'=c1,'Plasma Biochemistry'=c3,'metS2 Criteria'=c2)###une las tablas
tot2
```

    ## 
    ## --------Summary descriptives table by 'e4'---------
    ## 
    ## _________________________________________________________________________ 
    ##                         00          01          11      p.overall p.trend 
    ##                       N=1952       N=422       N=20                       
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
    ## Anthropometric and Lifestyle:
    ##     Age            50.0 (7.78)  50.3 (7.48) 49.1 (8.48)   0.691    0.691  
    ##     bmi            27.4 (1.39)  27.4 (1.35) 27.1 (1.31)   0.381    0.214  
    ##     SMQSTAT:                                              0.463    0.681  
    ##         Non-smoker 665 (34.3%)  148 (35.2%)  7 (35.0%)                    
    ##         Ex-smoker  576 (29.7%)  131 (31.2%)  3 (15.0%)                    
    ##         Current    700 (36.1%)  141 (33.6%) 10 (50.0%)                    
    ##     ALDQWK         66.1 (59.4)  64.2 (57.2) 53.0 (67.0)   0.532    0.357  
    ##     WKTYPE:                                               0.277    0.281  
    ##         H          1656 (84.8%) 345 (81.9%) 18 (90.0%)                    
    ##         S          296 (15.2%)  76 (18.1%)   2 (10.0%)                    
    ## Plasma Biochemistry:
    ##     LBXTC           216 (36.0)  219 (36.9)  227 (41.9)    0.096    0.036  
    ##     LBXTR           149 (106)    156 (103)   193 (160)    0.094    0.065  
    ##     LBXHDD         52.6 (10.7)  51.3 (10.6) 52.6 (10.9)   0.073    0.041  
    ##     LBXGLU         98.1 (16.1)  99.5 (18.1) 98.2 (22.5)   0.289    0.153  
    ##     LBXIN          6.72 (5.01)  7.21 (6.30) 7.43 (4.40)   0.254    0.100  
    ##     LBXGH          5.47 (0.46)  5.47 (0.53) 5.50 (0.61)   0.965    0.872  
    ##     HOMA           1.67 (1.44)  1.84 (1.95) 1.90 (1.50)   0.124    0.043  
    ## metS2 Criteria:
    ##     waist.ms:                                             0.343    0.900  
    ##         No         1577 (80.8%) 347 (82.2%) 14 (70.0%)                    
    ##         Yes        375 (19.2%)  75 (17.8%)   6 (30.0%)                    
    ##     hta.ms         1115 (57.1%) 242 (57.3%) 13 (65.0%)    0.777    0.719  
    ##     glu.ms         730 (37.4%)  180 (42.7%)  7 (35.0%)    0.126    0.085  
    ##     tg.ms2         777 (39.8%)  183 (43.4%) 11 (55.0%)    0.168    0.078  
    ##     hdl.ms2        159 (8.15%)  49 (11.6%)   3 (15.0%)    0.034    0.013  
    ##     num.ms2:                                              0.301    0.022  
    ##         0          337 (17.3%)  68 (16.1%)   2 (10.0%)                    
    ##         1          611 (31.3%)  119 (28.2%)  5 (25.0%)                    
    ##         2          588 (30.1%)  125 (29.6%)  6 (30.0%)                    
    ##         3+         416 (21.3%)  110 (26.1%)  7 (35.0%)                    
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

### Obese

``` r
t1=compareGroups(e4~LBXAGE+bmi+SMQSTAT+ALDQWK+WKTYPE,data=tE,subset=bmi>=30)
c1=createTable(t1,show.p.trend=T,hide.no='No')

t3=compareGroups(e4~LBXTC+LBXTR+LBXHDD+LBXGLU+LBXIN+LBXGH+HOMA,
                 data=tE,subset=bmi>=30)
c3=createTable(t3,show.p.trend=T,hide.no='No')

t2=compareGroups(e4~waist.ms+hta.ms+glu.ms+tg.ms2+hdl.ms2+num.ms2,data=tE,subset=bmi>=30)
c2=createTable(t2,show.p.trend=T,hide.no='0')

tot3=rbind('Anthropometric and Lifestyle'=c1,'Plasma Biochemistry'=c3,'metS2 Criteria'=c2)###une las tablas
tot3
```

    ## 
    ## --------Summary descriptives table by 'e4'---------
    ## 
    ## ________________________________________________________________________ 
    ##                        00          01          11      p.overall p.trend 
    ##                       N=828       N=174       N=13                       
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
    ## Anthropometric and Lifestyle:
    ##     Age            51.5 (6.67) 51.2 (6.64) 50.5 (8.52)   0.758    0.473  
    ##     bmi            32.7 (2.51) 32.7 (2.36) 32.4 (2.66)   0.915    0.770  
    ##     SMQSTAT:                                             0.291    0.592  
    ##         Non-smoker 274 (33.3%) 54 (31.0%)   8 (61.5%)                    
    ##         Ex-smoker  296 (35.9%) 67 (38.5%)   2 (15.4%)                    
    ##         Current    254 (30.8%) 53 (30.5%)   3 (23.1%)                    
    ##     ALDQWK         69.0 (64.6) 72.1 (67.3) 62.3 (72.5)   0.785    0.755  
    ##     WKTYPE:                                              0.628    0.513  
    ##         H          717 (86.6%) 155 (89.1%) 11 (84.6%)                    
    ##         S          111 (13.4%) 19 (10.9%)   2 (15.4%)                    
    ## Plasma Biochemistry:
    ##     LBXTC          215 (38.7)  222 (39.9)  222 (29.8)    0.053    0.020  
    ##     LBXTR           178 (127)   187 (114)   236 (186)    0.186    0.123  
    ##     LBXHDD         49.8 (10.6) 48.6 (9.71) 40.7 (13.0)   0.004    0.008  
    ##     LBXGLU         105 (25.1)  103 (18.0)  109 (19.3)    0.481    0.593  
    ##     LBXIN          9.90 (5.96) 9.73 (5.78) 9.85 (6.22)   0.948    0.774  
    ##     LBXGH          5.72 (0.78) 5.60 (0.48) 5.73 (0.52)   0.169    0.117  
    ##     HOMA           2.62 (1.92) 2.57 (1.91) 2.77 (1.94)   0.920    0.934  
    ## metS2 Criteria:
    ##     waist.ms:                                            0.670    0.487  
    ##         No         104 (12.6%) 25 (14.4%)   2 (15.4%)                    
    ##         Yes        724 (87.4%) 149 (85.6%) 11 (84.6%)                    
    ##     hta.ms         603 (72.8%) 124 (71.3%) 12 (92.3%)    0.285    0.695  
    ##     glu.ms         444 (53.6%) 90 (51.7%)   9 (69.2%)    0.468    0.883  
    ##     tg.ms2         467 (56.4%) 101 (58.0%)  9 (69.2%)    0.612    0.433  
    ##     hdl.ms2        126 (15.2%) 31 (17.8%)   3 (23.1%)    0.454    0.272  
    ##     num.ms2:                                             0.742    0.902  
    ##         0          19 (2.29%)   4 (2.30%)   1 (7.69%)                    
    ##         1          73 (8.82%)  16 (9.20%)   1 (7.69%)                    
    ##         2          200 (24.2%) 43 (24.7%)   1 (7.69%)                    
    ##         3+         536 (64.7%) 111 (63.8%) 10 (76.9%)                    
    ## ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

``` r
totr=cbind('Normal Weight'=tot1,'Overweight'=tot2,'Obese'=tot3)
```

ODDS RATIO
==========

### MetS~APOE4 MODELO 1. Todos los BMIs juntos

``` r
or1=glm(metS2~e,data=tE,family='binomial')
m1=cbind(round(exp(coef(or1)),2),round(exp(confint(or1)),2))
```

    ## Waiting for profiling to be done...

``` r
or2=glm(metS2~e+LBXAGE+SMQSTAT+WKTYPE,data=tE,family='binomial')##MODELO 1. Todos los BMIs juntos
m2=cbind(round(exp(coef(or2)),2),round(exp(confint(or2)),2))
```

    ## Waiting for profiling to be done...

``` r
rbind('sin ajustar'=m1[2,],'ajustado'=m2[2,])
```

    ##                  2.5 % 97.5 %
    ## sin ajustar 1.14  0.97   1.35
    ## ajustado    1.16  0.97   1.37

### MetS~APOE4 MODELO 2. BMIs separados

``` r
#####Al ajustar por LBXAGE y WKTYPE se pierde la significacion
#####ajustando por SMQSTAT~sin ajustar
or3=glm(metS2~e+SMQSTAT+LBXAGE+WKTYPE,data=tE,subset=bmi<25,family='binomial')
m3=cbind(round(exp(coef(or3)),2),round(exp(confint(or3)),2))
```

    ## Waiting for profiling to be done...

``` r
or4=glm(metS2~e+SMQSTAT+LBXAGE+WKTYPE,data=tE,subset=bmi>=25&bmi<30,family='binomial')
m4=cbind(round(exp(coef(or4)),2),round(exp(confint(or4)),2))
```

    ## Waiting for profiling to be done...

``` r
or5=glm(metS2~e+SMQSTAT+LBXAGE+WKTYPE,data=tE,subset=bmi>=30,family='binomial')
m5=cbind(round(exp(coef(or5)),2),round(exp(confint(or5)),2))
```

    ## Waiting for profiling to be done...

``` r
rbind('normopeso'=m3[2,],'sobrepeso'=m4[2,],'obeso'=m5[2,])
```

    ##                2.5 % 97.5 %
    ## normopeso 0.96  0.48   1.81
    ## sobrepeso 1.33  1.04   1.69
    ## obeso     1.05  0.74   1.48

RISK RATIO
==========

``` r
tE$bmiTert=factor(tE$bmi3cat,levels=c("norm","overweight","obese"),labels=c('1st','2nd','3rd'))
```

### Normalweight

``` r
z1=subset(tE,tE$bmiTert=='1st')
tab1=with(z1,table(e4,metS2))
round(prop.table(tab1,1),3)*100
```

    ##     metS2
    ## e4      0    1
    ##   00 93.3  6.7
    ##   01 93.8  6.2
    ##   11 87.5 12.5

``` r
cmh_test(metS2~e4,data=z1,scores=list(e4=0:2))##chi2 con modificacion de cmh
```

    ## 
    ##  Asymptotic Linear-by-Linear Association Test
    ## 
    ## data:  metS2 by e4 (00 < 01 < 11)
    ## Z = -0.069526, p-value = 0.9446
    ## alternative hypothesis: two.sided

``` r
riskratio.small(tab1)[2]
```

    ## Warning in chisq.test(xx, correct = correction): Chi-squared approximation
    ## may be incorrect

    ## $measure
    ##     risk ratio with 95% C.I.
    ## e4    estimate     lower     upper
    ##   00 1.0000000        NA        NA
    ##   01 0.9189815 0.4901685  1.722932
    ##   11 1.8379630 0.2884959 11.709380

### Overweight

``` r
z2=subset(tE,tE$bmiTert=='2nd')
tab2=with(z2,table(e4,metS2))
round(prop.table(tab2,1),3)*100
```

    ##     metS2
    ## e4      0    1
    ##   00 78.7 21.3
    ##   01 73.9 26.1
    ##   11 65.0 35.0

``` r
cmh_test(metS2~e4,data=z2,scores=list(e4=0:2))
```

    ## 
    ##  Asymptotic Linear-by-Linear Association Test
    ## 
    ## data:  metS2 by e4 (00 < 01 < 11)
    ## Z = -2.5004, p-value = 0.0124
    ## alternative hypothesis: two.sided

``` r
riskratio.small(tab2)[2]
```

    ## Warning in chisq.test(xx, correct = correction): Chi-squared approximation
    ## may be incorrect

    ## $measure
    ##     risk ratio with 95% C.I.
    ## e4   estimate     lower    upper
    ##   00 1.000000        NA       NA
    ##   01 1.220805 1.0177721 1.464341
    ##   11 1.639209 0.8966512 2.996711

### Obese

``` r
z3=subset(tE,tE$bmiTert=='3rd')####utilizar solo si se definen 3 categorias de peso
tab3=with(z3,table(e4,metS2))
round(prop.table(tab3,1),3)*100
```

    ##     metS2
    ## e4      0    1
    ##   00 35.3 64.7
    ##   01 36.2 63.8
    ##   11 23.1 76.9

``` r
cmh_test(metS2~e4,data=z3,scores=list(e4=0:2))
```

    ## 
    ##  Asymptotic Linear-by-Linear Association Test
    ## 
    ## data:  metS2 by e4 (00 < 01 < 11)
    ## Z = -0.24306, p-value = 0.808
    ## alternative hypothesis: two.sided

``` r
riskratio.small(tab3)[2]
```

    ## Warning in chisq.test(xx, correct = correction): Chi-squared approximation
    ## may be incorrect

    ## $measure
    ##     risk ratio with 95% C.I.
    ## e4   estimate     lower    upper
    ##   00 1.000000        NA       NA
    ##   01 0.985249 0.8715150 1.113826
    ##   11 1.188034 0.8784164 1.606784

``` r
#################### 3 figura
tab3=rbind(prop.table(tab1,1),prop.table(tab2,1),prop.table(tab3,1))
tabss=data.frame(bmi=c(rep('Norm.',3),rep('Overweight',3),rep('Obese',3)),
                 gene=rownames(tab3),
                 YES=as.numeric(tab3[,2]),
                 NO=as.numeric(tab3[,1]))
tabss$gene=factor(tabss$gene,levels=c('00','01','11'),labels=c('_,_','_,E4','E4,E4'))
tabss$bmi=factor(tabss$bmi,levels=c('Norm.','Overweight','Obese'))
tabss=gather(tabss,"metS2",'n',YES:NO)

c=ggplot(tabss,aes(x=gene,y=n,fill=metS2))+
  geom_bar(stat='identity',position='fill')+
  scale_fill_manual(values = c("black", "grey"))+
  ylab('MetS (frequency)')+
  coord_cartesian(ylim=c(0,1))+
  
  facet_wrap(~bmi)+
  
  theme_bw() +
  theme(axis.title.y = element_text(face="bold", size=rel(1.5),vjust=1),
        axis.title.x = element_text(face="bold", size=rel(1.5),vjust=-1),
        axis.text = element_text(size = rel(1.3),colour='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(colour=NA, size=16, face="bold"),
        legend.text = element_text(colour="black", size = 15),
        strip.text = element_text(size = 20, colour = "black"),
        aspect.ratio=1.2)

c
```

![](scriptE4-bmi3cat_files/figure-markdown_github/unnamed-chunk-12-1.png)

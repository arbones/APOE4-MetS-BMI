Data frame sin purgar
---------------------

### tipo de trabajo

``` r
wt <- read.delim("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/datos.asociados/R/WKTYPE_oct.txt")
```

### primera cesion de datos (Oct/2012)

``` r
oct <- read.delim("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/datos.asociados/R/datosAsociadosMuestrasCesionJMArbones20120110.txt")

oct2=merge(oct,wt,by='CODALIC1',all.x=T)
oct3=subset(oct2,select=-c(CODPLAC2,CODALIC2))
```

### segunda cesion de datos (Julio/2013)

``` r
julio <- read.delim("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/datos.asociados/R/datosAsociadosMuestrasCesionJMArbones_20130731.txt")
```

### Merge las dos dataframes

### t1--\>data frame sin purgar

``` r
t1=rbind(oct3,julio)
```

Carga data frame ya purgada
---------------------------

NO HAY RECORDS DEL PROCESO DE PURGA !!!!!

``` r
t1p <- read.csv("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/datos.asociados/R/t1p.txt", sep="", stringsAsFactors=FALSE)

###t1p-->t1 purgado###
t1p$BMXWAIST[t1p$BMXWAIST>150]=NA
t1p$BMXWT[t1p$BMXWT>500]=NA

###t1m-->t1p modificado con BMIs. Numerico y categorias###
###nuevas variables
t1m=within(t1p,{bmi=BMXWT/(BMXHT/100)^2})
t1m$bmi3cat[t1m$bmi<25]='norm'
t1m$bmi3cat[t1m$bmi>=25 & t1m$bmi<30]='overweight'
t1m$bmi3cat[t1m$bmi>=30]='obese'
t1m$bmi3cat=factor(t1m$bmi3cat,levels=c('norm','overweight','obese'))
```

CRITERIOS DE SINDROME METABOLICO
================================

### t1ms--\>t1m modificado con MetS

### criterios de MetS (AHA/NHLBI-2004)\>\>

### <http://circ.ahajournals.org/content/112/17/2735.full>

``` r
t1ms <- read.csv("~/Dropbox/SyncBriefcase/LAB/IACS/AWHS/datos.asociados/R/t1ms.txt", sep="")
```

### ms\>\>\>hipolipemiantes contabilizados en el criterio de HDL

``` r
t1ms$waist.ms=ifelse(t1ms$BMXWAIST>=102,1,0)
t1ms$hta.ms=ifelse(t1ms$BPDSYAV<130 & t1ms$BPDDIAV<85 & (is.na(t1ms$MEQBLPR) | t1ms$MEQBLPR==0),0,1)
t1ms$glu.ms=ifelse(t1ms$LBXGLU<100& (is.na(t1ms$MEQDIAB) | t1ms$MEQDIAB==0),0,1)
t1ms$tg.ms=ifelse(t1ms$LBXTR>=150,1,0)
t1ms$hdl.ms=ifelse(t1ms$LBXHDD>=40 & (is.na(t1ms$MEQLIPID) | t1ms$MEQLIPID==0),0,1)
###Numero de criterios
t1ms$num.ms=t1ms$waist.ms+t1ms$hta.ms+t1ms$glu.ms+t1ms$tg.ms+t1ms$hdl.ms
###MetS>>yes/no (1/0) 
t1ms$metS=ifelse(t1ms$num.ms>=3,1,0)
```

### ms3\>\>\>hipolipemiantes contabilizados en el criterio de TRIGLICERIDOS

``` r
t1ms$tg.ms2=ifelse(t1ms$LBXTR<=150& (is.na(t1ms$MEQLIPID) | t1ms$MEQLIPID==0),0,1)
t1ms$hdl.ms2=ifelse(t1ms$LBXHDD<40,1,0)
###Numero de criterios
t1ms$num.ms2=t1ms$waist.ms+t1ms$hta.ms+t1ms$glu.ms+t1ms$tg.ms2+t1ms$hdl.ms2
###MetS>>yes/no (1/0) 
t1ms$metS2=ifelse(t1ms$num.ms2>=3,1,0)
```

### ms3\>\>\>hipolipemiantes contabilizados en el criterio de TRIGLICERIDOS

``` r
t1ms$tg.ms3=ifelse(t1ms$LBXTR>=150,1,0)
t1ms$hdl.ms3=ifelse(t1ms$LBXHDD<40,1,0)
###Numero de criterios
t1ms$num.ms3=t1ms$waist.ms+t1ms$hta.ms+t1ms$glu.ms+t1ms$tg.ms3+t1ms$hdl.ms3
###MetS>>yes/no (1/0) 
t1ms$metS3=ifelse(t1ms$num.ms3>=3,1,0)
######################################################
```

``` r
######################################################log2
###Buscar y eliminar duplicados
duplicated(r008$id)
apoe2=r008[!duplicated(r008$id),]


###Merge con t1ms
tE=merge(t1ms,apoe2,by='CODALIC1',all.x=T)

###categorizar APOE
###E2,E3,E4
tE$apoCat[tE$genot=='E2E2'| tE$genot=='E2E3']='E2'
tE$apoCat[tE$genot=='E3E3']='E3'
tE$apoCat[tE$genot=='E3E4'| tE$genot=='E4E4' |tE$genot=='E4/??']='E4'

tE$apoCat=as.factor(tE$apoCat)

###E3,noE3
tE$apoCat2[tE$genot=='E3E3']='E3'

###mas outliers quitados
tE$bmi[tE$bmi>50]=NA
tE$BMXWAIST[tE$BMXWAIST<60]=NA

###Numero de alteraciones metabolicas
#=metS menos la condicion de cintura
tE$alter=tE$hta.ms+tE$glu.ms+tE$tg.ms+tE$hdl.ms
```

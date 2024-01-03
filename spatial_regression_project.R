library(rgdal)
library(sp)
library(spdep)
library(geosphere)
library(dplyr)
library(spatialreg)
library(coda)

setwd("~/spatial_regression_project")

###Wczytywanie podstawowych danych
#Podstawowa mapa NUTS 3
Francja <- readOGR(".", "NUTS_RG_20M_2021_3035")
Francja <- spTransform(Francja, "+proj=longlat")
Francja@data$NUTS_ID_char <- as.character(Francja@data$NUTS_ID)
Francja@data$country <- substr(Francja@data$NUTS_ID_char, 1, 2) 
Francja <- Francja[Francja@data$country == "FR", ]

#Mapa NUTS 2
Francja <- Francja[nchar(Francja@data$NUTS_ID_char)==4, ]

#Mapa NUTS 2 bez terytoriów zamorskich
Francja <- Francja[!grepl("^FRY", Francja@data$NUTS_ID_char), ]

#Wczytywanie danych o gestosci torów kolejowych
trans <- read.csv("tran_r_net__custom_5215623_page_linear.csv")
trans <- trans[!grepl("^FRY", trans$geo), ]
trans <- select(trans, 'geo','OBS_VALUE')
names(trans)[names(trans)=='OBS_VALUE']<-'Tory'

#Laczenie danych i tworzenie mapy
dane <- merge(y = trans, x = Francja, by.y = "geo", by.x = "NUTS_ID_char")
kolorek <- colorRampPalette(c("white", "blue"), bias = 1)
spplot(dane, zcol = 'Tory', colorkey = TRUE, col.regions = kolorek(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Dlugosc torów w kilometrach na 1000 kilometrów kwadratowych we Francji wg. NUTS 2")

#Wczytywanie dodatkowych danych
tour <- read.csv("C:/Users/User/Desktop/Szkola/Ekonometria_przestrzenna/Pd/tour_cap_nuts2_page_linear.csv")
tour <- select(tour, 'geo','OBS_VALUE')
names(tour)[names(tour)=='OBS_VALUE']<-'Turysci'
dane <- merge(y = tour, x = dane, by.y = "geo", by.x = "NUTS_ID_char")

motorway <- read.csv("C:/Users/User/Desktop/Szkola/Ekonometria_przestrzenna/Pd/tran_r_net__custom_5335089_page_linear.csv")
motorway <- select(motorway, 'geo','OBS_VALUE')
names(motorway)[names(motorway)=='OBS_VALUE']<-'Drogi'
dane <- merge(y = motorway, x = dane, by.y = "geo", by.x = "NUTS_ID_char")

PKB <- read.csv("C:/Users/User/Desktop/Szkola/Ekonometria_przestrzenna/Pd/nama_10r_2gdp_page_linear.csv")
PKB <- select(PKB, 'geo','OBS_VALUE')
names(PKB)[names(PKB)=='OBS_VALUE']<-'PKB'

perCapita <- read.csv("C:/Users/User/Desktop/Szkola/Ekonometria_przestrzenna/Pd/nama_10r_3popgdp_page_linear.csv")
perCapita <- select(perCapita, 'geo','OBS_VALUE')
names(perCapita)[names(perCapita)=='OBS_VALUE']<-'Ludnosc'
PKBperCapita <- merge(y = PKB, x = perCapita, by.y = "geo", by.x = "geo")
PKBperCapita$PKBpc<-PKBperCapita$PKB/PKBperCapita$Ludnosc
PKBperCapita <- select(PKBperCapita, 'geo','PKBpc')
dane <- merge(y = PKBperCapita, x = dane, by.y = "geo", by.x = "NUTS_ID_char")

###Tworzenie roznych macierzy W
#Normalizacja wartoscia wlasna
one_neigh <- poly2nb(dane, queen = T)
norm <- nb2listw(one_neigh, style = "B", zero.policy = TRUE)
norm <- listw2mat(norm)
eig <- eigen(norm)
max_eig <- max(eig$values)
norm <- norm/max_eig 
norm <- mat2listw(norm, style='B') #finalna macierz

#Odwrócone kwadraty odleglosci
centroids <- coordinates(dane)
distance <- distm(centroids, fun = distCosine) / 1000
rownames(distance) <- dane@data$FID
colnames(distance) <- dane@data$FID

for (i in 1:nrow(distance)) {
  for (j in 1:ncol(distance)) {
    if (distance[i,j] >= 200) {
      distance[i,j] <- 0
    }
  }
}

gamma <- 2
odw_odl <- 1 / (distance ^ gamma)
odw_odl[odw_odl==Inf] <- 0 
odw_odl <- mat2listw(odw_odl, style='B') #finalna macierz

#Odleglosc euklidesowa
dane@data$Turysci <- (dane@data$Turysci-mean(dane@data$Turysci))/sd(dane@data$Turysci)
dane@data$Drogi <- (dane@data$Drogi-mean(dane@data$Drogi))/sd(dane@data$Drogi)
dane@data$PKBpc <- (dane@data$PKBpc-mean(dane@data$PKBpc))/sd(dane@data$PKBpc)

dist_euk <- matrix(0, nrow = nrow(distance), ncol = nrow(distance))

for(i in 1:nrow(distance)) {
  for(j in 1:nrow(distance)) {
    dist_euk[i, j] <- sqrt(((dane@data$Turysci[i]-dane@data$Turysci[j])^2)+((dane@data$Drogi[i]-dane@data$Drogi[j])^2)+((dane@data$PKBpc[i]-dane@data$PKBpc[j])^2))
  }
} 

dist_euk <- mat2listw(dist_euk, style='B')#finalna macierz

###Testowanie na obecnosc efektow przestrzennych
#Uzyta macierz W: macierz znormalizowana najwieksza wartoscia wlasna (norm)
model1 <- lm(dane@data$Tory ~ dane@data$Turysci+dane@data$PKBpc) #2 zmienne objasniajace ze wzgledu na próbe 22 regionów
summary(model1)

#Ocena wizualna
dane$reszty <- model1$residuals
paleta <- colorRampPalette(c("green", "white", "blue"), bias = 0.95)
spplot(dane, zcol = "reszty", colorkey = TRUE, col.regions = paleta(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty")

#Globalny test Morana
lm.morantest(model1, norm, zero.policy = TRUE, alternative="greater") 
#H0 odrzucona dla alfa=~0.01 czyli istnieje ujemna autokorelacja przestrzenna.

moran.plot(dane@data$reszty, norm, ylab = "Opóznienie przestrzenne reszt: W*e", xlab = "Reszty: e", pch = 20, main = "Wykres Morana", col = "blue", zero.policy = TRUE)

#Lokalny test Morana
LM <- localmoran(dane@data$reszty, norm, zero.policy = TRUE)
LM[,5] <- p.adjustSP(LM[,5], one_neigh, method = "bonferroni")
LM
#Tylko w jednym obszarze mozna przy poziomie alfa=0.1 odrzucic H0 o braku autokorelacji przestrzennej

#Test Gearyego
geary.test(dane@data$reszty, norm, zero.policy = TRUE)
#H0 odrzucona dla alfa=0.01 przy statystyce C>0 czyli istnieje ujemna autokorealcja przestrzenna

#Test liczby polaczen
summary(dane@data$PKBpc)
#mediana=-0.219 => -0.2 - punkt wyjscia
prog <- -0.2
dane$PKBbin <- ifelse(dane$PKBpc >= prog, 1, 0)
joincount.test(as.factor(dane$PKBbin>0), listw=norm, zero.policy = TRUE)

#df<-data.frame(progi=c(), pvalues=c())
#df[11, 'progi']<-0.4
#df[11, 'pvalue']<- 0.6941
#dput(df)

#Po wpisaniu do ramki danych wyników 11 kolejnych testów dla róznych progów, otrzymano ramke danych:
pvalues<- structure(list(pvalue = c(0.8697, 0.8208, 0.8208, 0.7209, 0.9678, 0.9236, 0.9123, 0.8047, 0.594, 0.2881, 0.6941), progi = c(-0.2, -0.3, -0.4, -0.5, -0.6, -0.1, 0, 0.1, 0.2, 0.3, 0.4)), row.names = c(NA, 11L), class = "data.frame")
pvalues<- pvalues[order(pvalues$progi),]
plot(pvalues$progi, pvalues$pvalue, main='Wartosci p-value testu liczby polaczen dla róznych progów podzialu zmiennej PKBpc', xlab='Progi', ylab='P-values', type='l')

###Dodanie dodatkowej zmiennej - punkty POI w kazdym regionie
#Dane zostaly pobrane dla kazdego regionu oddzielnie oraz wypakowane do jednego folderu o nazwie "france".
#Przedstawiona petla dziala na folderach o nazwach "nazwa-regionu-latest-free.shp" które zawieraja dane shapefile dla odpowiedniego regionu
#i zwraca ramke danych z nazwa regionu w jednej kolumnie i liczba odpowiednich punktów poi w drugiej.
setwd("~france")
lista_plikow<-list.files(getwd())
pois<-data.frame('nazwa'=character(), 'liczba_obiektow'=numeric())

for (path in lista_plikow) {
  setwd(paste("~france/", path, sep=''))
  temp<-readOGR(".","gis_osm_pois_free_1")
  temp <- temp[grepl('attraction|museum|monument|memorial|art|viewpoint|castle', temp@data$fclass), ]
  len<-length(temp@data$code)
  pois<-rbind(pois, c(path, len))
}

#Poniewaz nazwy regionów sa zapisane w formie takiej, jak ustalilo geofabrik (niektóre nazwy skrócone, brak francuzkich symboli) wymagana jest reczna zmiana kodowania regionów.
pois[1]<-c('FRF1','FRI1','FRK1','FRD1','FRC1','FRH0','FRB0','FRF2','FRM0','FRC2','FRD2','FR10','FRJ1','FRI2','FRF3','FRJ2','FRE1','FRG0','FRE2','FRI3','FRL0','FRK2')
names(pois)<-c('code', 'number_of_pois')
dane <- merge(y = pois, x = dane, by.y = "code", by.x = "NUTS_ID_char")
dane@data$number_of_pois<-as.numeric(dane@data$number_of_pois)
spplot(dane, zcol = 'number_of_pois', colorkey = TRUE, col.regions = kolorek(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent'), fontsize = list(text = 20)),
       main = "Liczba punktów poi dla turystów we Francji wg. regionów NUTS 2")

model2<-lm(dane$Turysci~dane$number_of_pois)
summary(model2)

###Modele SAR, SEM i SLX
#Uzyta listW: maciez odleglosci znormalizowana najwieksza wartoscia wlasna
setwd("~/spatial_regression_project")

#SAR
SAR<-lagsarlm(data=dane, Turysci~PKBpc+number_of_pois+Urbanizacja, zero.policy = T, listw=norm)
summary(SAR)
dane$resSAR <- SAR$residuals
moran.test(dane$resSAR, listw = norm, zero.policy = T)

paleta2 <- colorRampPalette(c("green", "white", "blue"), bias = 1.15)
spplot(dane, zcol = "resSAR", colorkey = TRUE, col.regions = paleta2(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SAR")

#SEM
SEM<-GMerrorsar(data=dane, Turysci~PKBpc+number_of_pois+Urbanizacja, zero.policy = T, listw=norm)
summary(SEM)
dane$resSEM <- SEM$residuals
moran.test(dane$resSEM, listw = norm, zero.policy = T)

spplot(dane, zcol = "resSEM", colorkey = TRUE, col.regions = paleta2(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SEM")

#SLX
SLX <- lmSLX(dane$Turysci~dane$PKBpc+dane$number_of_pois+dane$Urbanizacja, zero.policy = T, listw=norm)
summary(SLX)
dane$resSLX<-SLX$residuals
lm.morantest(SLX, listw = norm, zero.policy = T)
lm.LMtests(SLX, listw = norm, test = "all", zero.policy = T)

paleta3 <- colorRampPalette(c("green", "white", "blue"), bias = 1.28)
spplot(dane, zcol = "resSLX", colorkey = TRUE, col.regions = paleta3(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SLX")

#Testy LM
modelLM<-lm(dane$Turysci~dane$PKBpc+dane$number_of_pois+dane$Urbanizacja)
lm.LMtests(modelLM, listw = norm, test = "all", zero.policy = T)

###Modele SARAR, SDM i SDEM
#SARAR
SARAR<-sacsarlm(data=dane, Turysci~PKBpc+number_of_pois+Urbanizacja, zero.policy = T, listw=norm)
summary(SARAR)
dane$resSARAR <- SARAR$residuals
moran.test(dane$resSARAR, listw = norm, zero.policy = T)

spplot(dane, zcol = "resSARAR", colorkey = TRUE, col.regions = paleta2(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SARAR")

#SDM
SDM<-lagsarlm(data=dane, Turysci~PKBpc+number_of_pois+Urbanizacja, zero.policy = T, listw=norm, type = "Durbin")
summary(SDM)
dane$resSDM <- SDM$residuals
moran.test(dane$resSDM, listw = norm, zero.policy = T)

paleta4 <- colorRampPalette(c("green", "white", "blue"), bias = 0.9)
spplot(dane, zcol = "resSDM", colorkey = TRUE, col.regions = paleta4(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SDM")

#SDEM
W <- listw2mat(norm)
X <- cbind(dane$PKBpc, dane$number_of_pois, dane$Urbanizacja)
WX <- W %*% X
lag.PKBpc <- WX [, 1]
lag.number_of_pois <- WX [, 2]
lag.Urbanizacja <- WX [, 3]
SDEM<-errorsarlm(data=dane, Turysci~PKBpc+number_of_pois+Urbanizacja+lag.PKBpc+lag.number_of_pois+lag.Urbanizacja, zero.policy = T, listw=norm)
summary(SDEM)
dane$resSDEM <- SDEM$residuals
moran.test(dane$resSDEM, listw = norm, zero.policy = T)

spplot(dane, zcol = "resSDEM", colorkey = TRUE, col.regions = paleta4(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Reszty przestrzenne modelu SDEM")

#Wybór najlepszego modelu
#Testy LM
modelLM<-lm(dane$Turysci~dane$PKBpc+dane$number_of_pois+dane$Urbanizacja)
lm.LMtests(modelLM, listw = norm, test = "all", zero.policy = T)

#Test LR
lL0 <- logLik(SEM) #SEM H0
lL1 <- logLik(SDM) #SDM H1
LR.Sarlm(lL0, lL1)

#Kryterium Akaike'a
AIC<-data.frame('Model'=c('SAR', 'SEM', 'SLX', 'SARAR', 'SDM', 'SDEM'),'AIC'=c(AIC(SAR),AIC(SEM),AIC(SLX),AIC(SARAR),AIC(SDM),AIC(SDEM)))
AIC

###WYbrany najlepszy model - SAR
###Interpretacja efektow przestrzennych
#Efekt bezposredni, posredni, calkowity
impacts.SAR <- impacts(SAR, listw = norm, zstats = TRUE, R = 200)
summary(impacts.SAR)

HPDinterval(impacts.SAR, prob=0.95)

#Wizualizacja efektów dla zmiennej PKBpc
#Wybrany region: Provence-Alpes-Cote d'Azur (FRL0)

#Macierz Sk :
SAR$coefficients
#Beta dla PKBpc - 0.676233990
Sk<-0.676233990*solve(diag(22)-as.numeric(SAR$rho)*nb2mat(norm$neighbours, zero.policy = T))
colnames(Sk)<-Francja$NUTS_ID
rownames(Sk)<-Francja$NUTS_ID

#Kolumna z efektami dla analizowanej prowincji
dane$FRL0<-round(Sk[,'FRL0'], 4)

#Mapa
paleta5 <- colorRampPalette(c("green", "blue"))
par(cex.main=1.8, cex.lab=1.8)
spplot(dane, zcol = "FRL0", colorkey = TRUE, col.regions = paleta5(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = list(label="Efekt mnozników przetrzennych dla zmiennej PKBpc (przyblizony do 4 miejsca po przecinku)", cex=1.8),
       sp.layout=list('sp.text', coordinates(dane), dane$FRL0, cex=1.8))

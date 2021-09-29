library(TransPhylo)
library(ape)
library(lubridate)
library(coda)
library(igraph)
library(reshape2)
library(tidyverse)
library(dplyr)
library(geosphere)
library(pointdensityP)

setwd("x")
dir.create("output")

# 1.Reconstruct transmission tree from BEAST time-scaled tree
# 2.Compute infection chain lengths (ICL) between samples from undirected transmission network
# 3.Compute movement pathlengths (MPL) between samples from animal movement network
# 4.Check correlation between ICL and MPL
# 5.Estimate farm-level R
# 6.Prepare network data for ERGMs

# 1.Reconstruct transmission tree from BEAST time-scaled tree
rm(list=ls())
set.seed(0)
t <- read.nexus("beastmcc.tree")
d <- read.csv("metadata.csv")
d$sampling_date <- as.Date(d$sampling_date)
m <- read.csv("movement.csv")
m$date <- as.Date(m$date)
m <- m %>% mutate(time = decimal_date(date))
l <- read.csv("locationall.csv")
tgen <- 14.5 #specify generation time of the pathogen(days)
dateT <- Inf #Ongoing outbreak (If the outbreak is ended, use the time at which observation of cases stopped)
rad <- 10 #farm proximity (km) that the disease is able to spread eg. via aerosol
  #1.1 Convert Beast Nexus tree to Newick 
write.tree(t, append = TRUE, file = "beastmcc.nwk")
ptree <- ptreeFromPhylo(read.tree("beastmcc.nwk"),dateLastSample=decimal_date(max(d$sampling_date)))
  #1.2 Compute gamma distribution and execute MCMC
mt <- tgen / 365.24 #convert to yearly time scale
s <- mt # Exponential (SIR) model
shape <- mt^2 / s^2
rate <- mt / s^2
scale <- 1/rate
c(shape, scale)
w.shape <- shape
w.scale <- scale
record<-inferTTree(ptree,mcmcIterations=50000,w.shape=w.shape,w.scale=w.scale,dateT=dateT)
record<-inferTTree(ptree,mcmcIterations=50000,w.shape=w.shape,w.scale=w.scale,dateT=dateT)
  #1.3 Diagnostics (trace & ESS)
pdf(file = "output/TTMCMCtrace.pdf") 
plotTraces(record)
dev.off() 
mcmc=convertToCoda(record)
sink(file = "output/TTMCMCess.txt")
effectiveSize(mcmc)
sink()
  #1.4 Plot transmission tree
lastIteration<-record[[length(record)]]
pdf(file = "output/transmissiontree.pdf") 
plotCTree(lastIteration$ctree)
dev.off() 
  #1.5 Get consensus transmission tree edges list
cons <- consTTree(record)
source <- cons[["ttree"]][,3]
target <- 1:nrow(cons[["ttree"]])
ttedge <- data.frame(source,target)

# 2.Compute infection chain lengths (ICL) between samples from undirected transmission network
r <- 1:length(cons[["nam"]])
n <- length(cons[["nam"]])
  #2.1 Compute shortest path length (SPL)
graph <- graph_from_data_frame(ttedge, directed = FALSE, vertices = NULL)
matrix <- shortest.paths(graph, v=V(graph), to=V(graph), mode = "out")
mat <- subset(melt(matrix), value!=0)
filtmat <- filter(mat, Var1<=n & Var1!=0 &  Var2<=n & Var2!=0 & value!=Inf)
filtmat2 <- filtmat[!duplicated(apply(filtmat[,c('Var1', 'Var2')],1,function(x) paste(sort(x),collapse=''))),]
old <- 1:length(cons[["nam"]])
new <- cons[["nam"]]
filtmat2$Var1[filtmat2$Var1 %in% old] <- new[match(filtmat2$Var1, old)]
filtmat2$Var2[filtmat2$Var2 %in% old] <- new[match(filtmat2$Var2, old)]
  #2.2 Assign source-target of edges based on sampling date
filtmat2$id <- 1:nrow(filtmat2) 
merge1 <- merge(x = filtmat2, y = d[ , c("name", "farm_name", "sampling_date")], by.x = "Var1", by.y = "name", all.x = TRUE, sort = FALSE)
merge1 <- merge1[order(merge1$id), ] %>% rename(v1_date = "sampling_date", v1_farm = "farm_name")
merge2 <- merge(x = merge1, y = d[ , c("name", "farm_name", "sampling_date")], by.x = "Var2", by.y = "name", all.x = TRUE, sort = FALSE)
merge2 <- merge2[order(merge2$id), ] %>% rename(v2_date = "sampling_date", v2_farm = "farm_name")
mut <- merge2 %>% mutate(ICL=value-1, sample.source = ifelse(v1_date <= v2_date, Var1, Var2), sample.target = ifelse(v2_date >= v1_date, Var2, Var1), origin = ifelse(v1_date <= v2_date, v1_farm, v2_farm), destination = ifelse(v2_date >= v1_date, v2_farm, v1_farm),Torigin = ifelse(v1_date <= v2_date, decimal_date(v1_date), decimal_date(v2_date)), Tdestination = ifelse(v2_date >= v1_date, decimal_date(v2_date), decimal_date(v1_date)), T1 = if_else(v1_date <= v2_date, v1_date, v2_date),T2 = if_else(v2_date >= v1_date, v2_date, v1_date))
mut.drop <- mut[,!(names(mut) %in% c("id","value", "Var1", "Var2", "v1_date", "v2_date"))]
col_order <- c("sample.source", "sample.target", "origin", "destination", "T1", "T2", "Torigin", "Tdestination","ICL")
ttedge.sample <- mut.drop[, col_order]

# 3.Compute movement pathlengths (MPL) between samples from animal movement network
ttedge.sample$origin <- as.character(ttedge.sample$origin)
ttedge.sample$destination <- as.character(ttedge.sample$destination)
m$tail <- as.character(m$tail)
m$head <- as.character(m$head)
x <- c(1:nrow(ttedge.sample))
mm.path <- c()
mm.path <- data.frame()
for (i in x){
  # Filter movement data between estimated Torigin & Tdestination 
  mfilt <- filter(m, time > ttedge.sample$Torigin[i] & time < ttedge.sample$Tdestination[i])
  # Graph
  graph <- graph_from_data_frame(mfilt, directed = TRUE, vertices = NULL) #convert filtered movement data to directed graph
  matrix <- shortest.paths(graph, v=V(graph), to=V(graph), mode = "out") #find all the shortest paths
  s <- subset(melt(matrix), value!= Inf)  #exclude Inf value
  distance <- filter(s, Var1 == ttedge.sample$origin[i] & Var2 == ttedge.sample$destination[i]) #select the match origin and destination
  dist <- c(distance$value[distance$value >= 0], NA)[1] #get a number of intermediates, if no path exists, turn NA
  mm.path.i <- data.frame(dist = dist)
  mm.path <- rbind(mm.path,mm.path.i)
  }
ttedge.mm <- cbind(ttedge.sample,mm.path$dist) %>% rename(MPL = "mm.path$dist")

# 4.Check correlation between ICL and MPL
lmd <- ttedge.mm[!is.na(ttedge.mm$MPL), ]
lm <- lm(ICL ~ MPL, data = lmd)
pdf(file = "output/normtestICL~MPL.pdf") 
par(mfrow=c(2,2)); plot(lm)
dev.off() 
sink(file = "output/lmICL~MP.txt")
cor.test(lmd$ICL, lmd$MPL, method = "pearson")
summary(lm)
sink()

#plot ICL against MPL data points
newx <- seq(min(lmd$MPL), max(lmd$MPL), length.out=100)
preds <- predict(lm, newdata = data.frame(MPL=newx), 
                 interval = 'confidence')
pdf(file = "output/ICL~MPL.pdf")
par(mfrow=c(1,1))
plot(lmd$ICL ~ lmd$MPL, ylab = "Infection chain length", xlab = "Movement pathlength", cex.lab = 1.25, cex.axis = 1.1)
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = '#deebf7', border = NA)
points(lmd$ICL ~ lmd$MPL,col="red", pch = 19)
abline(lm, col)
dev.off() 

# 5.Estimate farm-level R
  #5.1 Set a threshold for farm-to-farm direct transmission pair candidate selection using median ICL at MPL of 1
pdf(file = "output/bpICL~MPL.pdf")
par(mfrow=c(1,1))
bp <- boxplot(ICL ~ MPL, data = lmd, xlab = "Movement path length", ylab = "Infection chain length")
text(x = col(bp$stats) - .5, y = bp$stats, labels = bp$stats)
abline(h=median(lmd$ICL[lmd$MPL==1]), col="red")
dev.off() 
threshold <- median(lmd$ICL[lmd$MPL==1])
rcand <- filter(ttedge.mm, ICL < threshold & Tdestination-Torigin < 1) #ICL less than threshold and transmission event less than 1 year interval
  #5.2 Find the most possible transmission link for each destination and summarize 
translink <- rcand %>% group_by(sample.target) %>% filter(ICL == min(ICL)) %>% ungroup()
rval <- translink %>% count(sample.source) %>% rename(Rvalue = "n")
sink(file = "output/R-value.txt")
summary(rval)
sink()

# 6.Prepare network data for ERGMs

  #6.1 Dyads
# define transmission links
dtranslink <- translink
dtranslink$trans_link <- 1
dtranslink <- dtranslink[c("sample.source", "sample.target","trans_link")]
dmerge <- left_join(ttedge.mm,dtranslink, by = c("sample.source", "sample.target"))
dmerge[is.na(dmerge)] <- 0
# calculate time length (days)
dmerge$timelength_d <- dmerge$T2 - dmerge$T1
dmerge$timelength_d <- as.numeric(dmerge$timelength_d)
# calculate distance (km) between farms
dmerge$id <- 1:nrow(dmerge) 
coord1 <- merge(x = dmerge, y = d[ , c("name", "lat", "lon")], by.x = "sample.source", by.y = "name", all.x = TRUE, sort = FALSE)
coord1 <- coord1[order(coord1$id), ] %>% rename(source.lat = "lat", source.lon = "lon")
coord2 <- merge(x = coord1, y = d[ , c("name", "lat", "lon")], by.x = "sample.target", by.y = "name", all.x = TRUE, sort = FALSE)
coord2 <- coord2[order(coord2$id), ] %>% rename(target.lat = "lat", target.lon = "lon")
coord <- coord2[c("sample.source", "sample.target","trans_link","MPL","timelength_d","source.lon","source.lat","target.lon","target.lat")]
rownames(coord) <- NULL
#The shortest path between two points on an ellipsoid = the geodesic
dist <- coord %>% mutate(distance_km = distGeo(coord[,6:7], coord[,8:9])/1000)
dyad <- dist[c("sample.source", "sample.target","trans_link","MPL","timelength_d","distance_km")] %>% rename(v1="sample.source", v2="sample.target", pathlength="MPL")
write.csv(dyad, file = "output/dyad_ergm.csv")

  #6.2 Vertices
# season of outbreak from sampling date
dattr <- d
season <- c(
  "01" = "Winter", "02" = "Winter",
  "03" = "Spring", "04" = "Spring", "05" = "Spring",
  "06" = "Summer", "07" = "Summer", "08" = "Summer",
  "09" = "Fall", "10" = "Fall", "11" = "Fall",
  "12" = "Winter"
)
dattr$season <- season[format(dattr$sampling_date, "%m")]
# compute point density of a farm within rad km radius
pd <- pointdensity(df = l, ,lat_col = "lat",lon_col = "lon", grid_size = 1, radius = rad)
pdmerge <- left_join(l,pd, by = c("lat", "lon"))
dattr$id <- 1:nrow(dattr) 
attrmerge <- merge(x = dattr[ , c("id", "name", "farm_name", "sampling_date", "farm_type", "herd_size", "season")], y = pdmerge[ , c("farm_name", "count")], by.x = "farm_name", by.y = "farm_name", all.x = TRUE, sort = FALSE)
attrmerge <- attrmerge[order(attrmerge$id), ] %>% rename(p_density = "count")
# compute in and out degree of each farm during outbreak period (6 months)
md <- m
a <- c(unique(year(md$date)))
b <- a + 0.5
z <- data.frame(sort(c(a,b))) %>% rename(x="sort.c.a..b..") %>% mutate(y=x+0.5)
df <- c()
df <- data.frame()
for (i in 1:nrow(z)){
    # Filter movement data every half year
    mdfilt <- filter(md, time > z$x[i] & time < z$y[i])
    # Graph
    net.md <- graph_from_data_frame(mdfilt, directed = TRUE, vertices = NULL) %>% set_vertex_attr("year", value = z$x[i])
    df.i <- data.frame(farm_name = V(net.md)$name,
                      year = vertex_attr(net.md, "year"),
                      m_indeg = degree(net.md, mode = "in", normalized = TRUE),
                      m_outdeg = degree(net.md, mode = "out", normalized = TRUE))
    df <- rbind(df,df.i)
  }
attrhy <- attrmerge %>% mutate(decdate = decimal_date(sampling_date)) %>% mutate(year = ifelse((decdate %% 1)*10 < 5, year(sampling_date), year(sampling_date)+0.5))
attrhymerge <- left_join(attrhy,df, by = c("farm_name", "year"))
attrhymerge[is.na(attrhymerge)] <- 0
attr <- attrhymerge[c("name", "farm_name","farm_type","herd_size","p_density","season","m_indeg","m_outdeg")] 
write.csv(attr, file = "output/attr_ergm.csv")

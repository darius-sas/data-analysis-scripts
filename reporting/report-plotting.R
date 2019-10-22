library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(reshape2)

UDsmell <- "Unstable"
CDsmell <- "Cyclic"
HLsmell <- "Hublike"
GCsmell <- "God Comp."

cleanData <- function(df.smells, df.sizes){
  df.smells$smellTypeLabel <- plyr::revalue(df.smells$smellType, 
                             c("cyclicDep"=CDsmell, 
                               "hubLikeDep"=HLsmell, 
                               "unstableDep"=UDsmell, 
                               "godComponent"=GCsmell))
  df.smells$affectedComponentType <- plyr::revalue(df.smells$affectedComponentType, 
                                                   c("0"="package"))
  df.smells$versionDate <- as.Date(df.smells$versionDate)
  df.smells <- df.smells %>%
    mutate(typeOfSmell = paste(smellTypeLabel, affectedComponentType, sep=" on "))
  
  df.sizes$versionDate <- as.Date(df.sizes$versionDate)
  return(list(df.smells=df.smells, df.sizes=df.sizes))
}


plotNumOfSmellsOverTimeByType <- function(df.smells){
  df.count <- df.smells %>%
    group_by(versionPosition, versionDate, typeOfSmell) %>%
    tally()
  ggplot(df.count, aes(versionDate, n, group = typeOfSmell, color = typeOfSmell)) + 
    geom_line() + 
    scale_x_date(breaks = pretty) +
    labs(x="Time", y="Count", title = "Number of architectural smells over time", color = "Arch. Smells")
}

plotNumOfSmellsLastVersionByType <- function(df.smells){
  maxRow = which.max(df.smells$versionPosition);
  maxVersion = df.smells$versionPosition[maxRow]
  maxVersionDate = df.smells$versionDate[maxRow]
  df.count <- df.smells %>% filter(versionPosition == maxVersion) %>%
    group_by(versionPosition, typeOfSmell) %>%
    tally()
  ggplot(df.count, aes(typeOfSmell, n, fill = typeOfSmell, label = n)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(stat='identity', vjust=-1) +
    labs(x="Number of architectural smells by type", y = "Count", 
         title = paste("Number of smells on ", maxVersionDate, sep=""), 
         fill="Arch. Smells")
}

plotComponentCountPerVersion <- function(df.sizes){
  df.melt <- melt(df.sizes, id.vars = c("project", "version", "versionDate", "versionPosition"))
  df.melt$variable <- plyr::revalue(df.melt$variable, 
                                                   c("nPackages"="Packages", "nClasses"="Classes"))
  ggplot(df.melt, aes(x=versionDate, y=value, color=variable)) +
    geom_line() +
    scale_x_date(breaks = pretty) +
    labs(title=paste("Number of classes and packages per version"),
         x = "Versions",
         y = "Count",
         color="")
}

getTop5OldestSmellsByType <- function(df.smells){
  maxRow <- which.max(df.smells$versionPosition)
  lastVersion <- df.smells$version[maxRow]
  lastVersionDate <- df.smells$versionDate[maxRow]
  
  df.smells <- df.smells %>% filter(version == lastVersion)
  
  df.smells.curr <- df.smells %>% 
    group_by(typeOfSmell) %>%
    arrange(desc(age)) %>%
    top_n(5, age) %>% 
    select(uniqueSmellID)
  
  oldestSmells <- right_join(df.smells, df.smells.curr, by="uniqueSmellID")  %>%
    arrange(typeOfSmell.x, desc(age)) %>%
    select(age, affectedElements, typeOfSmell.x, size, numOfEdges, shape)
  
  return(oldestSmells)
}













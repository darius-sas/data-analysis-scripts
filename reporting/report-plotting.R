library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(tidyr)
library(survival)
library(survminer)

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
    theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
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

getTopOldestSmellsByType <- function(df.smells, n_top = 3){
  maxRow <- which.max(df.smells$versionPosition)
  lastVersion <- df.smells$version[maxRow]
  lastVersionDate <- df.smells$versionDate[maxRow]
  
  df.smells.cv <- df.smells %>% filter(version == lastVersion)
  
  df.top.smells <- df.smells.cv %>% 
    group_by(typeOfSmell) %>%
    filter(shape != "tiny" | is.na(shape)) %>%
    arrange(desc(age)) %>%
    top_n(n_top, age) %>% 
    select(uniqueSmellID, typeOfSmell)
  
  oldestSmells <- right_join(df.smells.cv, df.top.smells, by="uniqueSmellID")  %>%
    arrange(typeOfSmell.x, desc(age)) %>%
    transmute("Type of Smell" = typeOfSmell.x,
              Age = age,
              "Affected Elements" = stringr::str_replace_all(affectedElements, "\\$", "->"),
              Shape = shape,
              "Num. of Class/Pack." = size,
              "Num. of Depend." = numOfEdges)

  
  return(oldestSmells)
}


getComponentsAffectedByMostSmells <- function(df.comps, df.smells, n_top = 5){
  lastVersion <- df.smells$version[which.max(df.smells$versionPosition)]
  df.comps.lv <- df.comps %>% filter(version == lastVersion)
  
  df.comps.join <- right_join(df.smells, df.comps.lv, 
                            by=c("version"="version", "uniqueSmellID"="affectedByUniqueSmellId")) %>%
                    select(name, type, smellTypeLabel)
  
  df.comps.count <- df.comps.join %>%
    group_by(name, type, smellTypeLabel) %>%
    tally() %>%
    spread(smellTypeLabel, n, fill = 0) %>%
    mutate(Total = !!sym(UDsmell) + !!sym(HLsmell) + !!sym(CDsmell) + !!sym(GCsmell)) %>%
    group_by(type) %>%
    top_n(n_top, Total) %>%
    arrange(type, Total, name) %>%
    transmute(Name = name,Total = Total, 
              Cyclic = Cyclic, 
              "God Comp." = !!sym(GCsmell), 
              Unstable = Unstable,
              Hublike = Hublike)
  
  return(df.comps.count)
}


plotSurvivalProbabilities <- function(df.smells, strata = "typeOfSmell", base.size = 12){
    surv <- computeSurvivalAnalysis(df.smells, strata)
    p<-ggsurvplot(surv$model, data = surv$data,
                           surv.median.line = "v", short.panel.labs = T,
                           ggtheme = theme_grey(base_size = base.size),
                           title = "Survival probability of smells", facet.by = "project") +
      theme(legend.position = "right") +
      guides(color=guide_legend(ncol = 1, title = "Arch. Smells"))
   return(p)
}


computeSurvivalAnalysis<- function(df, strata = "smellType"){
  df.proj <- df %>% group_by(project) %>%
    summarise(firstVersion = version[which.min(versionPosition)], 
              lastVersion = version[which.max(versionPosition)],
              n.versions = max(versionPosition),
              n.smells = length(unique(uniqueSmellID)))
  df.smel <- df %>%
    group_by(uniqueSmellID, project) %>%
    summarise(lastVersion = max(versionPosition))
  df.smel <- left_join(df.smel, df.proj[,c("project","n.versions")], by="project")
  df.smel$presentInLastVersion <- df.smel$lastVersion == df.smel$n.versions
  
  df.dup <- df[!duplicated(df[, c("project", "uniqueSmellID")]), cbind(c("project", "uniqueSmellID", "age", "versionPosition"), strata)]
  df.surv <- left_join(df.smel, df.dup, by=c("project", "uniqueSmellID"))
  
  model.formula <- as.formula(paste("Surv(time = age, event = !presentInLastVersion) ~", strata, "+ project"))
  fitmodel = survfit(model.formula, data = df.surv)
  fitmodel$call$formula <- model.formula
  
  return(list(model=fitmodel, data = df.surv))
}










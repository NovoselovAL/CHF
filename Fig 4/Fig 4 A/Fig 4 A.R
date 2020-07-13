
###########################################################################################################################################################
#                                                                                                                                                         #
#                                                                Script for "figure 4 a" of paper                                                         # 
#  Microbiome changes in chronic heart failure with preserved ejection fraction patients correlate with fibrosis markers: description of Russian cohort   #
# R version 3.6.1 (2019-07-05)                                                                                                                            #
# Tutorial on drawing a correlation map using ggplot2                                                                                                     #
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)                                                                                           #
###########################################################################################################################################################

library(ggplot2)

# Set working directory
setwd("C:/R project/")
getwd()


abund_table<-read.csv("C:/path to/abundance table_4a.csv",row.names=1,header=TRUE, sep=";")
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
meta_table<- read.table("C:/path to/Data SCFA_TMAO",row.names=1, header=TRUE, sep="\t")

# You can use sel_env_label to specify the labels for the panel, in case of "figure 4 a" we have to use SCFA and TMAO

sel_env_label <- list(
  'Propionic.acid'="Propionic acid",
  'Hydroxide.pyruvate'="Hydroxide pyruvate",
  'Butyric.acid'="Butyric acid",
  'Isobutyric.acid'="Isobutyric acid",
  'Methylbutyric.acid'="Methylbutyric acid",
  'Isovalerianic.acid'="Isovalerianic acid",
  'Valeric.acid'="Valeric acid",
  'TMAO'="TMAO"
  
)

sel_env_label<-t(as.data.frame(sel_env_label))
sel_env_label<-as.data.frame(sel_env_label)
colnames(sel_env_label)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)


#######
#Apply normalisation (either use relative or log-relative transformation)

x<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))

x<-x[,order(colSums(x),decreasing=TRUE)]
#Extract list of top N Taxa
N<-15
taxa_list<-colnames(x)[1:N]

#
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-data.frame(meta_table)

#Get grouping information
grouping_info<-read.csv("C:/path to/metadata Groups.csv",row.names=1,header=TRUE, sep=";")

groups<-grouping_info[,1]

#You can use kendall, spearman, or pearson below:
method<-"kendall"


#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for signifant values
df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
df$Env<-factor(df$Env,as.character(df$Env))

#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label[as.character(value),"Trans"])
}

p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
p<-p+facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",labeller=Env_labeller)

ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df) +
  geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  geom_text(aes(label=Significance), color="black", size=10)+labs(y=NULL, x=NULL, fill=method) +
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",labeller=Env_labeller) 
\dontrun{
# Creat data frame for mcmc-like
library(rJava)      # Needed for tabulizer
library(tabulizer)  # Handy tool for PDF Scraping
#
#Data set A: Alstonville
crop_ler <- extract_tables(
    file   = system.file(paste0("external/Allen 1986.pdf"), package="contactsimulator"),
    pages = 11,
    method = "decide",
    output = "data.frame",
    guess = FALSE)

crop_ler_tbl <- crop_ler %>%
    pluck(1) %>%
    as_tibble()


diff_new<- function(x){
  c(30,diff(x))
}
# Show first 6 rows
crop_ler_tbl %>% head() %>% knitr::kable()
# Extract the part of the part of data
crop_data<- crop_ler_tbl%>%slice(3:n())%>%as.data.frame()
crop_data[,1]=gsub(" ", "", crop_data[,1])
crop_data[,1]=gsub("\\s+|[-]", "", crop_data[,1])
crop_data[,1]=gsub("^(.{6})", "\\1:", crop_data[,1])
crop_data_update<- crop_data%>% separate(Month.P.1.ant,c("Month","Age"), sep = ":")%>%dplyr::select(-c(X,X.1))
crop_data_update<- as.data.frame(apply(crop_data_update,2,function(x){
   x<- gsub("\\s+","",x)
   return(x)
}))

crop_data_update[,2:ncol(crop_data_update)]<- as.data.frame(apply(crop_data_update[,2:ncol(crop_data_update)],2,function(x)as.numeric(gsub('\\s+|[-]', '',x))))
crop_data_update$LER.9[25]<- 0.1193
crop_data_update$LER.5[27]<- 0.0733
crop_data_update$LER.15[27]<- 0.0633
# Reshape from large to long
names(crop_data_update)[4:21]=c("Tmin",paste0("Cultiv.",1:17))
crop_data_update_reshape<- crop_data_update%>%gather(Cultivars,LER,Cultiv.1:Cultiv.17)
crop_data_update_ler<- crop_data_update_reshape%>%group_by(Cultivars)%>%mutate(Age_=diff_new(Age))
crop_data_update_ler$leaf<- crop_data_update_ler$LER*crop_data_update_ler$Age_
crop_data_update_ler$Date<- sapply(crop_data_update_ler$Month,function(x){
  x<- gsub("\\s+","",x)
  x<- gsub("[,]", "19", x)
  return(paste0("15",x))
})
crop_data_update_ler$Date<-  as.Date(crop_data_update_ler$Date,"%d%b%Y")
crop_data_update_ler$day<- as.numeric(difftime(as.Date(crop_data_update_ler$Date),as.Date(paste0(year(crop_data_update_ler$Date),"-01-01")),units = "days"))
L<- lapply(1:nrow(crop_data_update_ler), function(ii){
   xx<- lapply(1:500,function(ll){
     sp<- sample(90000:100000,1)
     x1<- rBTFinv4(1,crop_data_update_ler$day[ii],samp[sp,1],samp[sp,2],ceiling(crop_data_update_ler$leaf[ii]))
     return(x1)
   })
   xx<- unlist(xx)
   return(c(quantile(xx,c(0.025,0.5,0.975)),mean(xx)))
 })
LL<- do.call(rbind,L)
LL<- LL%>%as.data.frame()
names(LL)<- c("min","median","max","mean")
crop_data_all<- cbind(as.data.frame(crop_data_update_ler),LL)
crop_data_all$Date<- gsub("197-","1978-",crop_data_all$Date)
crop_data_all%>%ggplot(aes(Date,Age_sum)) + geom_point(aes(color="red")) +
  theme_bw() +
  scale_x_discrete(breaks=crop_data_all$Date[seq(1,length(crop_data_all$Date),3)],labels=crop_data_all$Date[seq(1,length(crop_data_all$Date),3)])+
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) + facet_wrap(.~Cultivars,nrow=6)+
  xlab("Date of observation") + ylab("Time taken for leaves to grow (days)")


crop_data_all%>%ggplot(aes(Age,Age_)) + geom_line(aes(color="red")) + geom_line(aes(Age,median,color="black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.85,0.05), text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  scale_color_manual(values=c('black','red'), labels=c("Median","Obsevation"))+
  guides(color=guide_legend(" ")) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha=0.2) + facet_wrap(.~Cultivars_f,nrow=6)+
  xlab("Time after planting (days)") + ylab("Time taken for leaves to grow (days)")

# Data set C


crop_ler_b <- extract_tables(
   file   = system.file(paste0("external/Allen 1986.pdf"), package="contactsimulator"),
    pages = c(12,12),
    area = list(c(72.57776,158.17891,555.27855,419.42923),
                c(72.57776,421.47025,456.28917,675.57700)),
    method = "decide",
    output = "data.frame",
    guess = FALSE)

Lis<- vector(mode = "list", length = 2)
for(i in 1:2){
  crop_ler_tbl_b <- crop_ler_b %>%
    pluck(i) %>%
    as_tibble()
  crop_data_b <- crop_ler_tbl_b%>%slice(5:n())%>%dplyr::select(-c("Sub..Variety.","X"))%>%as.data.frame()%>%
  na.omit()
  crop_data_b<- crop_data_b[!apply(crop_data_b == "", 1, all),]
  Lis[[i]]<- crop_data_b
}
crop_data_b_merge<- do.call(rbind,Lis)
# Show first 6 rows
locations<- "NSW"
locations<- paste(locations,1:5,sep = "-")
crop_data_b_merge$Subset<- rep(locations,times=c(17,3,15,14,14))
crop_data_b_update<- as.data.frame(apply(crop_data_b_merge,2,function(x){
   x<- gsub("\\s+","",x)
   return(x)
}))
names(crop_data_b_update)<- names(crop_data_b_merge)
crop_data_b_update$Date<- gsub("[?]","p",crop_data_b_update$Date)
crop_data_b_update$Date<- gsub("[-]","",crop_data_b_update$Date)
crop_data_b_update$Date<- gsub("t1","M",crop_data_b_update$Date)
vars<- c("Plant","Leaf.emergence")
crop_data_b_clean<- crop_data_b_update%>%mutate_each_(funs(as.character),vars)%>%mutate_each_(funs(as.numeric),vars)
crop_data_b_clean<- crop_data_b_clean%>%group_by(Subset)%>%mutate(Age=diff_new(Plant))%>%as.data.frame()
crop_data_b_clean$leaf<- crop_data_b_clean$Leaf.emergence*crop_data_b_clean$Age
crop_data_b_clean$Date<- sapply(crop_data_b_clean$Date,function(x){
  x<- gsub("[,]", "19", x)
  x<- gsub("\\.","",x)
  return(paste0("15",x))
})
 crop_data_b_clean[32,1]<- "15Dec1980"
 crop_data_b_clean[45,1]<- "15Dec1979"
 crop_data_b_clean$Month<-  as.Date(crop_data_b_clean$Date,"%d%b%Y")
 crop_data_b_clean$day<- as.numeric(difftime(as.Date(crop_data_b_clean$Month),as.Date(paste0(year(crop_data_b_clean$Month),"-01-01")),units = "days"))

 L<- lapply(1:nrow(crop_data_b_clean), function(ii){
   xx<- lapply(1:500,function(ll){
     sp<- sample(90000:100000,1)
     x1<- rBTFinv4(1,crop_data_b_clean$day[ii],samp[sp,1],samp[sp,2],ceiling(crop_data_b_clean$leaf[ii]))
     return(x1)
   })
   xx<- unlist(xx)
   return(c(quantile(xx,c(0.025,0.5,0.975)),mean(xx)))
 })

LL<- do.call(rbind,L)
LL<- LL%>%as.data.frame()
names(LL)<- c("min","median","max","mean")
crop_data_all<- cbind(crop_data_b_clean,LL)

crop_data_all%>%ggplot(aes(Plant,Age)) + geom_line(aes(color="red")) + geom_line(aes(Plant,median,color="black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.85,0.05), text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  scale_color_manual(values=c('black','red'), labels=c("Median","Obsevation"))+
  guides(color=guide_legend(" ")) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha=0.2) + facet_wrap(.~Subset,nrow=3, scales = "free")+
  xlab("Time after planting (days)") + ylab("Time taken for leaves to grow (days)")

crop_data_all%>%ggplot(aes(Plant,Age)) + geom_line(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_ribbon(aes(ymin=min, ymax=max),alpha=0.1) + facet_wrap(.~Subset,nrow=3, scales = "free")+
  xlab("Time after planting (days)") + ylab("Time taken for leaves to grow (days)")

NSW1<- crop_data_all%>%filter(Subset=="NSW-1")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

NSW2<- crop_data_all%>%filter(Subset=="NSW-2")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

NSW3<- crop_data_all%>%filter(Subset=="NSW-3")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

NSW4<- crop_data_all%>%filter(Subset=="NSW-4")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

NSW5<- crop_data_all%>%filter(Subset=="NSW-5")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

plot_grid(NSW1, NSW2, NSW3, NSW4, NSW5, labels = "AUTO")

#------------------------------------------------------------------
# Data set B
#-----------------------------------------------------------------

crop_ler_c <- extract_tables(
    file   = system.file(paste0("external/Allen 1986.pdf"), package="contactsimulator"),
    pages = c(13,13),
    area = list(c(78.32072, 151.04116, 559.18646, 393.52901 ),
                c(82.43068, 393.52901, 546.85657, 636.01687)),
    method = "decide",
    output = "data.frame",
    guess = FALSE)

Lis<- vector(mode = "list", length = 2)
for(i in 1:2){
  crop_ler_tbl_c <- crop_ler_c %>%
    pluck(i) %>%
    as_tibble()
  crop_data_c <- crop_ler_tbl_c%>%slice(5:n())%>%dplyr::select(-c("Sub..Variety.","X"))%>%as.data.frame()%>%
  na.omit()
  crop_data_c<- crop_data_c[!apply(crop_data_c == "", 1, all),]
  Lis[[i]]<- crop_data_c
}
names(Lis[[2]])[3]<- "Tmax"
crop_data_c_merge<- do.call(rbind,Lis)
# Show first 6 rows
locations<- "RSA"
locations<- paste(locations,1:4,sep = "-")
crop_data_c_merge$Subset<- rep(locations,times=c(22,12,34,12))
crop_data_c_update<- as.data.frame(apply(crop_data_c_merge,2,function(x){
   x<- gsub("\\s+","",x)
   return(x)
}))
names(crop_data_c_update)<- names(crop_data_c_merge)
crop_data_c_update[,2:5]<- as.data.frame(apply(crop_data_c_update[,2:5],2,function(x){
   x<- gsub("D","0",x)
   return(x)
}))

vars<- c("Plant","Leaf.emergence")
crop_data_c_clean<- crop_data_c_update%>%mutate_each_(funs(as.character),vars)%>%mutate_each_(funs(as.numeric),vars)
crop_data_c_clean<- crop_data_c_clean%>%group_by(Subset)%>%mutate(Age=diff_new(Plant))%>%as.data.frame()
crop_data_c_clean$leaf<- crop_data_c_clean$Leaf.emergence*crop_data_c_clean$Age
crop_data_c_clean$Date<- sapply(crop_data_c_clean$Date,function(x){
  x<- gsub("\\.","",x)
  x<- gsub("[,]", "19", x)
  x<- gsub("\\.","",x)
  return(paste0("15",x))
})
crop_data_c_clean[33,1]<- "15May1983"
crop_data_c_clean[79,1]<- "15May1983"
crop_data_c_clean$Month<-  as.Date(crop_data_c_clean$Date,"%d%b%Y")
 crop_data_c_clean$day<- as.numeric(difftime(as.Date(crop_data_c_clean$Month),as.Date(paste0(year(crop_data_c_clean$Month),"-01-01")),units = "days"))

 L<- lapply(1:nrow(crop_data_c_clean), function(ii){
   xx<- lapply(1:500,function(ll){
     sp<- sample(90000:100000,1)
     x1<- rBTFinv4(1,crop_data_c_clean$day[ii],samp[sp,1],samp[sp,2],ceiling(crop_data_c_clean$leaf[ii]))
     return(x1)
   })
   xx<- unlist(xx)
   return(c(quantile(xx,c(0.025,0.5,0.975)),mean(xx)))
 })


LL<- do.call(rbind,L)
LL<- LL%>%as.data.frame()
names(LL)<- c("min","median","max","mean")
crop_data_all<- cbind(crop_data_c_clean,LL)

crop_data_all%>%ggplot(aes(Plant,Age)) + geom_line(aes(color="red")) + geom_line(aes(Plant,median,color="black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.85,0.35), text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  scale_color_manual(values=c('black','red'), labels=c("Median","Obsevation"))+
  guides(color=guide_legend(" ")) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha=0.2) + facet_wrap(.~Subset,nrow=3, scales = "free")+
  xlab("Time after planting (days)") + ylab("Time taken for leaves to grow (days)")


RSA1<- crop_data_all%>%filter(Subset=="RSA-1")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)

#------------------------------------------------------------------
# Data set B
#-----------------------------------------------------------------

crop_ler_d <- extract_tables(
    file   = system.file(paste0("external/Allen 1986.pdf"), package="contactsimulator"),
    pages = c(15,15,15),
    area = list(c(89.02734, 150.39971, 506.46328, 355.02517),
                c(89.02734, 361.16393, 504.41702, 562.72001),
                c(85.95796, 566.81252, 449.16815, 771.43798 )),
    method = "decide",
    output = "data.frame",
    guess = FALSE)

Lis<- vector(mode = "list", length = 3)
for(i in 1:3){
  crop_ler_tbl_d <- crop_ler_d %>%
    pluck(i) %>%
    as_tibble()
  if(i==1){
    crop_data_d <- crop_ler_tbl_d%>%slice(5:n())%>%as.data.frame()%>%
  na.omit()
  }
  else{
    crop_data_d <- crop_ler_tbl_d%>%slice(4:n())%>%as.data.frame()%>%
  na.omit()
  }

  crop_data_d<- crop_data_d[!apply(crop_data_d == "", 1, all),]
  Lis[[i]]<- crop_data_d
}

for(i in 1:3){

  if(i==1 | i==3){
    Lis[[i]]$Date<- paste(Lis[[i]]$Sub..Date,Lis[[i]]$X,sep=" ")
    Lis[[i]]<- Lis[[i]]%>%mutate(Date=gsub("\\bQL\\w+","",Date))%>%
      mutate(Date=gsub("\\bRI\\w+","",Date))%>%
      mutate(Date=gsub("\\bcon\\w+","",Date))%>%
      dplyr::select(-starts_with("X"))%>%
      separate(Plant.Tmax,into=c("Plant","Tmax"),sep = 3)%>%dplyr::select(-starts_with("Sub"))
  }
  else{
    Lis[[i]]<- Lis[[i]]%>%rename(Date="Sub..Date")%>%mutate(Date=gsub("\\bQL\\w+","",Date))%>%
      mutate(Date=gsub("\\bRI\\w+","",Date))%>%
      mutate(Date=gsub("\\bcon\\w+","",Date))
  }

}
crop_data_c_merge<- do.call(rbind,Lis)
crop_data_c_merge<- as.data.frame(gsub(" ", "", as.matrix(crop_data_c_merge)))
crop_data_c_merge<- as.data.frame(gsub("l", "1", as.matrix(crop_data_c_merge)))
crop_data_c_merge<- as.data.frame(gsub("o", "0", as.matrix(crop_data_c_merge)))
# Show first 6 rows
locations<- c("QLD","RIC")
locations<- c(paste(locations[1],1:7,sep = "-"),paste(locations[2],1:4,sep = "-"))
crop_data_c_merge$Subset<- rep(locations,times=c(7,8,5,5,4,8,4,12,12,13,7))
crop_data_c_update<- as.data.frame(apply(crop_data_c_merge,2,function(x){
   x<- gsub("\\s+","",x)
   return(x)
}))
names(crop_data_c_update)<- names(crop_data_c_merge)

vars<- c("Plant","Leaf.emergence")
crop_data_d_clean<- crop_data_c_update%>%mutate_each_(funs(as.character),vars)%>%mutate_each_(funs(as.numeric),vars)
crop_data_d_clean<- crop_data_d_clean%>%group_by(Subset)%>%mutate(Age=diff_new(Plant))%>%as.data.frame()
crop_data_d_clean<- crop_data_d_clean%>%mutate(Age=ifelse(Age>31,30,Age))
crop_data_d_clean$leaf<- crop_data_d_clean$Leaf.emergence*crop_data_c_clean$Age
crop_data_d_clean$Date<- sapply(crop_data_d_clean$Date,function(x){
  x<- gsub("\\.","",x)
  x<- gsub("[,]", "19", x)
  x<- gsub("\\.","",x)
  return(paste0("15",x))
})

crop_data_d_clean$Date<- gsub("N0","No",crop_data_d_clean$Date)
crop_data_d_clean$Date<- gsub("R!C1","",crop_data_d_clean$Date)
crop_data_d_clean$Date<- gsub("Ju1","Jul",crop_data_d_clean$Date)
crop_data_d_clean$Date<- gsub("Oec","Dec",crop_data_d_clean$Date)

crop_data_d_clean$Month<-  as.Date(crop_data_d_clean$Date,"%d%b%Y")
 crop_data_d_clean$day<- as.numeric(difftime(as.Date(crop_data_d_clean$Month),as.Date(paste0(year(crop_data_c_clean$Month),"-01-01")),units = "days"))

 L<- lapply(1:nrow(crop_data_d_clean), function(ii){
   xx<- lapply(1:500,function(ll){
     sp<- sample(90000:100000,1)
     x1<- rBTFinv4(1,crop_data_d_clean$day[ii],samp[sp,1],samp[sp,2],ceiling(crop_data_d_clean$leaf[ii]))
     return(x1)
   })
   xx<- unlist(xx)
   return(c(quantile(xx,c(0.025,0.5,0.975)),mean(xx)))
 })



LL<- do.call(rbind,L)
LL<- LL%>%as.data.frame()
names(LL)<- c("min","median","max","mean")
crop_data_all<- cbind(crop_data_d_clean,LL)

crop_data_all%>%filter(grepl("QLD",Subset))%>%ggplot(aes(Plant,Age)) + geom_line(aes(color="red")) + geom_line(aes(Plant,median,color="black"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.50,0.20), text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  scale_color_manual(values=c('black','red'), labels=c("Median","Obsevation"))+
  guides(color=guide_legend(" ")) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha=0.2) + facet_wrap(.~Subset,nrow=3, scales = "free")+
  xlab("Time after planting (days)") + ylab("Time taken for leaves to grow (days)")


RSA1<- crop_data_all%>%filter(Subset=="RSA-1")%>%ggplot(aes(as.factor(Month),Age)) + geom_point(aes(color="red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.5,
                 position=position_dodge(.9)) +
  xlab("Date of observation") + ylab("Length of period of leaf growth")+
  facet_grid(.~Subset)


#--------------------------------------------------------------------------------------------------------------------
#                                      Simulation of the incubation period
#--------------------------------------------------------------------------------------------------------------------

Inc_pe=array(0,c(500,365))
for(i in 1:365){
  xx<- lapply(1:500,function(ll){
     sp<- sample(90000:100000,1)
     x1<- rBTFinv3(1,i,samp[sp,1],samp[sp,2],4)
     return(x1)
   })
  Inc_pe[,i]=unlist(xx)# From distribution
}

                   # Using Rob Allen
Ip<- sapply(1:365, function(t1) uniroot(fu,c(-1000,1000),t1=t1,l=2)$root) - 1:365
Lp<- sapply(1:365, function(t1) uniroot(fu,c(-1000,1000),t1=t1,l=3.7)$root) - 1:365
                   # plots

Inc_ci=t(apply(Inc_pe,2,quantile,c(0.025,0.05,.5,0.95,.975)))%>%as.data.frame()
Inc_ci$mean<- apply(Inc_pe,2,mean)
names(Inc_ci)[1:5]<- c("fquan","squan","med","tquan","laquan")
Inc_ci$rob_ip<- Ip
Inc_ci$rob_lp<- Lp

ggplot(Inc_ci,aes(1:365,mean))+  geom_line(aes(1:365,mean,colour="Posterior Mean"))+
  geom_line(aes(1:365,med, color="Posterior Median"))+
  geom_line(aes(1:365,rob_ip, color="R Allen incubation period"))+
  geom_ribbon(aes(ymin=fquan,ymax=laquan),alpha=0.1) +
  geom_ribbon(aes(ymin=squan,ymax=tquan),alpha=0.3)+ ylab("Days after inoculation") + xlab("Days since 1 January")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.87,.93),
        legend.title =element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(name='', values=c('Posterior Median'='red', 'Posterior Mean'='black', "R Allen incubation period"="blue"))+
  ggtitle("95% confidence interval of the time it takes for 2 leaves to emerge")


ggplot(Inc_ci,aes(1:365,mean))+  geom_line(aes(1:365,mean,colour="Posterior Mean"))+
  geom_line(aes(1:365,med, color="Posterior Median"))+
  geom_line(aes(1:365,rob_lp, color="R Allen latent period"))+
  geom_ribbon(aes(ymin=fquan,ymax=laquan),alpha=0.1) +
  geom_ribbon(aes(ymin=squan,ymax=tquan),alpha=0.3)+ ylab("Days after inoculation") + xlab("Days since 1 January")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.87,.93),
        legend.title =element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(name='', values=c('Posterior Median'='red', 'Posterior Mean'='black', "R Allen latent period"="green"))+
  ggtitle("95% confidence interval of the time it takes for 4 leaves to emerge")

###############################################################################################################
#                                             MCMC and simualtion
##############################################################################################################

data_tran<- function(data_f){
  if("Leaf.emergence"%in%names(data_f)){
    names(data_f)[which(names(data_f)=="Leaf.emergence")]<- "LER"
  }

  if("Subset"%in%names(data_f)){
    names(data_f)[which(names(data_f)=="Subset")]<- "Cultivars"
  }

  if("Plant"%in%names(data_f)){
    names(data_f)[which(names(data_f)=="Age")]<- "Age_"
    names(data_f)[which(names(data_f)=="Plant")]<- "Age"
  }
  return(data_f)
}

data_f<- data_tran(crop_data_c_clean)
df<- data_f%>%mutate(t_1=day - Age_,t_2=day, obs=ceiling(leaf))%>%as.data.frame()



xx<- mcmc_leaf(df[,c(11,12,13)],0.062,0.903,samp = 100000)

y<- lapply(1:nrow(xx[90000:100000,]),function(ii){
  #ss<- Sum_Leaf_emergence_rate(df$t_1,df$t_2,xx[ii,1],xx[ii,2])
  ss<- xx[ii,1]*(1+ xx[ii,2]*cos(2*pi/365*(df$day)))
  return(ss)
})

yy<- do.call(cbind,y)
yy1<- t(apply(yy,1,quantile,c(0.025,0.5,0.975)))
df_sim<- data.frame(lower_q=yy1[,1],median_q=yy1[,2],upper_q=yy1[,3],obs=df$LER)

df_<- cbind(as.data.frame(data_f),df_sim)

df_%>%filter(Subset=="NSW-1")%>%ggplot(aes(Age,median_q)) +geom_point() + geom_line() + geom_point(aes(Age,leaf),col="red") + geom_line(aes(Age,leaf),col="red") +
  geom_ribbon(aes(ymin=lower_q,ymax=upper_q),alpha=0.1) + geom_point(aes(Age,ss1),col="blue") + geom_line(aes(Age,ss1),col="blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 20, family = "Tahoma",face = "bold"))+
  xlab("Date since planting") + ylab("Number of leaf")

ss1<- 0.062*(df1$t_2-df1$t_1 + 0.903/(2*pi/365)*(sin(2*pi/365*(df1$t_2-15))-sin(2*pi/365*(df1$t_1-15))))


LL1<- list(do.call(rbind,list(L1[[2]],L1[[3]],L1[[4]])))
LL1[[2]]<- as.data.frame(L1[[1]])
names(LL1[[1]])<- c("Month","Age","Tmax","Tmin","LER","Cultivars","Age_","leaf","Date","day")
DF<- do.call(rbind,LL1)
DF1<- DF%>%mutate(Latitude=ifelse(grepl("Cultiv",Cultivars),-28.85,ifelse(grepl("RSA",Cultivars),-25.12,ifelse(grepl("RIC",Cultivars),5.6289,-17.5226))),
                  Longitude=ifelse(grepl("Cultiv",Cultivars),153.4496,ifelse(grepl("RSA",Cultivars),31.08,ifelse(grepl("RIC",Cultivars),-4.0823,146.0285))))%>%
  mutate(Locations=ifelse(grepl("Cultiv",Cultivars),"Alstonville",ifelse(grepl("RSA",Cultivars),"Burgershall",ifelse(grepl("RIC",Cultivars),"Azaguie","Alstonville"))))

#write.table(unique(DF1[,c("Date","Tmax","Tmin","Latitude","Longitude","Locations")]),"/media/sf_D_DRIVE/BBTV_PROJECT/Exported-data-after-cleaning/LER.csv",row.names = F)

# LER
df_$culti_<- factor(df_$Cultivars, levels=unique(df_$Cultivars))
df_%>%filter(grepl(paste(paste0("Cultiv.",13:17,"$"),collapse = "|"),Cultivars))%>%
  ggplot(aes(Age,median_q)) + geom_line() + geom_line(aes(Age,LER),col="red") +
  geom_ribbon(aes(ymin=lower_q,ymax=upper_q),alpha=0.1) +
  theme_bw() + facet_wrap(~culti_, ncol = 3, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none", text = element_text(size = 10, family = "Tahoma",face = "bold"))+
  xlab("Date since planting") + ylab("Number of leaf")


df_%>%filter(!grepl(paste(paste0("Cultiv.",1:17,"$"),collapse = "|"),Cultivars))%>%
  ggplot(aes(median_q,LER,colour=as.factor(Cultivars))) + geom_point()





}

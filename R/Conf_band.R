#' Visualisation of the crdible band
#'
#' Posterior credible band.
#'
#'\code{Conf_Band} Provide a 95% credible band.
#'
#' @param OBS Trajectory. See \code{\link{Trajectory_all}}
#' @param samp Tinning axis label
#' @param labl Axis label
#' @param t0bs limit of the observed time. This is used if user one want both observed and predicted band on the same graph,
#' @return A 95% credible band for the obsrvation. By default it considers the removals.
#'
#' @examples
#' data(rast2)
#' data(Dat_obs)
#' data(Param)
#' #Conf_Band(rast2,Param,Dat_obs=Dat_obs)
#' @export
Conf_Band<-function(OBS,samp=NULL,labl=NULL,ylbl=NULL, tobs=NULL, limt=NULL){
  #Extrct infos from the raster
  df3<- t(apply(OBS[,3:ncol(OBS)],1,quantile,c(0.025,.5,.975)))
  colnames(df3)<- c("min","median","max")
  df3<-as.data.frame(df3,row.names = FALSE)
  df3$times1<- OBS$times
  df3$obs<-   OBS$obs

  lab=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

  if(is.null(ylbl)){
    ylbl="Incidence"
  }

  if(is.null(limt)){
    limt <- nrow(df3)
  }

  if(is.null(tobs)){
    if(is.null(samp) & is.null(labl)){
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times[1:limt],median[1:limt]), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min[1:limt],ymax=max[1:limt]),alpha=0.3) + xlab("Date") + ylab(ylbl) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_colour_discrete( breaks=c("0", "1"),labels=c("Plan", "Back"))+
      ggplot2::theme(text = element_text(size = 20, family = "Tahoma",face = "bold"))

  }
  else{
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times1[1:limt],median[1:limt]), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min[1:limt],ymax=max[1:limt]),alpha=0.3) + ggplot2::scale_x_discrete(limits=times1[samp],labels=labl[samp])+
      xlab("Date") + ylab(ylbl) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_colour_discrete( breaks=c("0", "1"),labels=c("Plan", "Back"))+
      ggplot2::theme(text = element_text(size = 20, family = "Tahoma",face = "bold"))


  }
  }

  else{
    if(is.null(samp) & is.null(labl)){
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times[1:limt],median[1:limt]), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min[1:limt],ymax=max[1:limt], fill=times1<=tobs),alpha=0.3) + xlab("Date") + ylab(ylbl) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_colour_discrete( breaks=c("0", "1"),labels=c("Plan", "Back"))+
      ggplot2::theme(text = element_text(size = 20, family = "Tahoma",face = "bold"))+
        scale_fill_manual(values=c("pink", "gray"), name="fill")

  }
  else{
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times1[1:limt],median[1:limt]), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min[1:limt],ymax=max[1:limt]),alpha=0.3) +
      ggplot2::geom_ribbon(aes(ymin=min[1:limt],ymax=max[1:limt], fill=times1<=tobs),alpha=0.3) +
      ggplot2::scale_x_discrete(limits=times1[samp],labels=labl[samp])+
      xlab("Date") + ylab(ylbl) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_colour_discrete( breaks=c("0", "1"),labels=c("Plan", "Back"))+
      ggplot2::theme(text = element_text(size = 20, family = "Tahoma",face = "bold")) +
       scale_fill_manual(values=c("pink", "gray"), name="fill", labels = c("Out of sample", "In sample"))+
      theme(legend.title = element_blank(), legend.position="top")
      # scale_fill_discrete( labels = c("Out of sample", "In sample"))

  }
  }






}


#' @export
Conf_Band_year<-function(OBS,df1=NULL, ylabl=NULL){
  #Extrct infos from the raster
  lab=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  latentt<- c("beta[1]==0","omega==2*pi")
  trans1=c("Cyclic ","Gamma ","Exponential ")
  if(is.null(df1)){
    df<- OBS
  }
  else{
    df3<- do.call(rbind,OBS)
    colnames(df3)<- c("min","median","max")
    df<-as.data.frame(df3,row.names = FALSE)
    df$times<- df1$times
    df$obs<-   df1$obs
    df$year<- df1$year
  }

  if(is.null(ylabl)){
    ylabl<- "Number removed"
  }

  df<- df%>%arrange(year)

  year = as.character(unique(df$year))
  #dat<- data.frame(times=unique(df$times),month=lab)
  #df<- dplyr::left_join(df,dat,by="times")

  ggplot(df,aes(times,obs,group=year))+
    geom_point()+ geom_line(aes(times,median), color="red")+ geom_line()+
    geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3) +
    facet_wrap(~year, ncol=3, scales = "free") +
    scale_x_continuous(breaks=c(30,58,89,119,150,180,211,242,272,303,333,364),labels=unique(df$month))+
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 10, family = "Tahoma",face = "bold")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Months") + ylab(ylabl)






}

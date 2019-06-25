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
Conf_Band<-function(OBS,samp=NULL,labl=NULL,ylbl=NULL, tobs=NULL){
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

  if(is.null(tobs)){
    if(is.null(samp) & is.null(labl)){
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times,median), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3) + xlab("Date") + ylab(ylbl) +
      theme(text = element_text(size = 15))

  }
  else{
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times1,median), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3) + ggplot2::scale_x_discrete(limits=times1[samp],labels=labl[samp])+
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+ xlab("Date") + ylab(ylbl) +
      theme(text = element_text(size = 15))

  }
  }

  else{
    if(is.null(samp) & is.null(labl)){
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times,median), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min,ymax=max, fill=times1<=tobs),alpha=0.3) + xlab("Date") + ylab(ylbl) +
      theme(text = element_text(size = 15))+
        scale_fill_manual(values=c("pink", "gray"), name="fill")

  }
  else{
    ggplot2::ggplot(df3,aes(times1,obs))+
      ggplot2::geom_point()+ ggplot2::geom_line(aes(times1,median), color="red")+ ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3) +
      ggplot2::geom_ribbon(aes(ymin=min,ymax=max, fill=times1<=tobs),alpha=0.3) +
      ggplot2::scale_x_discrete(limits=times1[samp],labels=labl[samp])+
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+ xlab("Date") + ylab(ylbl) +
      theme(text = element_text(size = 15)) +
       scale_fill_manual(values=c("pink", "gray"), name="fill", labels = c("Out of sample", "In sample"))+
      theme(legend.title = element_blank(), legend.position="top")
      # scale_fill_discrete( labels = c("Out of sample", "In sample"))

  }
  }






}


#' @export
Conf_Band_year<-function(OBS,df1){
  #Extrct infos from the raster
  df3<- do.call(rbind,OBS)
  colnames(df3)<- c("min","median","max")
  df<-as.data.frame(df3,row.names = FALSE)
  df$times<- df1$times
  df$obs<-   df1$obs
  df$year<- df1$year

  lab=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  latentt<- c("beta[1]==0","omega==2*pi")
  trans1=c("Cyclic ","Gamma ","Exponential ")
  year = c("2011","2012","2013","2014","2015","2016","2017")

  ggplot(df,aes(times,obs,group=year))+
    geom_point()+ geom_line(aes(times,median), color="red")+ geom_line()+
    geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3) +
    facet_wrap(~year, ncol=3) + scale_x_discrete(limits=times,labels=lab)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Months") + ylab("Number removed")






}

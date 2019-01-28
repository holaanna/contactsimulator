#' @export
baseline_alt_v0<- function(farm_pos_cat,vis_int_per_cat,farm_inf,farm_inf_on,min_tim_cont,indx_rem,indx_con,mm_1,nn_1){


  #   if(farm_inf[mm_1,nn_1]==0 & (min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])>=360*2){ # Either A or B
  #
  #     if((min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])>10000  ){   # Never recorded
  #       eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "A"))
  #       indx<- which(vis_int_per_cat$cat=="A")
  #
  #     }
  #     else{
  #       eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "B"))
  #       indx<- which(vis_int_per_cat$cat=="B")
  #
  #     }
  #
  #     eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
  #     eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
  #
  #
  #   }
  #   else{# Either C, D or E
  #
  #     k<- (farm_inf_on[mm_1,nn_1]>0) + 1
  #     indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
  #     # show(k)
  #     switch(k,
  #            {  # No infection during the visit
  #
  #
  #                eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
  #
  #
  #            },
  #
  #            { # >=1 infection
  #              kk<- (farm_inf[mm_1,nn_1]>1) + (farm_inf[mm_1,nn_1]>10) + 1
  #              switch(kk,
  #                     { # 1 infection
  #                       eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "C"))
  #                     },
  #                     { # 1-10 infections
  #                       eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "D"))
  #                     },
  #                     { # >10 infections
  #                       eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "E"))
  #                     }
  #              )
  #
  #              indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
  #              eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
  #              eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
  #
  #              eval.parent(substitute(farm_pos_cat$tim_lst_pos[indx_con]<- min_tim_cont))
  #              eval.parent(substitute(farm_pos_cat$nb_round[indx_con]<- 360*2/vis_int_per_cat$vis_int[indx] -1))
  #
  #
  #            }
  #
  #     )
  #
  #     eval.parent(substitute(farm_inf[mm_1,nn_1]<- 0))
  #
  #   }
  #
  #
  # }










  if((min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])<0 ){   # Never recorded
    if(farm_inf[mm_1,nn_1]==0){
      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "A"))
      indx<- which(vis_int_per_cat$cat=="A")
      eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
      eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
    }
    else{
      kk<- (farm_inf[mm_1,nn_1]>1) + (farm_inf[mm_1,nn_1]>10) + 1
      switch(kk,
             { # 1 infection
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "C"))
             },
             { # 1-10 infections
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "D"))
             },
             { # >10 infections
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "E"))
             }
      )
      indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
      eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
      eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))

      eval.parent(substitute(farm_pos_cat$tim_lst_pos[indx_con]<- min_tim_cont))
    }

  }

  if((min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])>=360*2){ # Lat posive more than two years
    if(farm_inf[mm_1,nn_1]==0){
      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "B"))
      indx<- which(vis_int_per_cat$cat=="B")
      eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
      eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
    }
    else{
      kk<- (farm_inf[mm_1,nn_1]>1) + (farm_inf[mm_1,nn_1]>10) + 1
      switch(kk,
             { # 1 infection
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "C"))
             },
             { # 1-10 infections
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "D"))
             },
             { # >10 infections
               eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "E"))
             }
      )
      indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
      eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
      eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))
      eval.parent(substitute(farm_pos_cat$tim_lst_pos[indx_con]<- min_tim_cont))
    }

  }

  if((min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])<=360 & (min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])>0){
    k<- (farm_inf[mm_1,nn_1]>0) + 1
    # show(k)
    switch(k,
           {  # No infection during the visit
             # show(farm_pos_cat$vis[indx_con])
             indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
             eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
             eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))

           },

           { # >=1 infection
             kk<- (farm_inf[mm_1,nn_1]>1) + (farm_inf[mm_1,nn_1]>10) + 1
             switch(kk,
                    { # 1 infection
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "C"))
                    },
                    { # 1-10 infections
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "D"))
                    },
                    { # >10 infections
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "E"))
                    }
             )

             indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
             eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
             eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))

             eval.parent(substitute(farm_pos_cat$tim_lst_pos[indx_con]<- min_tim_cont))

           }

    )

    eval.parent(substitute(farm_inf[mm_1,nn_1]<- 0))

  }


  if((min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])< 2*360 & (min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con])>360){
    k<- (farm_inf[mm_1,nn_1]>0) + 1
    # show(k)
    switch(k,
           {  # No infection during the visit
             # show(farm_pos_cat$vis[indx_con])
             indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
             eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
             eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))


           },

           { # >=1 infection
             kk<- (farm_inf[mm_1,nn_1]>1) + (farm_inf[mm_1,nn_1]>10) + 1
             switch(kk,
                    { # 1 infection
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "C"))
                    },
                    { # 1-10 infections
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "D"))
                    },
                    { # >10 infections
                      eval.parent(substitute(farm_pos_cat$cat[indx_con]<- "E"))
                    }
             )

             indx<- which(vis_int_per_cat$cat==farm_pos_cat$cat[indx_con])
             eval.parent(substitute(farm_pos_cat$vis_int[indx_con]<- vis_int_per_cat$vis_int[indx]))
             eval.parent(substitute(farm_pos_cat$vis[indx_con]<- min_tim_cont+ vis_int_per_cat$vis_int[indx]))

             eval.parent(substitute(farm_pos_cat$tim_lst_pos[indx_con]<- min_tim_cont))

           }

    )

    eval.parent(substitute(farm_inf[mm_1,nn_1]<- 0))

  }

  # show(farm_pos_cat$vst)

}


#

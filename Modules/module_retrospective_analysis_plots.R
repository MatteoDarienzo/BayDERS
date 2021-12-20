


# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
#######################################################################################################################
plot.time.shifts.step1 <- function(dir, 
                                               gaug.1 , rec.1, 
                                               data.annotate.off, 
                                               data.annotate.gaug.1,
                                               data.annotate.rec.1,
                                               data.annotate.step1,
                                               color, 
                                               df.limni, limni.labels, grid_limni.ylim,
                                               pdf.ts.1, pdf.ts.2, pdf.ts.3,
                                               save_file_name) {
  ########################################################################################################
  title1 <- "Segmentation of gaugings"
  title2 <- "Recession analysis"
  title3 <- "Final selection of shift times"  #Combined results of STEP 1"
  title4 <- "Official dates of RC update"
  
  message("- Reading and plotting combined results of Gauging segmentation + Recession analysis.")
  message("- Generating the plot, Please, Wait ... ")
  message(paste0("- Plot available also at ", dir, "/",save_file_name ))
  


  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;  names(gaug.1) = c("h","Q", "uQ", "Period", "t", "t.real")
  if (!is.null(data.annotate.step1$t.adj)) {
    for (i in 1:length(gaug.1$t)) {
      if(gaug.1$t[i] <= data.annotate.step1$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = color[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.step1$t.adj)+1)) {
      for (i in 1:length(gaug.1$t)) {
        if ((gaug.1$t[i] <= tail(gaug.1$t,1)) & (gaug.1$t[i] >  data.annotate.step1$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = color[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$t)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = color[1]
      P_Gaug[i] = 1
    }
  }
  
  
  df.RC.step1 <- data.frame(gaug.1$h, gaug.1$Q, gaug.1$uQ, P_Gaug, gaug.1$t, c_Gaug)
  #df.RC.step1 <- gaug.1
  
  names(df.RC.step1) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step1, paste0(dir,"/data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.3, paste0(dir,"/pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step1, paste0(dir,"/shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  
  
  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (!is.null(df.limni)) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray80",size = 0.2)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot= t.plot + 
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step1, aes(x = t , y= h), size = 6, pch =21, fill= df.RC.step1$color) +
    scale_y_continuous(name   = limni.labels[2], expand = c(0,0), limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                       breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
    # geom_vline(xintercept = tail(df.limni$t_limni,1), color="red", size= 2.5)+
    # annotate("text", 
    #          x = tail(df.limni$t_limni,1) + 50,
    #          y = (grid_limni.ylim[2] - start.y.legend)/2,
    #          label= paste0("Initialization"), angle = 90,
    #          color = "red", size=10)
  
  
  
  
  
  #---------------------------------------
  #PLOT 2: gaugings segmentation results
  #---------------------------------------
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1)) {
    X1 =pdf.ts.1
    if (length(data.annotate.gaug.1$t.adj)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  # plot 2:
  t.plot2 <- ggplot() 
  if (!is.null(data.annotate.gaug.1)) {
      if (length(data.annotate.gaug.1$MAP)==1) {
        t.plot2 = t.plot2 + 
          geom_density(aes(x= X1.new, ..scaled.. ),
                       fill="blue", 
                       colour=NA, alpha=0.3)
      } else {
        
        t.plot2 = t.plot2 +
          geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                       fill= "blue", 
                       colour=NA, alpha=0.3)
      }
      t.plot2 = t.plot2 +
        geom_segment(data  = data.annotate.gaug.1, 
                     aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                     size = 0.7, color ="red")
    } else {
      if (is.null(df.limni)==FALSE) {
        t.plot2 = t.plot2 + 
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2,  
                   y = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2 , 
                   y = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      } else {
        t.plot2 = t.plot2 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x     = tail(g.s$X.tP.,1)/2, 
                   y     = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      }
    }
    
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
      #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
      # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
      #geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "blue", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title1)
    if (is.null(df.limni)==FALSE) {
      # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  
  
  
  
  
  

  
    
    
    
    
    
    
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (!is.null(data.annotate.rec.1)) {
      X2 =pdf.ts.2
      if (length(data.annotate.rec.1$t.adj)==1) {
        X2.new = X2
      } else {
        X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
        for (orderr in 1:(ncol(pdf.ts.2))) {
          X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
        }
      }
  }
    
  t.plot3 <- ggplot() 
  if (!is.null(data.annotate.rec.1)) {
      if (length(data.annotate.rec.1$MAP)==1) {
        t.plot3 = t.plot3 + 
          geom_density(aes(x= X2.new, ..scaled..),
                       fill="green", 
                       colour=NA, alpha=0.3)
      } else {
        
        t.plot3 = t.plot3 +
          geom_density(aes(x= X2.new$time, ..scaled.., group =X2.new$ord),  
                       fill= "green", 
                       colour=NA, alpha=0.3)
      }
      t.plot3 = t.plot3 +
        geom_segment(data  = data.annotate.rec.1, 
                     aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                     size = 0.7, color ="green")
  } else {
      if (is.null(df.limni)==FALSE) {
        t.plot3 = t.plot3 + 
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2,  
                   y = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2 , 
                   y = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      } else {
        t.plot3 = t.plot3 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x     = tail(g.s$X.tP.,1)/2, 
                   y     = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      }
  }  
  t.plot3 = t.plot3 +
    #geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="green")+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "green", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title2)
  if (is.null(df.limni)==FALSE) {
    # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }

    
    
    
  
  
  
  
  
  
  
  
  #-----------------------------------------
  #PLOT 4: combined shift times
  #-----------------------------------------
  if (!is.null(data.annotate.step1)) {
    X3 =pdf.ts.3
    if (length(data.annotate.step1$t.adj)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  t.plot4 <- ggplot() 
  if (!is.null(data.annotate.step1)) {
      if (length(data.annotate.step1$t.adj)==1) {
        t.plot4 = t.plot4 + 
          geom_density(aes(x= X3.new, ..scaled..),
                       fill="red", 
                       colour=NA, alpha=0.3)
      } else {
        
        t.plot4 = t.plot4 +
          geom_density(aes(x= X3.new$time, ..scaled.., group =X3.new$ord),  
                       fill= "red", 
                       colour=NA, alpha=0.3)
      }
      t.plot4 = t.plot4 +
        geom_segment(data  = data.annotate.step1, 
                     aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                     size = 0.7, color ="red")
    } else {
      if (is.null(df.limni)==FALSE) {
        t.plot4 = t.plot4 + 
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2,  
                   y = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2 , 
                   y = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      } else {
        t.plot4 = t.plot4 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label = "No shift detected (or not available)", color = "gray", size=8) +
          annotate("text",
                   x     = tail(g.s$X.tP.,1)/2, 
                   y     = - 0.4, 
                   label = "No shift adjusted (or not available)", color = "gray", size=8) 
      }
    }  
    
    
    t.plot4 = t.plot4 + 
      #geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="red")+
      #annotate("text", x =data.annotate.step1$t.adj, y = -0.8, label = seq(1, length(data.annotate.step1$t.adj),1), color = "black", size = 3) +
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "red", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title3)
    if (is.null(df.limni)==FALSE) {
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    } else {
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
 
  
    
    
    
    
  
  #----------------------
  #PLOT 6: official dates
  #----------------------
  t.plot6 <- ggplot() 
  t.plot6 = t.plot6 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text         = element_text(size=20)
          ,axis.title       = element_text(size=30,face="bold")
          ,plot.title       = element_text(size = 15, color= "black", vjust=-6)
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x      = element_blank()
          ,axis.text.y      = element_blank()
          ,axis.ticks.y     = element_blank()
          ,axis.ticks.x     = element_blank()) +
    ggtitle(title4)
  if (is.null(df.limni)==FALSE) {
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  if (exists("data.annotate.off")==TRUE)  {
    if (!is.null(data.annotate.off)) {
      t.plot6 = t.plot6 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
    } else {
      if (is.null(df.limni) == FALSE) {
        t.plot6 = t.plot6 + 
          annotate("text",
                   x    = tail(df.limni$t_limni,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
      } else {
        t.plot6 = t.plot6 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
        
      }
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot7 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot6, 
                       ncol = 1, nrow = 4, rel_heights = c(1,1,1,1)) 
  t.plot8 = plot_grid( t.plot, 
                       t.plot7, 
                       ncol = 1, nrow = 2,  rel_heights = c(1,1))
  # ggsave(t.plot8, filename =paste0(dir,"/STEP1_time_series.png"), 
  #                 device  = "png", width = 20, height =16, dpi = 400, units = "in")
  
  
  
  
  t.plot2bis = t.plot2  + theme(plot.title = element_text(face="bold", size = 25, color= "blue"))
  t.plot3bis = t.plot3  + theme(plot.title = element_text(face="bold", size = 25, color= "green"))
  t.plot4bis = t.plot4  + theme(plot.title = element_text(face="bold", size = 25, color= "red"))                              
  t.plot6bis = t.plot6  + theme(plot.title = element_text(face="bold", size = 25, color= "black", vjust=0.5))                              
  t.plot7bis = plot_grid( t.plot2bis,  
                          t.plot3bis, 
                          t.plot4bis,
                          t.plot6bis, 
                          ncol = 1, nrow = 4, 
                          rel_heights = c(1,1,1,0.8)) 
  t.plot8bis = plot_grid( t.plot, 
                          t.plot7bis, 
                          ncol = 1,
                          nrow = 2,  
                          rel_heights = c(1.2,1))
  pdf(paste0(dir, save_file_name), 20,13 ,useDingbats=F)
  print(t.plot8bis)
  dev.off()
  
  # ggsave(t.plot8bis, filename =paste0(dir,"/STEP1_time_series.png"), 
  #                  device  = "png", width = 20, height =15, dpi = 400, units = "in")
  # 
  # 
  
  return(t.plot8)
  #------------------------------------------------------------------------
}























# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
#######################################################################################################################
plot.time.shifts.file <- function(dir, 
                                  gaug.1,  
                                  data.annotate.off, 
                                  data.annotate.file,
                                  data.annotate.chosen,
                                  color,  df.limni, limni.labels, grid_limni.ylim,
                                  pdf.ts.file, pdf.ts.chosen,
                                  save_file_name) {
  ########################################################################################################
  title1 <- "Proposed segmentation"
  title2 <- "Final selection of shift times"
  title3 <- "Official dates of RC update"
  
  message("Reading and plotting results provided in the text file.")
  message("Generating the plot, Please, Wait ... ")
  message(paste0("Plot available also at ", dir, "/",save_file_name ))
  
  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.chosen$t.adj)) {
    for (i in 1:length(gaug.1$t)) {
      if(gaug.1$t[i] <= data.annotate.chosen$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = color[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.chosen$t.adj)+1)) {
      for (i in 1:length(gaug.1$t)) {
        if ((gaug.1$t[i] <= tail(gaug.1$t,1)) & 
            (gaug.1$t[i] >  data.annotate.chosen$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = color[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = color[1]
      P_Gaug[i] = 1
    }
  }
  
  df.RC.step1 <- data.frame(gaug.1$h, gaug.1$Q, gaug.1$uQ, P_Gaug, gaug.1$t, c_Gaug)
  #df.RC.step1 <- gaug.1
  names(df.RC.step1) = c("h","Q", "uQ", "Period", "t","color")

  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (!is.null(df.limni)) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray80",size = 0.2)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot= t.plot + 
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step1, aes(x = t , y= h), size = 6, pch =21, fill= df.RC.step1$color) +
    scale_y_continuous(name   = limni.labels[2], expand = c(0,0), limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                       breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
  
  
  
  #---------------------------------------
  #PLOT 2: gaugings segmentation results
  #---------------------------------------
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.file)) {
    X1 =pdf.ts.file
    if (length(data.annotate.file$t.adj)==1) {
      X1.new = X1
    } else {
      X1.new = data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.file))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  # plot 2:
  t.plot2 <- ggplot() 
  if (!is.null(data.annotate.file)) {
    if (length(data.annotate.file$MAP)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1.new, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue",  colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 +
      geom_segment(data  = data.annotate.file, 
                   aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                   size = 0.7, color ="red")
  } else {
    if (is.null(df.limni)==FALSE) {
      t.plot2 = t.plot2 + 
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2,  
                 y = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2 , 
                 y = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    } else {
      t.plot2 = t.plot2 + 
        annotate("text",
                 x    = tail(g.s$X.tP.,1)/2,  
                 y    = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x     = tail(g.s$X.tP.,1)/2, 
                 y     = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    }
  }
  
  t.plot2 = t.plot2 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "blue", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title1)
  if (is.null(df.limni)==FALSE) {
    t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  
  

  #------------------------------------------
  #PLOT 3: Final chosen segmentation results
  #------------------------------------------
  if (!is.null(data.annotate.chosen)) {
    X2 =pdf.ts.chosen
    if (length(data.annotate.chosen$t.adj)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.chosen))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  
  t.plot3 <- ggplot() 
  if (!is.null(data.annotate.chosen)) {
    if (length(data.annotate.chosen$MAP)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2.new, ..scaled..),
                     fill="red",  colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.., group =X2.new$ord),  
                     fill= "red", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data  = data.annotate.chosen, 
                   aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                   size = 0.7, color ="red")
  } else {
    if (is.null(df.limni)==FALSE) {
      t.plot3 = t.plot3 + 
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2,  
                 y = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2 , 
                 y = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    } else {
      t.plot3 = t.plot3 + 
        annotate("text",
                 x    = tail(g.s$X.tP.,1)/2,  
                 y    = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x     = tail(g.s$X.tP.,1)/2, 
                 y     = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    }
  }  
  t.plot3 = t.plot3 +
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "green", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title2)
  if (is.null(df.limni)==FALSE) {
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  
  


  
  
  
  #----------------------
  #PLOT 6: official dates
  #----------------------
  t.plot4 <- ggplot() 
  t.plot4 = t.plot4 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text         = element_text(size=20)
          ,axis.title       = element_text(size=30,face="bold")
          ,plot.title       = element_text(size = 15, color= "black", vjust=-6)
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x      = element_blank()
          ,axis.text.y      = element_blank()
          ,axis.ticks.y     = element_blank()
          ,axis.ticks.x     = element_blank()) +
    ggtitle(title3)
  if (is.null(df.limni)==FALSE) {
    t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  if (exists("data.annotate.off")==TRUE)  {
    if (!is.null(data.annotate.off)) {
      t.plot4 = t.plot4 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
    } else {
      if (is.null(df.limni) == FALSE) {
        t.plot4 = t.plot4 + 
          annotate("text",
                   x    = tail(df.limni$t_limni,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
      } else {
        t.plot4 = t.plot4 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
        
      }
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot7 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4, 
                       ncol = 1, nrow = 3, rel_heights = c(1,1,1)) 
  t.plot8 = plot_grid( t.plot, 
                       t.plot7, 
                       ncol = 1, nrow = 2,  rel_heights = c(1,1))
  
  t.plot2bis = t.plot2  + theme(plot.title = element_text(face="bold", size = 25, color= "blue"))
  t.plot3bis = t.plot3  + theme(plot.title = element_text(face="bold", size = 25, color= "red"))                              
  t.plot4bis = t.plot4  + theme(plot.title = element_text(face="bold", size = 25, color= "black", vjust=0.5))                              
  t.plot7bis = plot_grid( t.plot2bis,  
                          t.plot3bis, 
                          t.plot4bis,
                          ncol = 1, nrow = 3, 
                          rel_heights = c(1,1,0.8)) 
  t.plot8bis = plot_grid( t.plot, 
                          t.plot7bis, 
                          ncol = 1,
                          nrow = 2,  
                          rel_heights = c(1.2,1))
  pdf(paste0(dir, save_file_name), 20,10 ,useDingbats=F)
  print(t.plot8bis)
  dev.off()
  
  
  return(t.plot8)
}



















# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
######################################################################################################
plot.reference.shift.times <- function(dir, 
                                   gaug.1 , rec.1, 
                                   data.annotate.off, 
                                   data.annotate.gaug.1,
                                   data.annotate.rec.1,
                                   data.annotate.step1,
                                   color, 
                                   df.limni, limni.labels, grid_limni.ylim,
                                   pdf.ts.1, pdf.ts.2, pdf.ts.3, save_file_name) {
  ########################################################################################################
  title1 <- "Segmentation of gaugings"
  title2 <- "Recession analysis"
  title3 <- "Final selection of shift times"  #Combined results of STEP 1"
  title4 <- "Official dates of RC update"
  title5 <- "Reference times of morphological shifts"
  
  message("Reading and plotting all results of effecive shift detection.")
  message("Generating the plot, Please, Wait ... ")
  message(paste0("Plot available also at ", dir, "/",save_file_name ))
  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.step1$t.adj)) {
    for (i in 1:length(gaug.1$t)) {
      if(gaug.1$t[i] <= data.annotate.step1$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = colo[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.step1$t.adj)+1)) {
      for (i in 1:length(gaug.1$t)) {
        if ((gaug.1$t[i] <= tail(gaug.1$t,1)) & 
            (gaug.1$t[i] >  data.annotate.step1$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = colo[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = colo[1]
      P_Gaug[i] = 1
    }
  }
  
  
  df.RC.step1 <- data.frame(gaug.1$h, gaug.1$Q, gaug.1$uQ, P_Gaug, gaug.1$t, c_Gaug)
  #df.RC.step1 <- gaug.1
  
  names(df.RC.step1) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step1, paste0(dir,"/STEP1_data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.3, paste0(dir,"/STEP1_pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step1, paste0(dir,"/STEP1_shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  
  
  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (!is.null(df.limni)) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray80",size = 0.2)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot= t.plot + 
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step1, aes(x = t , y= h), size = 6, pch =21, fill= df.RC.step1$color) +
    scale_y_continuous(name   = limni.labels[2], expand = c(0,0), limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                       breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
  # geom_vline(xintercept = tail(df.limni$t_limni,1), color="red", size= 2.5)+
  # annotate("text", 
  #          x = tail(df.limni$t_limni,1) + 50,
  #          y = (grid_limni.ylim[2] - start.y.legend)/2,
  #          label= paste0("Initialization"), angle = 90,
  #          color = "red", size=10)
  
  
  
  
  
  #---------------------------------------
  #PLOT 2: gaugings segmentation results
  #---------------------------------------
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1)) {
    X1 =pdf.ts.1
    if (length(data.annotate.gaug.1$t.adj)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  # plot 2:
  t.plot2 <- ggplot() 
  if (!is.null(data.annotate.gaug.1)) {
    if (length(data.annotate.gaug.1$MAP)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1.new, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 +
      geom_segment(data  = data.annotate.gaug, 
                   aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                   size = 0.7, color ="red")
  } else {
    if (is.null(df.limni)==FALSE) {
      t.plot2 = t.plot2 + 
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2,  
                 y = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2 , 
                 y = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    } else {
      t.plot2 = t.plot2 + 
        annotate("text",
                 x    = tail(g.s$X.tP.,1)/2,  
                 y    = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x     = tail(g.s$X.tP.,1)/2, 
                 y     = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    }
  }
  
  t.plot2 = t.plot2 + 
    # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
    #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
    # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
    #geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "blue", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title1)
  if (is.null(df.limni)==FALSE) {
    # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (!is.null(data.annotate.rec.1)) {
    X2 =pdf.ts.2
    if (length(data.annotate.rec.1$t.adj)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.2))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  
  t.plot3 <- ggplot() 
  if (!is.null(data.annotate.rec.1)) {
    if (length(data.annotate.rec.1$MAP)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2.new, ..scaled..),
                     fill="green", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.., group =X2.new$ord),  
                     fill= "green", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data  = data.annotate.rec.1, 
                   aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                   size = 0.7, color ="green")
  } else {
    if (is.null(df.limni)==FALSE) {
      t.plot3 = t.plot3 + 
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2,  
                 y = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2 , 
                 y = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    } else {
      t.plot3 = t.plot3 + 
        annotate("text",
                 x    = tail(g.s$X.tP.,1)/2,  
                 y    = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x     = tail(g.s$X.tP.,1)/2, 
                 y     = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    }
  }  
  t.plot3 = t.plot3 +
    #geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="green")+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "green", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title2)
  if (is.null(df.limni)==FALSE) {
    # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot3= t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  #-----------------------------------------
  #PLOT 4: combined shift times
  #-----------------------------------------
  if (!is.null(data.annotate.step1)) {
    X3 =pdf.ts.3
    if (length(data.annotate.step1$t.adj)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  t.plot4 <- ggplot() 
  if (!is.null(data.annotate.step1)) {
    if (length(data.annotate.step1$t.adj)==1) {
      t.plot4 = t.plot4 + 
        geom_density(aes(x= X3.new, ..scaled..),
                     fill="red", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot4 = t.plot4 +
        geom_density(aes(x= X2.new$time, ..scaled.., group =X2.new$ord),  
                     fill= "red", 
                     colour=NA, alpha=0.3)
    }
    t.plot4 = t.plot4 +
      geom_segment(data  = data.annotate.step1, 
                   aes(x = t.adj, y = -1, yend = 0, xend=t.adj), 
                   size = 0.7, color ="red")
  } else {
    if (is.null(df.limni)==FALSE) {
      t.plot4 = t.plot4 + 
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2,  
                 y = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x = tail(df.limni$t_limni,1)/2 , 
                 y = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    } else {
      t.plot4 = t.plot4 + 
        annotate("text",
                 x    = tail(g.s$X.tP.,1)/2,  
                 y    = 0.4, 
                 label = "No shift detected (or not available)", color = "gray", size=8) +
        annotate("text",
                 x     = tail(g.s$X.tP.,1)/2, 
                 y     = - 0.4, 
                 label = "No shift adjusted (or not available)", color = "gray", size=8) 
    }
  }  
  
  
  t.plot4 = t.plot4 + 
    #geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="red")+
    #annotate("text", x =data.annotate.step1$t.adj, y = -0.8, label = seq(1, length(data.annotate.step1$t.adj),1), color = "black", size = 3) +
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
    theme_set(theme_gray(base_size = 20))+
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(face = "bold", size = 15, color= "red", margin=margin(0,0,0,0))
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title3)
  if (is.null(df.limni)==FALSE) {
    t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  
  
  
  
  
  #----------------------
  #PLOT 6: official dates
  #----------------------
  t.plot6 <- ggplot() 
  t.plot6 = t.plot6 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text         = element_text(size=20)
          ,axis.title       = element_text(size=30,face="bold")
          ,plot.title       = element_text(size = 15, color= "black", vjust=-6)
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x      = element_blank()
          ,axis.text.y      = element_blank()
          ,axis.ticks.y     = element_blank()
          ,axis.ticks.x     = element_blank()) +
    ggtitle(title4)
  if (is.null(df.limni)==FALSE) {
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  
  
  if (exists("data.annotate.off")==TRUE)  {
    if (!is.null(data.annotate.off)) {
      t.plot6 = t.plot6 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
    } else {
      if (is.null(df.limni) == FALSE) {
        t.plot6 = t.plot6 + 
          annotate("text",
                   x    = tail(df.limni$t_limni,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
      } else {
        t.plot6 = t.plot6 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.4, 
                   label= "Official shift dates not available", color = "gray", size=8)
        
      }
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot7 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot6, 
                       ncol = 1, nrow = 4, rel_heights = c(1,1,1,1)) 
  t.plot8 = plot_grid( t.plot, 
                       t.plot7, 
                       ncol = 1, nrow = 2,  rel_heights = c(1,1))
  # ggsave(t.plot8, filename =paste0(dir,"/STEP1_time_series.png"), 
  #                 device  = "png", width = 20, height =16, dpi = 400, units = "in")
  
  
  
  
  t.plot2bis = t.plot2  + theme(plot.title = element_text(face="bold", size = 25, color= "blue"))
  t.plot3bis = t.plot3  + theme(plot.title = element_text(face="bold", size = 25, color= "green"))
  t.plot4bis = t.plot4  + theme(plot.title = element_text(face="bold", size = 25, color= "red"))                              
  t.plot6bis = t.plot6  + theme(plot.title = element_text(face="bold", size = 25, color= "black", vjust=0.5))                              
  t.plot7bis = plot_grid( t.plot2bis,  
                          t.plot3bis, 
                          t.plot4bis,
                          t.plot6bis, 
                          ncol = 1, nrow = 4, 
                          rel_heights = c(1,1,1,0.8)) 
  t.plot8bis = plot_grid( t.plot, 
                          t.plot7bis, 
                          ncol = 1,
                          nrow = 2,  
                          rel_heights = c(1.2,1))
  pdf(paste0(dir, save_file_name), 20,13 ,useDingbats=F)
  print(t.plot8bis)
  dev.off()
  
  # ggsave(t.plot8bis, filename =paste0(dir,"/STEP1_time_series.png"), 
  #                  device  = "png", width = 20, height =15, dpi = 400, units = "in")
  # 
  # 
  
  return(t.plot8)
  #------------------------------------------------------------------------
}































# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
##############################################################################################################
plot.time.shifts.step1.and.bt <- function(dir, dir.SPD.results,
                                   gaug.1 , rec.1, 
                                   data.annotate.off, 
                                   data.annotate.gaug.1,
                                   data.annotate.rec.1,
                                   data.annotate.step1,
                                   colo, 
                                   df.limni, limni.labels , start.y.legend, grid_limni.ylim, 
                                   pdf.ts.1, pdf.ts.2, pdf.ts.3, 
                                   limits.time
) {
  ###########################################################################################################
  title1 <- "Segmentation of gaugings"
  title2 <- "Recession analysis"
  title3 <- "Combined results of STEP 1"
  title4 <- "Official dates of RC update"
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1$MAP)) {
    X1 =pdf.ts.1
    if (ncol(pdf.ts.1)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.rec.1$MAP)) {
    X2 =pdf.ts.2
    if (ncol(pdf.ts.2)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.2))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step1$MAP)) {
    X3 =pdf.ts.3
    if (ncol(pdf.ts.3)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  
  
  #preparing gauging data per periods:
  nperiod.step1 = length(data.annotate.step1$t.adj)+1
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.step1$t.adj)) {
    for (i in 1:length(gaug.1$X.tP.)) {
      if(gaug.1$X.tP.[i] <= data.annotate.step1$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = colo[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:nperiod.step1) {
      for (i in 1:length(gaug.1$X.tP.)) {
        if ((gaug.1$X.tP.[i] <= tail(gaug.1$X.tP.,1)) & 
            (gaug.1$X.tP.[i] >  data.annotate.step1$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = colo[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = colo[1]
      P_Gaug[i] = 1
    }
  }
  df.RC.step1 <- data.frame(gaug.1$X.h., gaug.1$X.Q., gaug.1$X.uQ., P_Gaug, gaug.1$X.tP., c_Gaug)
  names(df.RC.step1) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step1, paste0(dir,"/STEP1_data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.3, paste0(dir,"/STEP1_pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step1, paste0(dir,"/STEP1_shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  #reading data of spd results (in particular the riverbed elevation b2 for each period:
  
  b1t= read.table(paste0(dir.SPD.STEP1.results, "/Results_MCMC_Cooked.txt"), header=TRUE)[,1:nperiod.step1]
  b2t = read.table(paste0(dir.SPD.STEP1.results, "/Results_MCMC_Cooked.txt"), header=TRUE)[, (nperiod.step1+3):(nperiod.step1*2+3)]
  
  
  time.adjust.before = c(0, data.annotate.step1$t.adj)
  time.adjust.plus = c(data.annotate.step1$t.adj, limits.time[2])
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "black",size = 0.3)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.time)
  } else {
    t.plot= t.plot + scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step1, aes(x = t , y= h), size = 5, pch =21, fill= df.RC.step1$color)+ #g.1$color)+
    #geom_vline(aes(xintercept = data.annotate.gaug.adjust.1$t.adj), color = "red", lwd =0.3, linetype = "dashed")+
    geom_segment(mapping=aes(x = time.adjust.before , 
                             y = read.res.rec$mu.res/100, 
                             xend = time.adjust.plus, 
                             yend = read.res.rec$mu.res/100), col="red") +
    # annotate("rect",  
    #          xmin = time.adjust.before, 
    #          xmax = time.adjust.plus,
    #          ymin=(read.res.rec$Q2.mu- read.res.rec$gamma[1])/100,  
    #          ymax=(read.res.rec$Q97.mu + read.res.rec$gamma[1])/100, 
    #          fill="red", alpha=0.5)+
    annotate("rect",  
             xmin= time.adjust.before, 
             xmax= time.adjust.plus,
             ymin= ((read.res.rec$Q2.mu)/100),  
             ymax=((read.res.rec$Q97.mu)/100), 
             fill="pink", alpha=0.9)
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(start.y.legend, grid_limni.ylim[2]), 
                       breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
  
  #------------------------------------------
  #PLOT 2: gaugings segmentation results
  #------------------------------------------
  if (exists("data.annotate.gaug.1")) {
    t.plot2 <- ggplot() 
    if (ncol(pdf.ts.1)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
      #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
      # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
      geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "blue", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title1)
    if (is.null(df.limni)==FALSE) {
      # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (exists("data.annotate.rec.1")) {
    t.plot3 <- ggplot() 
    if (ncol(pdf.ts.2)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2, ..scaled.. ),
                     fill="green", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.. , group =X2.new$ord),  
                     fill= "green", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="green")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "green", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title2)
    if (is.null(df.limni)==FALSE) {
      # t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 4: method B
  #-----------------------------------------
  if ((exists("data.annotate.step1"))) {
    t.plot4 <- ggplot() 
    if (ncol(pdf.ts.3)==1) {
      t.plot4 = t.plot4 + 
        geom_density(aes(x= X3, ..scaled..),
                     fill="red", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot4 = t.plot4 +
        geom_density(aes(x= X3.new$time, ..scaled.., group =X3.new$ord),  
                     fill= "red",
                     colour=NA, alpha=0.3)
    }
    t.plot4 = t.plot4 + 
      geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="red")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "red", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title3)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #----------------------
  #PLOT 6: official dates
  #----------------------
  t.plot6 <- ggplot() 
  t.plot6 = t.plot6 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(size = 25, color= "red", vjust=-6)
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title4)
  if (is.null(df.limni)==FALSE) {
    # t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
  } else {
    t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null(data.annotate.off)==FALSE) {
      t.plot6 <- t.plot6 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
      # geom_point(data = data.annotate.off, aes(x = xpotent, y = 0), color= "green", size = 4, shape =4, stroke=2 )
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot7 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot6, 
                       ncol = 1, nrow = 4, rel_heights = c(1,1,1,1)) 
  t.plot8 = plot_grid( t.plot, 
                       t.plot7, 
                       ncol = 1, nrow = 2,  rel_heights = c(2,1))
  ggsave(t.plot8, filename =paste0(dir,"/STEP1_time_series.png"), 
         device  = "png", width = 20, height =16, dpi = 400, units = "in")
  return(t.plot8)
  #------------------------------------------------------------------------
}

































##################################################################################################
#                                          STEP 2
##################################################################################################



# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
#######################################################################################################
plot.time.shifts.step2 <- function(dir, 
                                   gaug.1 , rec.1, 
                                   data.annotate.off, 
                                   data.annotate.gaug.1,
                                   data.annotate.rec.1,
                                   data.annotate.step1,
                                   data.annotate.ST,
                                   data.annotate.step2,
                                   colo, 
                                   df.limni, limni.labels , start.y.legend, grid_limni.ylim,
                                   pdf.ts.1, 
                                   pdf.ts.2,
                                   pdf.ts.3,
                                   pdf.ts.4,
                                   pdf.ts.5,
                                   limits.time
) {
  ###########################################################################################################
  title1 <- "Segmentation of gaugings"
  title2 <- "Recession analysis"
  title3 <- "Combined results of STEP 1"
  title4 <- "Sediment transport analysis (merged morphogenic events)"
  title5 <- "Combined results of STEP 2"
  title6 <- "Official dates of RC update"
  
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1$MAP)) {
    X1 =pdf.ts.1
    if (ncol(pdf.ts.1)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.rec.1$MAP)) {
    X2 =pdf.ts.2
    if (ncol(pdf.ts.2)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.2))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step1$MAP)) {
    X3 =pdf.ts.3
    if (ncol(pdf.ts.3)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.ST$MAP)) {
    X4 =pdf.ts.4
    if (ncol(pdf.ts.4)==1) {
      X4.new = X4
    } else {
      X4.new =data.frame(time= X4[,1], ord = rep(1, length(X4[,1])))
      for (orderr in 1:(ncol(pdf.ts.4))) {
        X4.new =rbind(X4.new, data.frame(time= X4[,orderr], ord = rep(orderr, length(X4[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step2$MAP)) {
    X5 =pdf.ts.5
    if (ncol(pdf.ts.5)==1) {
      X5.new = X5
    } else {
      X5.new =data.frame(time= X5[,1], ord = rep(1, length(X5[,1])))
      for (orderr in 1:(ncol(pdf.ts.5))) {
        X5.new =rbind(X5.new, data.frame(time= X5[,orderr], ord = rep(orderr, length(X5[,1]))))
      }
    }
  }
  
  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.step2$t.adj)) {
    for (i in 1:length(gaug.1$X.tP.)) {
      if(gaug.1$X.tP.[i] <= data.annotate.step2$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = colo[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.step2$t.adj)+1)) {
      for (i in 1:length(gaug.1$X.tP.)) {
        if ((gaug.1$X.tP.[i] <= tail(gaug.1$X.tP.,1)) & 
            (gaug.1$X.tP.[i] >  data.annotate.step2$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = colo[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = colo[1]
      P_Gaug[i] = 1
    }
  }
  df.RC.step2 <- data.frame(gaug.1$X.h., gaug.1$X.Q., gaug.1$X.uQ., P_Gaug, gaug.1$X.tP., c_Gaug)
  names(df.RC.step2) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step2, paste0(dir,"/STEP2_data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.5, paste0(dir,"/STEP2_pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step2, paste0(dir,"/STEP2_shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray40",size = 0.3)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.time)
  } else {
    t.plot= t.plot + scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step2, aes(x = t , y= h), size = 5, pch =21, fill= df.RC.step2$color)+ #g.1$color)+
    #geom_vline(aes(xintercept = data.annotate.gaug.adjust.1$t.adj), color = "red", lwd =0.3, linetype = "dashed")+
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(start.y.legend, grid_limni.ylim[2]), 
                       breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
  
  #------------------------------------------
  #PLOT 2: gaugings segmentation results
  #------------------------------------------
  if (exists("data.annotate.gaug.1")) {
    t.plot2 <- ggplot() 
    if (ncol(pdf.ts.1)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
      #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
      # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
      geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "blue", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title1)
    if (is.null(df.limni)==FALSE) {
      # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (exists("data.annotate.rec.1")) {
    t.plot3 <- ggplot() 
    if (ncol(pdf.ts.2)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2, ..scaled.. ),
                     fill="green", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.. , group =X2.new$ord),  
                     fill= "green", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="green")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "green", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title2)
    if (is.null(df.limni)==FALSE) {
      # t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 4: results STEP 1
  #-----------------------------------------
  if ((exists("data.annotate.step1"))) {
    t.plot4 <- ggplot() 
    if (ncol(pdf.ts.3)==1) {
      t.plot4 = t.plot4 + 
        geom_density(aes(x= X3, ..scaled..),
                     fill="red", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot4 = t.plot4 +
        geom_density(aes(x= X3.new$time, ..scaled.., group =X3.new$ord),  
                     fill= "red",
                     colour=NA, alpha=0.3)
    }
    t.plot4 = t.plot4 + 
      geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="red")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "red", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title3)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 5: results Sediment transport
  #-----------------------------------------
  if ((exists("data.annotate.ST"))) {
    t.plot5 <- ggplot() 
    if (ncol(pdf.ts.4)==1) {
      t.plot5 = t.plot5 + 
        geom_density(aes(x= X4, ..scaled..),
                     fill="violet", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot5 = t.plot5 +
        geom_density(aes(x= X4.new$time, ..scaled.., group =X4.new$ord),  
                     fill= "violet",
                     colour=NA, alpha=0.3)
    }
    t.plot5 = t.plot5 + 
      geom_segment(data = data.annotate.ST, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="violet")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "violet", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title4)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot5= t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot5 = t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 6: results STEP 2
  #-----------------------------------------
  if ((exists("data.annotate.step2"))) {
    t.plot6 <- ggplot() 
    if (ncol(pdf.ts.5)==1) {
      t.plot6 = t.plot6 + 
        geom_density(aes(x= X5, ..scaled..),
                     fill="black", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot6 = t.plot6 +
        geom_density(aes(x= X5.new$time, ..scaled.., group =X5.new$ord),  
                     fill= "black",
                     colour=NA, alpha=0.3)
    }
    t.plot6 = t.plot6 + 
      geom_segment(data = data.annotate.step2, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="black") + 
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1)) +
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank())+
      ggtitle(title5)
    if (is.null(df.limni)==FALSE) {
      #t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  
  #----------------------
  #PLOT 7: official dates
  #----------------------
  t.plot7 <- ggplot() 
  t.plot7 = t.plot7 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(size = 25, color= "black", vjust = -6)
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title6)
  if (is.null(df.limni)==FALSE) {
    # t.plot7 = t.plot7 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot7= t.plot7 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
  } else {
    t.plot7 = t.plot7 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null(data.annotate.off)==FALSE) {
      t.plot7 <- t.plot7 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
      # geom_point(data = data.annotate.off, aes(x = xpotent, y = 0), color= "green", size = 4, shape =4, stroke=2 )
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot8 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot5,
                       t.plot6, 
                       t.plot7, 
                       ncol = 1, nrow = 6, rel_heights = c(1,1,1,1,1,1)) 
  t.plot9 = plot_grid( t.plot, 
                       t.plot8, 
                       ncol = 1, nrow = 2,  rel_heights = c(1, 1))
  ggsave(t.plot9, filename =paste0(dir,"/STEP2_time_series.png"), 
         device  = "png", width = 20, height =16, dpi = 400, units = "in")
  return(t.plot9)
  #------------------------------------------------------------------------
}

































##################################################################################################
#                                          STEP 3:
##################################################################################################

plot.time.shifts.step3 <- function(dir, 
                                   gaug.1 , rec.1, 
                                   data.annotate.off, 
                                   data.annotate.gaug.1,
                                   data.annotate.rec.1,
                                   data.annotate.step1,
                                   data.annotate.ST,
                                   data.annotate.step2,
                                   data.annotate.step3,
                                   colo, 
                                   df.limni, limni.labels , start.y.legend, grid_limni.ylim,
                                   pdf.ts.1, 
                                   pdf.ts.2,
                                   pdf.ts.3,
                                   pdf.ts.4,
                                   pdf.ts.5,
                                   pdf.ts.6,
                                   limits.time
) {
  ###########################################################################################################
  title1 <- "Segmentation of gaugings"
  title2 <- "Recession analysis"
  title3 <- "Combined results of STEP 1"
  title4 <- "Sediment transport analysis (merged morphogenic events)"
  title5 <- "Combined results of STEP 2"
  title6 <- "Combined results of STEP 3"
  title7 <- "Official dates of RC update"
  
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1$MAP)) {
    X1 =pdf.ts.1
    if (ncol(pdf.ts.1)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.rec.1$MAP)) {
    X2 =pdf.ts.2
    if (ncol(pdf.ts.2)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.2))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step1$MAP)) {
    X3 =pdf.ts.3
    if (ncol(pdf.ts.3)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.ST$MAP)) {
    X4 =pdf.ts.4
    if (ncol(pdf.ts.4)==1) {
      X4.new = X4
    } else {
      X4.new =data.frame(time= X4[,1], ord = rep(1, length(X4[,1])))
      for (orderr in 1:(ncol(pdf.ts.4))) {
        X4.new =rbind(X4.new, data.frame(time= X4[,orderr], ord = rep(orderr, length(X4[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step2$MAP)) {
    X5 =pdf.ts.5
    if (ncol(pdf.ts.5)==1) {
      X5.new = X5
    } else {
      X5.new =data.frame(time= X5[,1], ord = rep(1, length(X5[,1])))
      for (orderr in 1:(ncol(pdf.ts.5))) {
        X5.new =rbind(X5.new, data.frame(time= X5[,orderr], ord = rep(orderr, length(X5[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step3$MAP)) {
    X6 =pdf.ts.6
    if (ncol(pdf.ts.6)==1) {
      X6.new = X6
    } else {
      X6.new =data.frame(time= X6[,1], ord = rep(1, length(X6[,1])))
      for (orderr in 1:(ncol(pdf.ts.6))) {
        X6.new =rbind(X6.new, data.frame(time= X6[,orderr], ord = rep(orderr, length(X6[,1]))))
      }
    }
  }
  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.step3$t.adj)) {
    for (i in 1:length(gaug.1$X.tP.)) {
      if(gaug.1$X.tP.[i] <= data.annotate.step3$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = colo[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.step3$t.adj)+1)) {
      for (i in 1:length(gaug.1$X.tP.)) {
        if ((gaug.1$X.tP.[i] <= tail(gaug.1$X.tP.,1)) & 
            (gaug.1$X.tP.[i] >  data.annotate.step3$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = colo[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = colo[1]
      P_Gaug[i] = 1
    }
  }
  df.RC.step3 <- data.frame(gaug.1$X.h., gaug.1$X.Q., gaug.1$X.uQ., P_Gaug, gaug.1$X.tP., c_Gaug)
  names(df.RC.step3) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step3, paste0(dir,"/STEP3_data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.6, paste0(dir,"/STEP3_pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step3, paste0(dir,"/STEP3_shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "black",size = 0.3)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.time)
  } else {
    t.plot= t.plot + scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step3, aes(x = t , y= h), size = 5, pch =21, fill= df.RC.step3$color)+ #g.1$color)+
    #geom_vline(aes(xintercept = data.annotate.gaug.adjust.1$t.adj), color = "red", lwd =0.3, linetype = "dashed")+
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(start.y.legend, grid_limni.ylim[2]), 
                       breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())
  
  #------------------------------------------
  #PLOT 2: gaugings segmentation results
  #------------------------------------------
  if (exists("data.annotate.gaug.1")) {
    t.plot2 <- ggplot() 
    if (ncol(pdf.ts.1)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
      #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
      # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
      geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "blue", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title1)
    if (is.null(df.limni)==FALSE) {
      # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (exists("data.annotate.rec.1")) {
    t.plot3 <- ggplot() 
    if (ncol(pdf.ts.2)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2, ..scaled.. ),
                     fill="green", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.. , group =X2.new$ord),  
                     fill= "green", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="green")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "green", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title2)
    if (is.null(df.limni)==FALSE) {
      # t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 4: results STEP 1
  #-----------------------------------------
  if ((exists("data.annotate.step1"))) {
    t.plot4 <- ggplot() 
    if (ncol(pdf.ts.3)==1) {
      t.plot4 = t.plot4 + 
        geom_density(aes(x= X3, ..scaled..),
                     fill="red", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot4 = t.plot4 +
        geom_density(aes(x= X3.new$time, ..scaled.., group =X3.new$ord),  
                     fill= "red",
                     colour=NA, alpha=0.3)
    }
    t.plot4 = t.plot4 + 
      geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="red")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "red", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title3)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 5: results Sediment transport
  #-----------------------------------------
  if ((exists("data.annotate.ST"))) {
    t.plot5 <- ggplot() 
    if (ncol(pdf.ts.4)==1) {
      t.plot5 = t.plot5 + 
        geom_density(aes(x= X4, ..scaled..),
                     fill="violet", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot5 = t.plot5 +
        geom_density(aes(x= X4.new$time, ..scaled.., group =X4.new$ord),  
                     fill= "violet",
                     colour=NA, alpha=0.3)
    }
    t.plot5 = t.plot5 + 
      geom_segment(data = data.annotate.ST, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="violet")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "violet", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title4)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot5= t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot5 = t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 6: results STEP 2
  #-----------------------------------------
  if ((exists("data.annotate.step2"))) {
    t.plot6 <- ggplot() 
    if (ncol(pdf.ts.5)==1) {
      t.plot6 = t.plot6 + 
        geom_density(aes(x= X5, ..scaled..),
                     fill="black", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot6 = t.plot6 +
        geom_density(aes(x= X5.new$time, ..scaled.., group =X5.new$ord),  
                     fill= "black",
                     colour=NA, alpha=0.3)
    }
    t.plot6 = t.plot6 + 
      geom_segment(data = data.annotate.step2, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="black") + 
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1)) +
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title5)
    if (is.null(df.limni)==FALSE) {
      #t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #-----------------------------------------
  #PLOT 7: results STEP 3
  #-----------------------------------------
  if ((exists("data.annotate.step3"))) {
    t.plot7 <- ggplot() 
    if (ncol(pdf.ts.6)==1) {
      t.plot7 = t.plot7 + 
        geom_density(aes(x= X6, ..scaled..),
                     fill="orange", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot7 = t.plot7 +
        geom_density(aes(x= X6.new$time, ..scaled.., group =X6.new$ord),  
                     fill= "orange",
                     colour=NA, alpha=0.3)
    }
    t.plot7 = t.plot7 + 
      geom_segment(data = data.annotate.step3, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="orange") + 
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1)) +
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 30, color= "orange", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title6)
    if (is.null(df.limni)==FALSE) {
      #t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot7= t.plot7 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot7 = t.plot7 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  #----------------------
  #PLOT 8: official dates
  #----------------------
  t.plot8<- ggplot() 
  t.plot8 = t.plot8 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(size = 25, color= "black", vjust=-6)
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title7)
  if (is.null(df.limni)==FALSE) {
    # t.plot8 = t.plot8 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot8= t.plot8 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
  } else {
    t.plot8 = t.plot8 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null(data.annotate.off)==FALSE) {
      t.plot8 <- t.plot8 +
        geom_point(data = data.annotate.off, aes(x = xeffect, y = 0), color= "black", size = 4, shape =4, stroke=2 )
      # geom_point(data = data.annotate.off, aes(x = xpotent, y = 0), color= "green", size = 4, shape =4, stroke=2 )
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot9 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot5,
                       t.plot6, 
                       t.plot7,
                       t.plot8,
                       ncol = 1, nrow = 7, rel_heights = c(1,1,1,1,1,1,1)) 
  t.plot10 = plot_grid(t.plot, 
                       t.plot9, 
                       ncol = 1, nrow = 2,  rel_heights = c(1,1))
  ggsave(t.plot10, filename =paste0(dir,"/STEP3_time_series.png"), 
         device  = "png", width = 20, height =20, dpi = 400, units = "in")
  return(t.plot8)
  #------------------------------------------------------------------------
}



































# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
#######################################################################################################################
plot.time.shifts.step2.retro <- function(dir, 
                                   gaug.1 , rec.1, 
                                   data.annotate.off, 
                                   data.annotate.gaug.1,
                                   data.annotate.rec.1,
                                   data.annotate.step1,
                                   data.annotate.sed.trasp,
                                   colo, 
                                   df.limni, limni.labels , start.y.legend, grid_limni.ylim,
                                   pdf.ts.1, pdf.ts.2, pdf.ts.3, 
                                   limits.time,
                                   hcritical
) {
  ########################################################################################################
  title1 <- "Effective rating shifts from the segmentation of gaugings"
  title2 <- "Effective rating shifts from the stage-recession analysis"
  title3 <- "Combined results of segmentation of gaugings and recession analysis"
  title5 <- "Potential rating shifts from the sediment transport proxy analysis"
  title6 <- "Official dates of RC update"
  #preparing datasets of pdf of shift times:
  if (!is.null(data.annotate.gaug.1$MAP)) {
    X1 =pdf.ts.1
    if (ncol(pdf.ts.1)==1) {
      X1.new = X1
    } else {
      X1.new =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts.1))) {
        X1.new =rbind(X1.new, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.rec.1$MAP)) {
    X2 =pdf.ts.2
    if (ncol(pdf.ts.2)==1) {
      X2.new = X2
    } else {
      X2.new =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])))
      for (orderr in 1:(ncol(pdf.ts.2))) {
        X2.new =rbind(X2.new, data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1]))))
      }
    }
  }
  if (!is.null(data.annotate.step1$MAP)) {
    X3 =pdf.ts.3
    if (ncol(pdf.ts.3)==1) {
      X3.new = X3
    } else {
      X3.new =data.frame(time= X3[,1], ord = rep(1, length(X3[,1])))
      for (orderr in 1:(ncol(pdf.ts.3))) {
        X3.new =rbind(X3.new, data.frame(time= X3[,orderr], ord = rep(orderr, length(X3[,1]))))
      }
    }
  }
  
  
  #preparing gauging data per periods:
  #gaugings:
  c_Gaug = 0; P_Gaug = 0;
  if (!is.null(data.annotate.step1$t.adj)) {
    for (i in 1:length(gaug.1$X.tP.)) {
      if(gaug.1$X.tP.[i] <= data.annotate.step1$t.adj[1]) {
        #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
        c_Gaug[i] = colo[1]
        P_Gaug[i] = 1
      }
    }
    for (j in 2:(length(data.annotate.step1$t.adj)+1)) {
      for (i in 1:length(gaug.1$X.tP.)) {
        if ((gaug.1$X.tP.[i] <= tail(gaug.1$X.tP.,1)) & 
            (gaug.1$X.tP.[i] >  data.annotate.step1$t.adj[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          c_Gaug[i] = colo[j]
          P_Gaug[i] = j
        }
      }
    }
  } else {
    for (i in 1:length(gaug.1$X.tP.)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      c_Gaug[i] = colo[1]
      P_Gaug[i] = 1
    }
  }
  df.RC.step1 <- data.frame(gaug.1$X.h., gaug.1$X.Q., gaug.1$X.uQ., P_Gaug, gaug.1$X.tP., c_Gaug)
  names(df.RC.step1) = c("h","Q", "uQ", "Period", "t","color")
  # save into files:
  write.table(df.RC.step1, paste0(dir,"/STEP1_data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(pdf.ts.3, paste0(dir,"/STEP1_pdf_ts.txt"),  sep ="\t", row.names=FALSE)
  write.table(data.annotate.step1, paste0(dir,"/STEP1_shift_times.txt"),  sep ="\t", row.names=FALSE)
  
  
  
  
  
  
  
  #plotting the stage record with gaugings
  #------------------------------------------
  #PLOT 1:
  #------------------------------------------
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray50",size = 0.3)+
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.time)
  } else {
    t.plot= t.plot + 
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  t.plot <- t.plot + 
    geom_point(data=df.RC.step1, aes(x = t , y= h), size = 6, pch =21, fill= df.RC.step1$color)+ #g.1$color)+
    #geom_vline(aes(xintercept = data.annotate.gaug.adjust.1$t.adj), color = "red", lwd =0.3, linetype = "dashed")+
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(start.y.legend, grid_limni.ylim[2]), 
                       breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none",
          plot.margin=unit(c(0.5,0.5,2,0),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.line.x = element_blank())+
    geom_line( aes(x= df.limni$t_limni, 
                   y = hcritical), color = "green", 
               size =1.5, linetype ="solid") +
    annotate("text",
             x= limits.time +20 , 
             y = hcritical[1] + 0.3,  
             label=TeX("$h_c$"), color = "green", size=9, hjust=0) + 
    geom_vline(xintercept = tail(df.limni$t_limni,1), color="red", size= 2.5, linetype ="dashed")+
    annotate("text", 
             x = tail(df.limni$t_limni,1) + 50,
             y = (grid_limni.ylim[2] - start.y.legend)/2,
             label= paste0("Initialization"), angle = 90,
             color = "red", size=10)
  
  
  
  
  #------------------------------------------
  #PLOT 2: gaugings segmentation results
  #------------------------------------------
  if (exists("data.annotate.gaug.1")) {
    t.plot2 <- ggplot() 
    if (ncol(pdf.ts.1)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X1, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X1.new$time, ..scaled.. , group =X1.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
      #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
      # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
      geom_segment(data = data.annotate.gaug.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "blue", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title1)
    if (is.null(df.limni)==FALSE) {
      # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  
  #------------------------------------------
  #PLOT 3: recession segmentation results
  #------------------------------------------
  if (exists("data.annotate.rec.1")) {
    t.plot3 <- ggplot() 
    if (ncol(pdf.ts.2)==1) {
      t.plot3 = t.plot3 + 
        geom_density(aes(x= X2, ..scaled.. ),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot3 = t.plot3 +
        geom_density(aes(x= X2.new$time, ..scaled.. , group =X2.new$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
    t.plot3 = t.plot3 +
      geom_segment(data = data.annotate.rec.1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title2)
    if (is.null(df.limni)==FALSE) {
      # t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  
  ###  NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO !!!!!!!!!!!!
  
  #-----------------------------------------
  #PLOT 4: method B
  #-----------------------------------------
  if ((exists("data.annotate.step1"))) {
    t.plot4 <- ggplot() 
    if (ncol(pdf.ts.3)==1) {
      t.plot4 = t.plot4 + 
        geom_density(aes(x= X3, ..scaled..),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot4 = t.plot4 +
        geom_density(aes(x= X3.new$time, ..scaled.., group =X3.new$ord),  
                     fill= "blue",
                     colour=NA, alpha=0.3)
    }
    t.plot4 = t.plot4 + 
      geom_segment(data = data.annotate.step1, aes(x = t.adj, y = -1, yend =0, xend=t.adj ), size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="blue", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title3)
    if (is.null(df.limni)==FALSE) {
      #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot4= t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  #------------------------------------------
  #PLOT 5: sediment trnapsort
  #------------------------------------------
  if (exists("data.annotate.sed.trasp")) {
    t.plot5 <- ggplot() +
      geom_segment(data = data.annotate.sed.trasp, aes(x = time, y = 0, yend =1, xend=time ),
                   size = 1.5, color ="blue")+
      scale_y_continuous(name = "", expand = c(0,0), limits = c(0,1))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
      theme_set(theme_gray(base_size = 20))+
      theme(axis.text=element_text(size=20)
            ,axis.title=element_text(size=30,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
            ,axis.text.x =element_blank()
            ,axis.text.y =element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_blank()) +
      ggtitle(title5)
      
    if (is.null(df.limni)==FALSE) {
      # t.plot3 = t.plot3 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      t.plot5 = t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
    } else {
      t.plot5 = t.plot5 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
    }
  }
  
  
  
  
  
  #----------------------
  #PLOT 6: official dates
  #----------------------
  t.plot6 <- ggplot() 
  t.plot6 = t.plot6 + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,plot.title = element_text(size = 15, color= "black", vjust=-6)
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2,0.9),"cm")
          ,axis.text.x =element_blank()
          ,axis.text.y =element_blank()
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()) +
    ggtitle(title6)
  if (is.null(df.limni)==FALSE) {
    # t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    t.plot6= t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
  } else {
    t.plot6 = t.plot6 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(g.1$X.tP.,1)))
  }
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null(data.annotate.off)==FALSE) {
      t.plot6 = t.plot6 +
        geom_point(#data = data.annotate.off, 
          aes(x = data.annotate.off$xeffect[which(data.annotate.off$xeffect<= tail(df.limni$t_limni,1))], y = 0),
          color= "black", size = 4, shape =4, stroke=2 )
      # geom_point(data = data.annotate.off, aes(x = xpotent, y = 0), color= "green", size = 4, shape =4, stroke=2 )
    }
  }
  
  
  #-----------------------------------------------------------------------
  #save plot:
  t.plot7 = plot_grid( t.plot2,  
                       t.plot3, 
                       t.plot4,
                       t.plot5,
                       t.plot6, 
                       ncol = 1, nrow = 5, rel_heights = c(1,1,1,1,1)) 
  t.plot8 = plot_grid( t.plot, 
                       t.plot7, 
                       ncol = 1, nrow = 2,  rel_heights = c(1,1))
  # ggsave(t.plot8, filename =paste0(dir,"/STEP1_time_series.png"), 
  #                 device  = "png", width = 20, height =16, dpi = 400, units = "in")
  
  
  
  
  t.plot2bis = t.plot2  + theme(plot.title = element_text(face="bold", size = 25, color= "black"))
  t.plot3bis = t.plot3  + theme(plot.title = element_text(face="bold",size = 25, color= "black"))
  t.plot4bis = t.plot4  + theme(plot.title = element_text(face="bold",size = 25, color= "black"))                              
  t.plot5bis = t.plot5  + theme(plot.title = element_text(face="bold",size = 25, color= "black"))  
  t.plot6bis = t.plot6  + theme(plot.title = element_text(face="bold",size = 25, color= "black")) #, vjust=-6))   
  
  
  t.plot7bis = plot_grid( t.plot2bis,  
                          t.plot3bis, 
                          #t.plot4bis,
                          t.plot5bis,
                          #t.plot6bis, 
                          ncol = 1, nrow = 3, 
                          rel_heights = c(1,1,0.7)) 
  t.plot8bis = plot_grid( t.plot, 
                          t.plot7bis, 
                          ncol = 1,
                          nrow = 2,  
                          rel_heights = c(1.5,0.7))
  pdf(paste0(dir,"/STEP2_time_series.pdf"), 20,14 ,useDingbats=F)
  print(t.plot8bis)
  dev.off()
  
  ggsave(t.plot8bis, filename =paste0(dir,"/STEP2_time_series.png"), 
         device  = "png", width = 20, height =14, dpi = 400, units = "in")
  
  
  
  return(t.plot8)
  #------------------------------------------------------------------------
}









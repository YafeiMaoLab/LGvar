exchange<-function(data){  
  data <- data[data$ref_start != 0, ]
  data[data$ref_start > data$ref_end, c('ref_start', 'ref_end')] <-  data[data$ref_start > data$ref_end, c('ref_end','ref_start')]
  data[data$query_start > data$query_end, c('query_start', 'query_end')] <-  data[data$query_start > data$query_end, c('query_end','query_start')] 
  return(data)
}
##merge SDR
intersect.unit<-function(data,chr_child){
  rm(inte)
  data<-data[data$query_chr==chr_child,]
  data<-exchange(data)
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
  overlapsque<-findOverlaps(ir1que,ir1que) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
  }
  if(length(which(overlapsque@from!=overlapsque@to))!=0){
    if(exists("inte")){
      inte<-rbind(inte,as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)]))
    }
    else{
      inte<-as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)])
    }
  }
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    for(j in 1:dim(inte)[1]){
      a<-c()
      a<-append(a,inte[j,1])
      a<-append(a,inte[j,2])
      for(k in 1:dim(inte)[1]){
        if(inte[k,1] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,1]))
        }
        if(inte[k,2] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,2]))
        }
      }
      if(j==1){
        x<-a
        data[min(x),]$ref_start<-min(data[x,]$ref_start)
        data[min(x),]$ref_end<-max(data[x,]$ref_end)
        data[min(x),]$query_start<-min(data[x,]$query_start)
        data[min(x),]$query_end<-max(data[x,]$query_end)
        data[x,]<-data[min(x),]
        data[x,]$anno<-"COMPLEX"
      }
      else{
        if(length(intersect(a,x))==0){
          x<-a
          data[min(x),]$ref_start<-min(data[x,]$ref_start)
          data[min(x),]$ref_end<-max(data[x,]$ref_end)
          data[min(x),]$query_start<-min(data[x,]$query_start)
          data[min(x),]$query_end<-max(data[x,]$query_end)
          data[x,]<-data[min(x),]
          data[x,]$anno<-"COMPLEX"
        }
      }
    }
    data<-distinct(data)
  }
  data<-data[order(data$ref_start),]
  return(data)
}


##----------------------function1 remove cen and telo alignments
alignintersect<-function(align1,align2,align3){
  del.list<-c()
  cat("-------Start to remove centromere and telomere alignments---------\n")
  
  ir1ref <- IRanges(start = align1$ref_start, end = align1$ref_end)
  ir1que  <- IRanges(start = align1$query_start, end = align1$query_end)
  ir2 <- IRanges(start =align2$ref_start , end =align2$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir2) 
  if(length(overlapsref@from)!=0){  
    xx<-cbind(overlapsref@from,
              align1[overlapsref@from,]$ref_chr,
              align1[overlapsref@from,]$ref_start,
              align1[overlapsref@from,]$ref_end,
              overlapsref@to,
              align2[overlapsref@to,]$ref_chr,
              align2[overlapsref@to,]$ref_start,
              align2[overlapsref@to,]$ref_end)
    xx<-as.data.frame(xx)
    for(i in 1:dim(xx)[1]){
      if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){  
        del.list<-append(del.list,xx$V1[i])
      }
    }
    
  }
  
  ir3 <- IRanges(start = align3$query_start, end = align3$query_end)
  overlapsque<-findOverlaps(ir1que,ir3)
  if(length(overlapsque@from)!=0){  
    xx<-cbind(overlapsque@from,
              align1[overlapsque@from,]$query_chr,
              align1[overlapsque@from,]$query_start,
              align1[overlapsque@from,]$query_end,
              overlapsque@to,
              align3[overlapsque@to,]$query_chr,
              align3[overlapsque@to,]$query_start,
              align3[overlapsque@to,]$query_end)
    xx<-as.data.frame(xx)
    for(i in 1:dim(xx)[1]){
      if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){  
        del.list<-append(del.list,xx$V1[i])
      }
    }
    
  }

  if(length(del.list)!=0)
  {align1<-align1[-as.numeric(del.list),]}
  return(align1)
}

reverse_xy<-function(data){
  selected_rows <- data$orient == '-'
  new_query_start <- data[selected_rows, ]$query_end
  new_query_end <- data[selected_rows, ]$query_start
  data[selected_rows, ]$query_start <- new_query_start
  data[selected_rows, ]$query_end <- new_query_end
  return(data)
}
#-------------------function4 split and cluster ----------------------
split_region<-function(pos.chr.region,cluster.id,clusterparas){
  pos.chr.region<-reverse_xy(pos.chr.region)  #reverse negative align
  df <- pos.chr.region %>%
    rowwise() %>%
    mutate(length = sqrt((ref_end - ref_start)^2 + (query_end- query_start)^2),
           segments = ceiling(length / 5000))
  df_segments <- df[rep(seq_len(nrow(df)), df$segments), ]
  df_segments <- df_segments %>%
    group_by(ref_start, ref_end, query_start, query_end) %>%
    mutate(segment_id = row_number(),
           segment_mid_x = ref_start + (ref_end - ref_start) * (segment_id - 0.5) / segments,
           segment_mid_y = query_start + (query_end - query_start) * (segment_id - 0.5) / segments)
  
  if(cluster.id==2){
    dbs.para<-10000
  }
  if(cluster.id==0){
    dbs.para<-clusterparas
    
  }
  if(cluster.id==1){
    df_segments[df_segments$orient=='-',]$segment_mid_y<--df_segments[df_segments$orient=='-',]$segment_mid_y
    df_segments[df_segments$orient=='-',]$query_start<--df_segments[df_segments$orient=='-',]$query_start
    df_segments[df_segments$orient=='-',]$query_end<--df_segments[df_segments$orient=='-',]$query_end
    pos.chr.region[pos.chr.region$orient=='-',]$query_start<--pos.chr.region[pos.chr.region$orient=='-',]$query_start
    pos.chr.region[pos.chr.region$orient=='-',]$query_end<--pos.chr.region[pos.chr.region$orient=='-',]$query_end
    dbs.para<-clusterparas}
  dbscan_result <- dbscan(df_segments[,c("segment_mid_x", "segment_mid_y")], eps = dbs.para, minPts =1)
  df_segments$cluster <- dbscan_result$cluster
  df_segments <- df_segments %>%
    group_by(ref_chr,ref_start, ref_end, ref_pos,query_chr,query_start, query_end,query_pos,orient) %>%
    summarise(cluster = names(which.max(table(cluster))))
  df_segments$query_start<-abs(df_segments$query_start)
  df_segments$query_end<-abs(df_segments$query_end)
  df_segments<-reverse_xy(df_segments) 
  return(df_segments)
}

#-------------------function5 plot ----------------------
dotplot_cluster <- function(plotpos, region, pdf_file = NULL) {
  plotpos <- reverse_xy(plotpos) 
  
  if (!missing(region)) {
    plotpos <- plotpos[plotpos$ref_end < region[2] & plotpos$ref_start > region[1], ]
  }
  
  has_cluster <- "cluster" %in% colnames(plotpos)
  
  plot <- ggplot(plotpos, aes(x = ref_start, y = query_start)) +
    geom_segment(
      aes(
        xend = ref_end,
        yend = query_end,
        color = if (has_cluster) paste("Cluster", cluster, "(", query_chr, ")") else NULL
      ),
      size = 0.8,
      arrow = arrow(length = unit(0.15, "cm"))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = if (has_cluster) "right" else "none",  
      legend.text = element_text(size = 10)
    ) +
    labs(
      title = "Dot Plot of Segments",
      x = "Reference Position",
      y = "Query Position"
    )
  
  if (has_cluster) {
    plot <- plot + scale_color_viridis_d(
      name = "Cluster (Query Chr)",  
      option = "plasma"
    )
  } else {
    plot <- plot + scale_color_manual(values = "black", guide = "none")
  }
  
  plotly_obj <- ggplotly(plot, tooltip = c("ref_start", "query_start", "ref_end", "query_end", if (has_cluster) "cluster", "query_chr"))
  
  if (!is.null(pdf_file)) {
    ggsave(
      filename = pdf_file,
      plot = plot,
      device = "pdf",
      width = 10,
      height = 7,
      units = "in"
    )
  }
  
  return(plotly_obj)
}


inte.minud<-function(endcluster1,id){
  if(id==1){
    if(length( which(endcluster1$orient=="-"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="-")),]
    }
  }
  if(id==2){
    if(length( which(endcluster1$orient=="+"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="+")),]
    }
  }
  return(endcluster1)
}



docall<-function(data){
  if(length(unique(data$query_chr)[unique(data$query_chr)!='0'])!=1){
    data<-do.call(rbind,mget(unique(data$query_chr)[unique(data$query_chr)!='0'], envir = .GlobalEnv))
  }
  else{
    data<-get( unique(data$query_chr)[unique(data$query_chr)!='0'])
  }
  return(data)
}

##----function extract inversion from cluster2
##!!!try to refine the breakpoints of inversions
optimized_inversions <- list()
inversion.extract<-function(invcluster,chrid){
  invcluster$query_start<-abs(invcluster$query_start)
  invcluster$query_end<-abs(invcluster$query_end)
  inversion<-invcluster[invcluster$orient=="-",] %>% group_by(cluster)%>%  #cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              orient=(names(which.max(table(orient))))) 
  colnames(inversion)[6]<-"query_start"
  colnames(inversion)[7]<-"query_end"
  inversion<-distinct(inversion) 
  #1.record + cluster
  write.table(invcluster[invcluster$orient == "+", c("ref_chr","ref_start","ref_end","ref_pos","query_chr","query_start","query_end","query_pos","orient","cluster")],
            file = "invcluster_forward.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # 2.record inversions
  inversion_reordered <- inversion[, c(2:ncol(inversion), 1)]  #move the first column to the end so we can do overlap

  write.table(inversion_reordered,
            file = "inversion_reverse.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  # 3.find overlap by bedtools intersect
  system("bedtools intersect -a inversion_reverse.bed -b invcluster_forward.bed -wa -wb -f 0.1 > overlap.bed")
  
  
  for (i in seq_len(nrow(inversion))) {
    b <- inversion[i, ]
    
    # 4.whether they have overlaps
    if (file.info("overlap.bed")$size == 0) {
      #breakpoints don't need to be refined, add to dataframe
      optimized_inversions[[length(optimized_inversions) + 1]] <- b
    } 
    else {
      overlaps <- read.table("overlap.bed", header = FALSE, stringsAsFactors = FALSE, comment.char = "")
  
      #find this inversion overlap
      current_overlaps <- overlaps %>%
        filter(overlaps$V1 == b$ref_chr & overlaps$V2 == b$ref_start & overlaps$V3 == b$ref_end)
      current_overlaps <- current_overlaps[order(current_overlaps$V11 - current_overlaps$V10, decreasing = TRUE), ]
      if (nrow(current_overlaps) > 1){
        update_ref_start = 0
        update_ref_end = 0
        update_que_start = 0
        update_que_end = 0
        for (j in seq_len(nrow(current_overlaps))) {
            a <- current_overlaps[j, ]
            #look for whether both breakpoint need to be refined
            if (a$V10 <= b$ref_start) {
              if (a$V11 >= b$ref_end) {
                optimized_inversions[[length(optimized_inversions) + 1]] <- b
                break
                }
              else {
                #(1)inversion_start >= forward align start
                update_ref_start = a$V11 #refine inversion ref_start as plus align ref end
                update_que_start = a$V15 #refine inversion que_start as plus align que end
              }
            } else if (a$V10 >= b$ref_start) {
              #(2)inversion_start <= forward align start
              update_ref_end = a$V10  #refine inversion ref_end as plus align ref start
              update_que_end = a$V14  #refine inversion ref_end as plus align que start
            } 
          }
          refine_ref_start <- ifelse(update_ref_start != 0, update_ref_start, b$ref_start)
          refine_ref_end <- ifelse(update_ref_end != 0, update_ref_end, b$ref_end)
          refine_que_start <- ifelse(update_que_start != 0, update_que_start, b$query_start)
          refine_que_end <- ifelse(update_que_end != 0, update_que_end, b$query_end)

          optimized_inversions[[length(optimized_inversions) + 1]] <- data.frame(
              cluster = b$cluster,
              ref_chr = b$ref_chr,
              ref_start = refine_ref_start,  
              ref_end = refine_ref_end,
              query_chr = b$query_chr,
              query_start = refine_que_start,
              query_end = refine_que_end,
              orient = "-"
            )
      }
      else if (nrow(current_overlaps) == 0){
        optimized_inversions[[length(optimized_inversions) + 1]] <- b
      }
      else {
        a = current_overlaps
        if (a$V10 <= b$ref_start) {
            ##if the inversion is involved in the plus align, it's fake inversion, skip
            if (a$V11 >= b$ref_end) {
                optimized_inversions[[length(optimized_inversions) + 1]] <- b
                }
            else {
            #only need to refine one breakpoint
            optimized_inversions[[length(optimized_inversions) + 1]] <- data.frame(
              cluster = b$cluster,
              ref_chr = b$ref_chr,
              ref_start = a$V11,  ##ref_start=overlap end
              ref_end = b$ref_end,  
              query_chr = b$query_chr,
              query_start = a$V15, ##query_start=overlap end
              query_end = b$query_end,
              orient = "-"
            )
          }
          } else if (a$V10 >= b$ref_start) {
              ##if the inversion contains plus align, don't refine, add it to dataframe
              if (a$V11 <= b$ref_end) {
                optimized_inversions[[length(optimized_inversions) + 1]] <- b
                }
              else {
                #only need to refine one breakpoint
                optimized_inversions[[length(optimized_inversions) + 1]] <- data.frame(
                  cluster = b$cluster,
                  ref_chr = b$ref_chr,
                  ref_start = b$ref_start,  
                  ref_end = a$V10,  ##ref_end=overlap start
                  query_chr = b$query_chr,
                  query_start = b$query_start,
                  query_end = a$V14, ##query_end=overlap start
                  orient = "-"
                )
              }
          }
        } 
      }
    }

  file.remove("invcluster_forward.bed")
  file.remove("inversion_reverse.bed")
  file.remove("overlap.bed")

  #combine results
  optimized_inversions_df <- do.call(rbind, optimized_inversions)
  return(optimized_inversions_df)
}

## merge inversion
clusterbigminus<-function(endcluster1,cluster.id,add){
  rm("xx")
  rm("inte")
  pos_end<-endcluster1 %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) 
  pos_end<-distinct(pos_end)
  pos_end<-pos_end[order(pos_end$ref_start),]

  pos_end<-distinct(pos_end)
  
cor1<-10
  for(chr in unique(endcluster1$query_chr)){
    x1<-pos_end[pos_end$query_chr==chr,]
    if(dim(x1)[1]==1){
      cor1<-1
      next
    }
    f1=1:dim(x1)[1]
    f2=x1$query_start
    if (var(f2) == 0) {
      f2 <- f2 + runif(length(f2), min = 1e-6, max = 1e-5)  
    }
    if (var(f1) == 0) {
      f1 <- f1 + runif(length(f1), min = 1e-6, max = 1e-5)  
    }
    cor1<-min(cor1,cor(f1, f2))
  }
  if(cor1<0.5){ 

  }
  if(cluster.id==0 & cor1>=0.5){
    newpos_end<-pos_end
    for(chr_child in unique(pos_end$query_chr)){
      minus<-which(pos_end$`(names(which.max(table(orient))))`=="-" & pos_end$query_chr==chr_child)
      if(length(minus)!=1 &length(minus)!=0){
        for(i in 2:length(minus)){
          if((minus[i]==minus[i-1]+1 |(minus[i]==minus[i-1]+2)) & pos_end[minus[i],]$query_end<=pos_end[minus[i-1],]$query_start ){ 
            if((minus[i]==minus[i-1]+2)){
              endcluster1[endcluster1$cluster==newpos_end[minus[i-1]+1,]$cluster,]$cluster<-newpos_end[minus[i-1],]$cluster
              newpos_end[minus[i-1]+1,]$cluster<-newpos_end[minus[i-1],]$cluster
            }
            endcluster1[endcluster1$cluster==newpos_end[minus[i],]$cluster,]$cluster<-newpos_end[minus[i-1],]$cluster
            newpos_end[minus[i],]$cluster<-newpos_end[minus[i-1],]$cluster
            colnames(newpos_end)[8]<-"orient"
            xx<-newpos_end %>% group_by(cluster)%>%  
              summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) 
            xx<-distinct(xx)
            
          }
        }
        
      }
    }
  }
  if(exists("xx")){
    return(list(endcluster1=endcluster1,pos_end=xx))
  }else{
    return(list(endcluster1=endcluster1,pos_end=pos_end))
  }
}

#-------------------function6 SDR ----------------------
reverse.region<-function(endcluster1,chrid,cluster.id,add){
  endcluster1$query_start<-abs(endcluster1$query_start)
  endcluster1$query_end<-abs(endcluster1$query_end)

  if(cluster.id==0){

    for(chr_child in unique(endcluster1$query_chr)){
      cal_len<-endcluster1[endcluster1$query_chr==chr_child,]
      cal_len$len<-cal_len$ref_end-cal_len$ref_start
      result <- cal_len %>%
        group_by(cluster) %>%
        summarise(total_length = sum(len))
      thereshod1=1000000
      cross.region<-cross.calcu(endcluster1)
      for(k in setdiff(result[result$total_length<thereshod1,]$cluster,cross.region)){

        if(k==1){
          next
        }

        for(j in result[result$total_length>=thereshod1,]$cluster){
          a<-which(endcluster1$cluster==k & endcluster1$query_chr==chr_child )
          b<-which(endcluster1$cluster==j & endcluster1$query_chr==chr_child)
          if(min(b)<min(a) & max(b)>max(a)){
            cat(k,j,"\n")
            endcluster1[(which(endcluster1$cluster==k & endcluster1$query_chr==chr_child)),]$cluster<-j
            endcluster1<-distinct(endcluster1)
            break
          }
          
        }
        
      }
      
      cal_len<-endcluster1[endcluster1$query_chr==chr_child,]
      cal_len$len<-cal_len$ref_end-cal_len$ref_start
      result <- cal_len %>%
        group_by(cluster) %>%
        summarise(total_length = sum(len))

      rm("delsytenic")
      rm("dupli")
      for(k in result[result$total_length<thereshod1,]$cluster){
        if(k==1){
          next
        }
        if(max(which(endcluster1$cluster==k))==dim(endcluster1)[1]){
          next
        }
        idd<-(which(endcluster1$cluster==k & endcluster1$query_chr==chr_child))
        if(min(idd)==1 |max(idd)==dim(endcluster1)[1]){
          next
        }
        if("-" %in%  unique(endcluster1[idd,]$orient)){next}
        if((endcluster1[min(idd)-1,]$orient=="-" & endcluster1[max(idd)+1,]$orient=="-") |(endcluster1[max(idd)+1,]$orient=="-" & endcluster1[min(idd),]$ref_end>endcluster1[max(idd)+1,]$ref_start)){
          if((endcluster1[min(idd),]$query_end<=endcluster1[min(idd)-1,]$query_start) &(endcluster1[max(idd),]$query_start>=endcluster1[max(idd)+1,]$query_end)){
            if(abs(endcluster1[min(idd),]$query_end-endcluster1[min(idd)-1,]$query_start)<abs(endcluster1[max(idd),]$query_start-endcluster1[max(idd)+1,]$query_end)){
              endcluster1[idd,]$cluster<-endcluster1[min(idd)-1,]$cluster
            }
            else{
              endcluster1[idd,]$cluster<-endcluster1[max(idd)+1,]$cluster
            }
          }

          else if((endcluster1[min(idd),]$query_start>=endcluster1[min(idd)-1,]$query_end) &(endcluster1[max(idd),]$query_end<=endcluster1[max(idd)+1,]$query_start)){
            if(abs(endcluster1[min(idd),]$query_start-endcluster1[min(idd)-1,]$query_end)<abs(endcluster1[max(idd),]$query_end-endcluster1[max(idd)+1,]$query_start)){
              endcluster1[idd,]$cluster<-endcluster1[min(idd)-1,]$cluster
            }
            else{
              endcluster1[idd,]$cluster<-endcluster1[max(idd)+1,]$cluster
            }
          }

          else{
            if(exists("delsytenic")){
              delsytenic<-rbind(delsytenic,endcluster1[idd,])
            }
            else{
              delsytenic<-endcluster1[idd,]
            }
            
            if(exists("dupli")){
              dupli<-rbind(dupli,transduplication_extract(endcluster1,endcluster1[idd,]))
            }else{dupli<-transduplication_extract(endcluster1,endcluster1[idd,])}
            endcluster1<-endcluster1[-idd,]
          }
        }
      }
    }

    
    for (iteration in 1:3) { 
      del.list<-c()
      for(chr_child in unique(endcluster1$query_chr)){
        for(k in unique(endcluster1[endcluster1$query_chr==chr_child,]$cluster)[-1]){
          datalist<-which(endcluster1$cluster==k  & endcluster1$query_chr==chr_child)
          if(endcluster1[min(datalist),]$ref_start<endcluster1[min(datalist)-1,]$ref_end){
            intervalue<-abs(endcluster1[min(datalist),]$ref_start-endcluster1[min(datalist)-1,]$ref_end)
            allval<-endcluster1[min(datalist),]$ref_end-endcluster1[min(datalist),]$ref_start
            allval2<-endcluster1[min(datalist)-1,]$ref_end-endcluster1[min(datalist)-1,]$ref_start
            if(intervalue/allval>0.5){
              del.list<-append(del.list,min(datalist))
            }
            else if(intervalue/allval2>0.5){
              del.list<-append(del.list,min(datalist)-1)
            }
          }
        }
        
      }
      if(length(del.list)!=0){
        if(exists("dupli")){
          dupli<-rbind(dupli,endcluster1[del.list,])
        }else{dupli<-endcluster1[del.list,]}
      }
    }
    if(!exists("dupli")){dupli<-"no"}

    if(cluster.id==1){
      invdata<-endcluster1[endcluster1$orient=="-",]
      for(k in unique(invdata$cluster)){
        if(k==1){
          next
        }
        temdata<-invdata[invdata$cluster==k,]
        mark<-which(endcluster1$cluster==k & endcluster1$orient=="-")
        if(endcluster1[min(mark),]$ref_start<endcluster1[min(mark)-1,]$ref_end | endcluster1[min(mark),]$query_start<endcluster1[min(mark)-1,]$query_end){
          endcluster1[mark,]$cluster<-endcluster1[min(mark)-1,]$cluster
        }
        if(endcluster1[max(mark),]$ref_end>endcluster1[max(mark)+1,]$ref_start | endcluster1[min(mark),]$query_end>endcluster1[min(mark)+1,]$query_start){
          endcluster1[mark,]$cluster<-endcluster1[max(mark)+1,]$cluster
        }
      }
    }
  }
  
  testminus<-clusterbigminus(endcluster1,cluster.id)
  endcluster1<-testminus$endcluster1
  pos_end<-testminus$pos_end
  pos_end<-pos_end[order(pos_end$ref_start),]

  if(nrow(endcluster1)!=1){
    middledata<-duplication_extract(endcluster1,pos_end) 
    pos_end<-middledata$pos_end
    
    if(exists("dupli")){
      if(is.character(dupli)){
        dupli<-middledata$dupli
      }else{
        dupli<-rbind(dupli,middledata$dupli)
      }
    }else{
      dupli<-middledata$dupli
    }
  }
  
  
  reall<-merge(pos_end,chrpc,by=c("ref_chr","query_chr"))
  reall<-reall[order(reall$ref_start),]
  
  for(chr_child in unique(endcluster1$query_chr)){
    chrchch<-reall[reall$query_chr==chr_child,]
    rm("inte")
    chrchch<-chrchch[chrchch$ref_start<=chrchch$ref_end,]
    ir1ref <- IRanges(start = chrchch$ref_start, end = chrchch$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1ref) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      inte$anno<-0
      for(k in 1:dim(inte)[1]){
        if(chrchch[inte[k,1],]$ref_end==chrchch[inte[k,2],]$ref_start){
          inte[k,]$anno<-"no"
        }
        if(chrchch[inte[k,1],]$ref_start>chrchch[inte[k,2],]$ref_start & chrchch[inte[k,1],]$ref_end<chrchch[inte[k,2],]$ref_end){
          inte[k,]$anno<-inte[k,1]
        }
        if(chrchch[inte[k,1],]$ref_start<chrchch[inte[k,2],]$ref_start & chrchch[inte[k,1],]$ref_end>chrchch[inte[k,2],]$ref_end){
          inte[k,]$anno<-inte[k,2]
        }
      }
      inte<-inte[inte$anno!="no",]
      rm("chrduplist")
      if(dim(inte[inte$anno!=0,])[1]!=0){
        chrduplist<-as.numeric(inte[inte$anno!=0,]$anno)
      }
      
      inte<-inte[!inte$V1 %in% inte[inte$anno!=0,]$anno,]
      inte<-inte[!inte$V2 %in% inte[inte$anno!=0,]$anno,]
      inte<-inte[inte$anno==0,]
      
      if(dim(inte)[1]!=0){
        for(k in 1:dim(inte)[1]){
          list<-c(inte[k,1],inte[k,2])
          for(j in 1:dim(inte)[1]){
            if(inte[j,1] %in% list |inte[j,2] %in% list){
              list<-append(list,c(inte[j,1],inte[j,2]))
            }
          }
          if(length(unique(diff(unique(list))))==1 & 1%in% unique(diff(unique(list))) ){
            rows_to_sort<-unique(list)
            sorted_rows <- chrchch[rows_to_sort[order(chrchch$ref_end[rows_to_sort])], ]
            chrchch[rows_to_sort,]<-sorted_rows
          }
        }
      }
      
      
    }
    if(exists("chrduplist")){
      chrchch<-chrchch[-chrduplist,]
    }
    if(dim(chrchch)[1]>1){
      for(m in 2:dim(chrchch)[1]){
        if((chrchch[m,]$`(names(which.max(table(orient))))` != '-') & (chrchch[m-1,]$`(names(which.max(table(orient))))` != '-') & (chrchch[m,]$ref_start < chrchch[m-1,]$ref_end)){
          chrchch[m,]$ref_start <- chrchch[m-1,]$ref_end
        }
      }
    } 
    
    if(length(unique(endcluster1$query_chr))!=1){
      if(min(which(reall$query_chr==chr_child))==1){
        new_row <- data.frame(ref_chr=chrid,ref_start = 0, ref_end = as.numeric(chrchch$ref_start[1]),query_chr=chr_child,query_start = 0, query_end = as.numeric(min(as.numeric(chrchch$query_start)))) 
      }else{
        new_row <- data.frame() 
      }
    }else{
      new_row <- data.frame(ref_chr=chrid,ref_start = 0, ref_end = as.numeric(chrchch$ref_start[1]),query_chr=chr_child,query_start = 0, query_end = as.numeric(min(as.numeric(chrchch$query_start)))) 
    }
    
    if(cluster.id==0){
      
      if(dim(chrchch)[1]>1){
        for (i in 1:dim(chrchch)[1]){
          if(((i != dim(chrchch)[1]) & (chrchch[i,]$`(names(which.max(table(orient))))` != '-') & (chrchch[i+1,]$`(names(which.max(table(orient))))` != '-')) | (i == dim(chrchch)[1])){

          nowvalue=chrchch$query_end[i]
          if(nowvalue==max(reall$query_end)){
            if(i==dim(chrchch)[1]){ 
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
            }
            else{
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
            }
          }else{
            rm("enda")
            query_endlist<-endcluster1[endcluster1$query_start>=nowvalue,]$query_start
            if(i!=dim(chrchch)[1]){
              if(chrchch[i+1,]$query_end<=chrchch[i,]$query_start){
                nowvalue=chrchch$query_start[i]
                query_endlist<-chrchch$query_end[which(chrchch$query_end<=nowvalue)]
                queendval<-query_endlist[which.min(abs(query_endlist -nowvalue))]
                if(length(queendval)!=0){
                  enda<-data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start =queendval, query_end = nowvalue)
                }
                
              }else{
                queendval<-query_endlist[which.min(abs(query_endlist -nowvalue))]
                if(length(queendval)!=0){
                  enda<-data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = queendval)
                }
                
              }
            }
            if(length(query_endlist)==0){
              next
            }
            if(i==dim(reall)[1]){
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(max(chrchch$query_end)), query_end = chrchch$query_len[i])
            }else if(i==dim(chrchch)[1]){
              xx<-data.frame()
            }
            else{
              if(exists("enda")){xx=enda}
              else{next}
              
            }
          }
          if(nrow(xx)!=0){
            colnames(xx)<-c("ref_chr","ref_start","ref_end","query_chr","query_start","query_end")
          }
          new_row <-rbind(new_row,xx)
        }
        }
      }
      else{
        if(dim(pos_end)[1]==1){
          for (i in 1:dim(chrchch)[1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
            colnames(xx)<-colnames(new_row)
            new_row <-rbind(new_row,xx)
          }
        }
        
      }
    }
    if(cluster.id==2){
      if(dim(chrchch)[1]>1){
        for (i in 2:dim(chrchch)[1]-1){
          if(((i != dim(chrchch)[1]-1) & (chrchch[i,]$`(names(which.max(table(orient))))` != '-') & (chrchch[i+1,]$`(names(which.max(table(orient))))` != '-')) | (i == dim(chrchch)[1]-1)){
          if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = as.numeric(chrchch$query_start[i+1]))
          }else if(as.numeric(chrchch$query_end[i])==chrchch$query_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = as.numeric(chrchch$query_start[i+1]))
          }else{
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end =as.numeric(chrchch$query_start[i+1]))
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
        }
      }
    }
    
    if(add=="mius"){
      for(j in 1:dim(chrchch)[1]){
        for( k in setdiff(1:dim(chrchch)[1],j)){
          if(chrchch[k,]$ref_start>=chrchch[j,]$ref_start & chrchch[k,]$ref_end<=chrchch[j,]$ref_end & chrchch[k,]$query_start>=chrchch[j,]$query_start & chrchch[k,]$query_end<=chrchch[j,]$query_end){
            chrchch[j,]$ref_start<-min(chrchch[c(k,j),]$ref_start)
            chrchch[j,]$ref_end<-max(chrchch[c(k,j),]$ref_end)
            chrchch[j,]$query_start<-min(chrchch[c(k,j),]$query_start)
            chrchch[j,]$query_end<-max(chrchch[c(k,j),]$query_end)
            endcluster1[endcluster1$cluster==chrchch[k,]$cluster,]$cluster<-chrchch[j,]$cluster
            chrchch[k,]<-chrchch[j,]  }
        }
      }
    }
    
    chrchch<-distinct(chrchch)
    if(cluster.id==3){
      if(dim(chrchch)[1]>1){
        for (i in 2:dim(chrchch)[1]-1){
          if(((i != dim(chrchch)[1]-1) & (chrchch[i,]$`(names(which.max(table(orient))))` != '-') & (chrchch[i+1,]$`(names(which.max(table(orient))))` != '-')) | (i == dim(chrchch)[1]-1)){
          if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }else if(as.numeric(chrchch$query_start[i])==chrchch$query_end[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }else{
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
        }
        
      }
      
    }
    if(nrow(new_row)!=0){
      colnames(new_row)<-c("ref_chr","ref_start","ref_end","query_chr","query_start","query_end")
      
    }
    assign(chr_child,new_row) 
  }
  allchr.reverse<-do.call(rbind,mget(unique(endcluster1$query_chr)))
  pointer<-which(allchr.reverse$ref_start==1)[which(allchr.reverse$ref_start==1)!=1]
  allchr.reverse[pointer,]$ref_start<-allchr.reverse[pointer-1,]$ref_start
  allchr.reverse[pointer-1,]$ref_end<-allchr.reverse[pointer,]$ref_end
  rm(list=unique(endcluster1$query_chr)) 
  
  
  if(cluster.id==2 | cluster.id==3){
    allchr.reverse<-allchr.reverse[-1,]
  }
  if(nrow(allchr.reverse)!=0){
    allchr.reverse$anno<-"SDR_NM"
  }
  if(cluster.id==0){
    if(exists("delsytenic")){
      return(list(reverse=allchr.reverse,endcluster1=endcluster1,minimaploc=reall,dup=dupli,delsytenic=delsytenic))
    }else{
      return(list(reverse=allchr.reverse,endcluster1=endcluster1,minimaploc=reall,dup=dupli))
    }
    
  }else{
    return(list(reverse=allchr.reverse,endcluster1=endcluster1))
  }
  
}

#-------------------function7 integrate dup ----------------------
repeat.integrate<-function(data,repeatid){
  repeat_region<-data.frame()
  data$query_start<-abs(data$query_start)
  data$query_end<-abs(data$query_end)
  data<-exchange(data)
  for(kk in 1:2){
    ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1ref) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
      if(exists("inte")){
        df_sorted <- t(apply(inte, 1, function(x) sort(x)))
        unique_df <- unique(df_sorted)
        inte <- as.data.frame(unique_df)
        deletelist<-c()
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            if(data[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
              deletelist<-append(deletelist,j)
            }
          }
          if(length(deletelist)!=0){
            inte<-inte[-deletelist,]
          }
        }
        
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            a<-c()
            a<-append(a,inte[j,1])
            a<-append(a,inte[j,2])
            for(k in 1:dim(inte)[1]){
              if(inte[k,1] %in% ((min(a)):(max(a)))){
                a<-unique(append(a,inte[k,1]))
              }
              if(inte[k,2] %in% ((min(a)):(max(a)))){
                a<-unique(append(a,inte[k,2]))
              }
            }
            x<-a
            x<-sort(x,decreasing = TRUE)
            for(mm in 1:(length(x)-1)){
              data[x[mm],]  
              data[x[mm+1],] 
              if(data[x[mm],]$ref_start<data[x[mm+1],]$ref_end & data[x[mm],]$ref_start>=data[x[mm+1],]$ref_start){
                v<-min(data[x[c(mm,mm+1)],]$ref_end)-data[x[mm],]$ref_start
                o<-data[x[mm],]$ref_end-data[x[mm],]$ref_start
                o2<-data[x[mm+1],]$ref_end-data[x[mm+1],]$ref_start
                if(o2==0 | o==0){
                  next
                }
                if(v/o>0.5){
                  repeat_region<-rbind(repeat_region,data[x[mm],])
                }
                if(v/o2>0.5){
                  repeat_region<-rbind(repeat_region,data[x[mm+1],])
                }
                if(dim(data[data$ref_start>=min(data[x[c(mm,mm+1)],]$ref_end) & data$ref_end<=-data[x[mm],]$ref_start,])[1]!=0){
                  repeat_region<-rbind(data[data$ref_start>=min(data[x[c(mm,mm+1)],]$ref_end) & data$ref_end<=-data[x[mm],]$ref_start,],)
                }
                data[x[mm+1],]$ref_end<-data[x[mm+1],]$ref_start
                
              }
            }
          }
        }
        
        
      }
    }
    if(kk==1){
      dup<-repeat_region
    }
    data<-distinct(data)
  }

  if(repeatid==1){
    i=2
    while (i<=dim(data)[1]){
      list<-which(data$query_start[i]<data[1:(i-1),]$query_end)
      
      if(length(list)!=0){
        if(data$orient[i]=="-" & all(data$orient[list] == '-')){
          i<-i+1
          next
        }
        if(length(list)>0.5*dim(data)[1]){
          data<-data[-i,]
          next
        }
        data[i,]$ref_start<-min(data[c(list,i),]$ref_start)
        data[i,]$ref_end<-max(data[c(list,i),]$ref_end)
        data[i,]$query_start<-min(data[c(list,i),]$query_start)
        data[i,]$query_end<-max(data[c(list,i),]$query_end)
        data<-data[-list,]
      }
      else{
        i<-i+1
      }
      
    }
  }
  
  data<-data[order(data$ref_start),]
  data<-distinct(data)
  if(nrow(dup)!=0){
    dup<-dup[c("ref_chr","ref_start","ref_end","query_chr","query_start","query_end","orient","cluster")]
  }
  
  return(list(afterdup=data,repeat.region=dup))
  
}

#-------------------function7 merge results ----------------------
endfilter<-function(all,chrid,chr_child,distance){
  data<-all
  data$reflen<-data$ref_end-data$ref_start
  data$querylen<-data$query_end-data$query_start
  distance <- as.integer(distance)
  if(nrow(data[(data$reflen<10000 & data$querylen<10000) & data$anno=="SDR_NM",])!=0){
    data[(data$reflen<10000 & data$querylen<10000) & data$anno=="SDR_NM",]$anno<-"SV_NM"
  }
  
  if(length(which(data$reflen<=distance & data$querylen>=10000 & data$reflen>=0))!=0){
    data[data$reflen<=distance & data$querylen>=10000 & data$reflen>=0,]$anno<-"SDR_INS"
  }
  if(length(which(data$reflen==0 & data$querylen<10000))!=0){
    data[data$reflen==0 & data$querylen<10000,]$anno<-"SV_INS"
  }
  if(length(which(data$querylen<=distance & data$reflen>=10000 & data$querylen>=0))!=0){
    data[data$querylen<=distance & data$reflen>=10000 & data$querylen>=0,]$anno<-"SDR_DEL"
  }
  if(length(which(data$querylen==0 & data$reflen<10000))!=0){
    data[data$querylen==0 & data$reflen<10000,]$anno<-"SV_DEL"
  }
  
  if(length(which(data$anno=="INV" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="INV" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_INV"
  }
  
  if(length(which(data$anno=="INV" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="INV" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_INV"
  }
  if(length(which(data$anno=="TRANS" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="TRANS" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_TRANS"
  }
  if(length(which(data$anno=="TRANS" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="TRANS" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_TRANS"
  }
  if(length(which(data$anno=="DUP" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="DUP" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_DUP"
  }
  if(length(which(data$anno=="DUP" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="DUP" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_DUP"
  }
  
  return(data)
}

insertsmall<-function(endcluster2before,storesmall,orientid){
  orientlist<-rle(endcluster2before[["orient"]] == orientid)
  if(orientid=="+"){ 
    reverseid=2
  }
  else{
    reverseid=3
  }
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
      insertpos=endcluster2before[start:end,]
      if(nrow(insertpos)==1){
        next
      }
      reverse_end<-reverse.region(insertpos,chrid,reverseid,"init")
      if(nrow(reverse_end$reverse)!=0){
        xx<-reverse_end$reverse
        if(reverseid==2){
          xx$orient<-"+"
        }else{
          xx$orient<-"-"
        }
        
        storesmall<-rbind(storesmall,xx[,colnames(storesmall)])
      }
      
    }
  }
  return(storesmall)
}

smalltransf<-function(endcluster0){
  x=table(endcluster0$cluster)
  thereshod<-quantile(sort(as.vector(x)),0.8)/5
  dellist<-c()
  for(k in names(x)[x<=thereshod]){
    if(k==1){
      next
    }
    if(max(endcluster0[endcluster0$cluster==k,]$ref_end)-min(endcluster0[endcluster0$cluster==k,]$ref_start)>2000000){
      next
    }
    numberslist<-rle(endcluster0$cluster)$values
    biglist<-names(x)[x>thereshod]
    tranloc<-which(numberslist==k)
    bigloc<-which(numberslist %in% biglist)
    positive_positions <- which(min(tranloc)-bigloc > 0)
    aclus<-bigloc[positive_positions][which.min(min(tranloc)-bigloc[positive_positions])]
    positive_positions <- which(bigloc-max(tranloc) > 0)
    bclus<-bigloc[positive_positions][which.min(bigloc[positive_positions]-min(tranloc))]
    s<-numberslist[aclus]
    e<-numberslist[bclus]
    a<-which(endcluster0$cluster==k)
    loc1<-max(intersect(1:min(a), which(endcluster0$cluster==s)))
    loc2<-min(intersect(max(a):dim(endcluster0)[1], which(endcluster0$cluster==e)))
    x1<-abs(endcluster0[min(a),]$query_start-endcluster0[loc1,]$query_start)
    x2<-abs(endcluster0[max(a),]$query_end-endcluster0[loc2,]$query_end)
    x3<-abs(endcluster0[min(a),]$ref_start-endcluster0[loc1,]$ref_end)
    x4<-abs(endcluster0[max(a),]$ref_end-endcluster0[loc2,]$ref_start)
    if(unique(is.na(x2)) |length(x2)!=1){
      if(x1/x3>20 ){
        dellist<-append(dellist,k)
      }
    }
    if( unique(is.na(x1)) |length(x1)!=1){
      if(x2/x4>20){
        dellist<-append(dellist,k)
      }
    }
    if(unique(!is.na(x1)) & unique(!is.na(x2)) &length(x2)==1  &length(x1)==1){
      if(x1/x3>20 |x2/x4>20){
        dellist<-append(dellist,k)
      }
    }
    
  }
  if(length(dellist)!=0){
    tran<-endcluster0[endcluster0$cluster %in% dellist,]
    tranbefore<-tran
    tran<-tran %>% group_by(cluster)%>%  
      summarise(ref_chr=ref_chr,ref_start=min(ref_start),
                ref_end=max(ref_end),
                query_chr=query_chr,
                query_starttem=min(pmin(query_end,query_start)),
                query_endtem=max(pmax(query_start,query_end)),
                (names(which.max(table(orient))))) 
    colnames(tran)[6]<-"query_start"
    colnames(tran)[7]<-"query_end"
    endcluster0<-endcluster0[!endcluster0$cluster %in% dellist,]
    colnames(tran)[8]<-"orient"
    return(list(tran=distinct(tran),endcluster0=endcluster0,tranbefore=tranbefore))
  }
  else{
    return(list(tran=NULL,endcluster0=endcluster0))
  }
  
}
clusterall<-function(data){
  data<-data %>% group_by(cluster)%>%  
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              (names(which.max(table(orient))))) 
  colnames(data)[6]<-"query_start"
  colnames(data)[7]<-"query_end"
  colnames(data)[8]<-"orient"
  data<-distinct(data)
  return(data)
}

smallcluster<-function(data,orientid){
  orientlist<-rle(endcluster2[["orient"]] == orientid)
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
    }
    else{
      start<-1
      end<-orientlist$lengths[i]
    }
    data[start,]$ref_start<-min(data[start:end,]$ref_start)
    data[start,]$ref_end<-max(data[start:end,]$ref_end)
    data[start,]$query_start<-min(data[start:end,]$query_start)
    data[start,]$query_end<-max(data[start:end,]$query_end)
    data[start,]$ref_pos<-(data[start,]$ref_start+data[start,]$ref_end)/2
    data[start,]$query_pos<-(data[start,]$query_start+data[start,]$query_end)/2
    data[start:end,]<-data[start,]
  }
  return(distinct(data))
}

complex<-function(data){
  if(dim(data)[1]>1){
    for(i in 2:dim(data)[1]){
      if(data[i,]$ref_start<data[i-1,]$ref_end){
        data[i,]$ref_start<-min(data[i,]$ref_start,data[i-1,]$ref_start)
        data[i,]$ref_end<-max(data[i,]$ref_end,data[i-1,]$ref_end)
        data[i,]$query_start<-min(data[i,]$query_start,data[i-1,]$query_start)
        data[i,]$query_end<-max(data[i,]$query_end,data[i-1,]$query_end)
        data[c(i,i-1),]<-data[i,]
        data[c(i,i-1),]$anno<-'COMPLEX'
      }
    }
  }
  data<-distinct(data)
  if(dim(data)[1]>1){
    for(i in 2:dim(data)[1]){
      if(data[i,]$query_start<data[i-1,]$query_end){
        data[i,]$ref_start<-min(data[i,]$ref_start,data[i-1,]$ref_start)
        data[i,]$ref_end<-max(data[i,]$ref_end,data[i-1,]$ref_end)
        data[i,]$query_start<-min(data[i,]$query_start,data[i-1,]$query_start)
        data[i,]$query_end<-max(data[i,]$query_end,data[i-1,]$query_end)
        data[c(i,i-1),]<-data[i,]
        data[c(i,i-1),]$anno<-'COMPLEX'
      }
    }
  }
  
  data<-distinct(data)
  list<-rle(data$anno)
  for (j in which(list$values=="COMPLEX")){
    if(j!=1){
      start<-sum(list$lengths[1:(j-1)])+1
      end<-start+list$lengths[j]-1
      
    }
    else{
      start<-1
      end<-list$lengths[j]
    }
    data[start,]$ref_start<-min(data[start:end,]$ref_start)
    data[start,]$ref_end<-max(data[start:end,]$ref_end)
    data[start,]$query_start<-min(data[start:end,]$query_start)
    data[start,]$query_end<-max(data[start:end,]$query_end)
    data[start:end,]<-data[start,]
  }
  return(distinct(data))
}
minus.next<-function(endcluster1){
  m=0
  for(clusid in unique(endcluster1[endcluster1$orient=="-",]$cluster)){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,200000)
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster$cluster<-paste(miuscluster$cluster, "1",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    m=m+1
  }
  for(clusid in unique(endcluster1[endcluster1$orient=="+",]$cluster)){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,200000)
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster$cluster<-paste(miuscluster$cluster, "1",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    m=m+1
  }
  endcluster1<-endcluster1[order(endcluster1$ref_start),]
  return(endcluster1)
}
cross.calcu<-function(endcluster1){
  del<-c()
  list<-unique(endcluster1$cluster)
  for(i in 2:(length(unique(endcluster1$cluster))-1)){
    a<-which(endcluster1$cluster==list[i])
    b<-which(endcluster1$cluster==list[i-1])
    if(max(b)>max(a) & min(b)>min(a) &  max(a)>min(b)){
      del<-append(del,list[i])
      del<-append(del,list[i-1])
    }
    c<-which(endcluster1$cluster==list[i+1])
    if(max(c)>max(a)& min(c)>min(a) & max(a)>min(c)){
      del<-append(del,list[i])
      del<-append(del,list[i+1])
    }
    
  }
  return(unique(del))
}
duplication_extract<-function(endcluster1,pos_end){
  data<-pos_end
  rm("inte")
  data<-data[data$ref_start<=data$ref_end,]
  endcluster1<-endcluster1[endcluster1$ref_start<=endcluster1$ref_end,]
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
  }
  duplilist<-data.frame()
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    if(nrow(inte)!=0){
      for(j in 1:dim(inte)[1]){
        
        if(data[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
          
          next
        }
        startregion<-max(data[inte[j,2],]$ref_start,data[inte[j,1],]$ref_start)
        endregion<-min(data[inte[j,2],]$ref_end,data[inte[j,1],]$ref_end)
        smalldup<-endcluster1[endcluster1$ref_start>=startregion & endcluster1$ref_end<=endregion,]
        
        if(nrow(smalldup)!=0){
          duplilist<-rbind(duplilist,smalldup) 
        }else{
          data[inte[j,2],]$ref_start=data[inte[j,1],]$ref_end
        }
        if(startregion>endregion){
          next
        }
        ir1ref <- IRanges(start = endcluster1$ref_start, end = endcluster1$ref_end)
        ir1que <- IRanges(start=startregion,end=endregion)
        overlapsref<-findOverlaps(ir1ref,ir1que) 
        for (m in overlapsref@from) {
          if(endcluster1[m,]$ref_end-endcluster1[m,]$ref_start==0){
            next
          }
          if((endregion-startregion)/(endcluster1[m,]$ref_end-endcluster1[m,]$ref_start)>0.5){
            duplilist<-rbind(duplilist,endcluster1[m,]) 
          }
        }
        
      }
    }
    
    duplilist<-distinct(duplilist)
    
  }
  
  return(list(pos_end=data,dupli=duplilist))}


complexinte<-function(data){
  data$cluster<-1:dim(data)[1]
  for(kk in 1:5){
    rm("inte")
    ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1ref) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      if(nrow(inte)!=0){
        for(j in 1:dim(inte)[1]){
          data[inte[j,2],]$cluster=data[inte[j,1],]$cluster
        }
        data<-clusterall(data)
      }
    }
    
    
  }
  return(data)
}

transduplication_extract<-function(endcluster1,pos_end){
  data<-pos_end
  duplilist<-data.frame()
  for(j in 1:dim(data)[1]){
    x<-endcluster1[which(endcluster1$ref_start>=data[j,]$ref_start & endcluster1$ref_end<=data[j,]$ref_end),]
    duplilist<-rbind(duplilist,duplication_extract(x,x)$dupli)
    endcluster1<-endcluster1[-which(endcluster1$ref_start>=data[j,]$ref_start & endcluster1$ref_end<=data[j,]$ref_end),]
  }
  
  rm("inte")
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  ir1req <- IRanges(start = endcluster1$ref_start, end = endcluster1$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir1req) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
  }
  
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    if(nrow(inte)!=0){
      for(j in 1:dim(inte)[1]){
        if(endcluster1[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
          next
        }
        startregion<-max(endcluster1[inte[j,2],]$ref_start,data[inte[j,1],]$ref_start)
        endregion<-min(endcluster1[inte[j,2],]$ref_end,data[inte[j,1],]$ref_end)
        smalldup<-endcluster1[endcluster1$ref_start>=startregion & endcluster1$ref_end<=endregion,]
        if(nrow(smalldup)!=0){
          duplilist<-rbind(duplilist,smalldup) 
        }
        
        if((endregion-startregion)/(data[inte[j,1],]$ref_end-data[inte[j,1],]$ref_start)>0.5){
          duplilist<-rbind(duplilist,data[inte[j,1],]) 
        }
        if((endregion-startregion)/(endcluster1[inte[j,2],]$ref_end-endcluster1[inte[j,2],]$ref_start)>0.5){
          add<-endcluster1[inte[j,2],colnames(duplilist)]
          duplilist<-rbind(duplilist,add) 
        }
        
      }
    }
    
    duplilist<-distinct(duplilist)
    
  }
  
  return(duplilist)}


complexinte<-function(data){
  data$cluster<-1:dim(data)[1]
  for(kk in 1:5){
    rm("inte")
    ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1ref) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      if(nrow(inte)!=0){
        for(j in 1:dim(inte)[1]){
          data[inte[j,2],]$cluster=data[inte[j,1],]$cluster
        }
        data<-clusterall(data)
      }
    }
    
    
  }
  return(data)
}

highdupregion<-function(pos.chr){
  pos.chr$ref_len<-pos.chr$ref_end-pos.chr$ref_start
  pos.chr$que_len<-pos.chr$query_end-pos.chr$query_start
  ir1ref <- IRanges(start = pos.chr$ref_start, end = pos.chr$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref)
  inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
  if(nrow(inte)!=0){
    inte$intersec=0
    inte$intersec <- apply(inte, 1, function(row) {
      q_len <- pos.chr$ref_len[row["queryHits"]]
      s_len <- pos.chr$ref_len[row["subjectHits"]]
      intersec <- (min(pos.chr[row["queryHits"], 3], pos.chr[row["subjectHits"], 3]) -
                     max(pos.chr[row["queryHits"], 2], pos.chr[row["subjectHits"], 2])) /
        max(q_len, s_len)
      return(intersec)
    })
    inte<-inte[inte$intersec>=0.1,]
    grouped_vectors <- tapply( inte$queryHits, inte$subjectHits, FUN = function(x) x)
    if(type(grouped_vectors)=="list"){
      for(group in names(grouped_vectors)){
        group_list<-grouped_vectors[[group]]
        for(group_id in group_list ){
          grouped_vectors[[group]]<-append(grouped_vectors[[group]],grouped_vectors[[as.character(group_id)]])
        }
        grouped_vectors[[group]]<-unique(grouped_vectors[[group]])
      }
      unique_list <- list()
      for (name in names(grouped_vectors)) {
        value <- grouped_vectors[[name]]

        if(!any(sapply(unique_list, function(x) all(value %in% x))) ){
          unique_list[[name]]<-value
        }
      }
      keep_indices <- c()
      for (i in seq_along(unique_list)) {
        is_subset <- FALSE
        for (j in seq_along(unique_list)) {
          if (i != j && all(unique_list[[i]] %in% unique_list[[j]])) {
            is_subset <- TRUE
            break
          }
        }
        if (!is_subset) {
          keep_indices <- c(keep_indices, i)
        }
      }
      result_sets <- unique_list[keep_indices]
      return(result_sets)
    }else{
      return("no")
    }
  }
  else{
    return("nono")
  }
  
}

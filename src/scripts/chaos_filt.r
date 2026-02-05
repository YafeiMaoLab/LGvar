options(warn = -1)

suppressPackageStartupMessages({library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
})

suppressMessages({
args <- commandArgs(trailingOnly = TRUE)
print(args)
alignintersect<-function(align1,align2,align3){
  del.list<-c()
  cat("init",dim(align1)[1],"\n")
  ir1ref <- IRanges(start = align1$V8, end = align1$V9)
  ir1que  <- IRanges(start = align1$V3, end = align1$V4)
  ir2 <- IRanges(start =align2$ref_start , end =align2$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir2)
  if(length(overlapsref@from)!=0){  
    del.list<-append(del.list,overlapsref@from)
  }
  if(!missing(align3)){
    ir3 <- IRanges(start = align3$query_start, end = align3$query_end)
    overlapsque<-findOverlaps(ir1que,ir3)
    if(length(overlapsque@from)!=0){
      del.list<-append(del.list,overlapsque@from)
    }
  }
  if(length(del.list)!=0)
  {align1<-align1[-del.list,]}
  cat("Removed：",dim(align1)[1],"\n")
  return(align1)
}

flit2saffire<-function(pos,name){
  col_names <- colnames(pos)
  col_names[c(1,3, 4, 8, 9,6,5)] <- c("query_chr","query_start", "query_end", "ref_start", "ref_end","ref_chr","orient")
  colnames(pos) <- col_names
  chrpc<-fread(args[1])
  chrpc<-distinct(chrpc)
  colnames(chrpc)<-c("ref_chr","ref_len","query_chr","query_len")
  saffire<-merge(pos,chrpc,by=c("ref_chr","query_chr"))
  saffire<-saffire[,c("ref_chr","ref_start","ref_end","ref_len","orient","query_chr","query_start","query_end","query_len")]
  saffire_name<-c("#reference_name","reference_start","reference_end","reference_length","strand",  "query_name","query_start","query_end","query_length","perID_by_matches","perID_by_events", "perID_by_all","matches","mismatches","deletion_events", "insertion_events", "deletions","insertions")
  saffire[,c("perID_by_matches","perID_by_events", "perID_by_all","matches","mismatches","deletion_events", "insertion_events", "deletions","insertions")]<-0
  colnames(saffire)<-saffire_name
  saffire$matches<-100
  saffire$mismatches<-100
  write.table(saffire, file = paste(args[2],name,sep = ""), sep = "\t", row.names = FALSE,quote = FALSE)
}


split_region<-function(pos.chr.region){
  col_names <- colnames(pos.chr.region)
  col_names[c(3, 4, 8, 9,5)] <- c("query_start", "query_end", "ref_start", "ref_end","orient")
  colnames(pos.chr.region) <- col_names

  x<-pos.chr.region[pos.chr.region$orient=="-",]
  if(nrow(x)!=0){
    
      selected_rows <- pos.chr.region$orient == '-'
      
      new_query_start <- pos.chr.region[selected_rows, ]$query_end
      new_query_end <- pos.chr.region[selected_rows, ]$query_start
      pos.chr.region[selected_rows, ]$query_start <- new_query_start
      pos.chr.region[selected_rows, ]$query_end <- new_query_end

  }
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
  
  dbscan_result <- dbscan(df_segments[,c("segment_mid_x", "segment_mid_y")], eps = as.numeric(args[5]), minPts =1)
  df_segments$cluster <- dbscan_result$cluster
  df_segments <- df_segments %>%
    group_by(ref_start, ref_end, query_start, query_end,sourse) %>%
    summarise(cluster = names(which.max(table(cluster))))
  merge<-df_segments %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_start=min(ref_start),ref_end=max(ref_end),query_start=min(query_start),query_end=max(query_end))
    #   all<-distinct(all)
  
  
  x<-table(df_segments$cluster)
  del<-as.numeric(names(x[which(x<=5)]))
  for(j in del){
    if(abs(merge[merge$cluster==j,]$ref_end-merge[merge$cluster==j,]$ref_start)>as.numeric(args[6])){
      del<-del[del!=j]
    }
  }
  delall<-df_segments[df_segments$cluster %in% del,]$sourse
  return(delall)
}
print(args[7])
if(args[7]=="cts"){
  cts=TRUE
# cta=FALSE
ctn=FALSE
}
if(args[7]=="ctn"){
  cts=FALSE
# cta=FALSE
ctn=TRUE
}


if(cts){  ##Only REF
  ## telomere
  result_reftelo<-fread(args[9])
  result_quetelo<-""
  ## centromere
  result_refcentr<-fread(args[8])
  result_quecentr<-""
}

if(ctn){
  result_reftelo<-""
  result_quetelo<-""
  result_refcentr<-""
  result_quecentr<-""
}
print(args[7])
pos<-fread(args[3],fill=TRUE )
chrnames <- unique(pos$V6)
numeric_part <- as.numeric(gsub("\\D", "", chrnames))
sorted_chrnames <- chrnames[order(numeric_part)]
for(chrid in sorted_chrnames){
  pos.chr<-pos[pos$V6==chrid,]
  print(chrid)
  align1<-pos.chr 
  if(cts){
    align2<-result_refcentr[result_refcentr$chr==chrid,]         
    intersect<-which(pos.chr$V8>=align2$ref_start & pos.chr$V9<=align2$ref_end)
    if(length(intersect)!=0){
      print("have")
      pos.chr<-pos.chr[-intersect,]
    }
    align2<-result_reftelo[result_reftelo$chr==chrid,]         
    for(j in dim(align2)[1]){
      intersect<-which(pos.chr$V8>=align2[j,]$ref_start & pos.chr$V9<=align2[j,]$ref_end)
      if(length(intersect)!=0){
      pos.chr<-pos.chr[-intersect,]
    }
    }
    assign(chrid,pos.chr)
  }

  
  if(ctn){
    assign(chrid,align1)
  }
}


pos<-do.call(rbind,mget(chrnames))
pos$sourse<-1:dim(pos)[1]
rm(list=chrnames)
flit2saffire(pos,"original.saffire") 

delall<-c()
for(chrid in sorted_chrnames){
  pos.chr<-pos[pos$V6==chrid,]
  if(nrow(pos.chr)!=0){del<-split_region(pos.chr)
  delall<-append(delall,del)}
  
}
xxx<-pos[!pos$sourse %in% delall,]

cat("After chaos：",dim(xxx)[1],"\n")
flit2saffire(xxx,"final.saffire")


write.table(xxx, file = args[4], sep = "\t", row.names = FALSE, col.names = FALSE,quote=FALSE)
print("Complete!!!!")
})

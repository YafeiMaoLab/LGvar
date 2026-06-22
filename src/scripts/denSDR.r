options(warn = -1)

suppressPackageStartupMessages({library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
})

suppressMessages({
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

pos <- read.delim(args[2], header = FALSE, col.names = c("ref_chr", "ref_start", "ref_end", "ref_pos", "query_chr", "query_start", "query_end", "query_pos", "orient"))
chrpc<-fread(args[3])
chrpc<-distinct(chrpc)
colnames(chrpc)<-c("ref_chr","ref_len","query_chr","query_len")

chrnames <- unique(pos$ref_chr)
numeric_part <- as.numeric(gsub("\\D", "", chrnames))
sorted_chrnames <- chrnames[order(numeric_part)] 

# cluster0paras<-700000  
# cluster1paras<-700000
cluster0paras<-args[5]  
cluster1paras<-args[5]
invparas <- args[6]

for(chrid in sorted_chrnames){

  store<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  storesmall<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss',orient="+")
  duplication<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster="0")
  inversion<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0)
  INV_INV<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0)
  smalltrans<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster="0")
  pos.chr<-pos[pos$ref_chr==chrid,]
  pos.chr<-pos.chr[order(pos.chr$ref_start),]
  highdupall<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0)
  storehighdup_sdr<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  
  list<-c()
  dupset<-highdupregion(pos.chr)
  if(length(dupset)!=1 ||length(dupset[[1]]>=5)){
    for (dupsetname in names(dupset)){
      if(length(dupset[[dupsetname]])>=5){
        vectore<-append(dupset[[dupsetname]],as.numeric(dupsetname))
        higndup<-pos.chr[vectore,]
        list<-append(list,vectore)
        higndup$cluster=1
        duplication<-rbind(duplication,higndup[,colnames(duplication)])
        highdupall<-rbind(highdupall,c(chrid,min(higndup$ref_start),max(higndup$ref_end)))
      
        dup_boundary<-pos.chr[min(vectore),]
        dup_boundary$anno="SDR_NM"
        if(min(vectore)!=1){
          dup_boundary$ref_end<-dup_boundary$ref_start
          dup_boundary$ref_start<-pos.chr[min(vectore)-1,]$ref_end
          if(pos.chr[min(vectore)-1,"orient"]=="-" & pos.chr[min(vectore),"orient"]=="-"){
            dup_boundary$query_start<-dup_boundary$query_end
            if(nrow(pos.chr[pos.chr$query_start>=dup_boundary$query_end,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_end,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
            
          }else{
            dup_boundary$query_end<-dup_boundary$query_start
            if(nrow(pos.chr[pos.chr$query_end<=dup_boundary$query_start,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_start,]$query_end
            dup_boundary$query_start<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            ##
            dup_boundary$query_start<-pos.chr[min(vectore)-1,]$query_end
            if(nrow(pos.chr[pos.chr$query_start>=dup_boundary$query_start,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_start,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
          }
          
        }
          
          dup_boundary<-pos.chr[max(vectore),]
          dup_boundary$anno="SDR_NM"
        if(max(vectore)!=dim(pos.chr)[1]){
          dup_boundary$ref_start<-dup_boundary$ref_end
          dup_boundary$ref_end<-pos.chr[max(vectore)+1,]$ref_start
          
          if(pos.chr[max(vectore)+1,"orient"]=="-" & pos.chr[max(vectore),"orient"]=="-"){
            dup_boundary$query_end<-dup_boundary$query_start
            if(nrow(pos.chr[pos.chr$query_end<=dup_boundary$query_start,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_start,]$query_end
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            

          }else{
            dup_boundary$query_start<-dup_boundary$query_end
            if(nrow(pos.chr[pos.chr$query_start>=dup_boundary$query_end,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_end,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
            ##
            dup_boundary$query_end<-pos.chr[max(vectore)+1,]$query_start
            if(nrow(pos.chr[pos.chr$query_end<=dup_boundary$query_end,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_end,]$query_end
            dup_boundary$query_start<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
          }
          
          
        }
        
        }
    }
    if(length(list)!=0){
      pos.chr<-pos.chr[-list,]
    }
  }
  storehighdup_sdr<-distinct(storehighdup_sdr)
  ## 1.Cluster，remove translocation in sytenic region
  endcluster0<-split_region(pos.chr,0,cluster0paras)
  
  aftertrans<-detect_translocations(endcluster0)

  for(i in unique((aftertrans$tranbefore)$cluster)){
    a<-aftertrans$tranbefore[(aftertrans$tranbefore)$cluster==i,]
    a$cluster<-as.character(1:dim(a)[1])
    if(sum(a$orient == "-") > 0.5* nrow(a)){
      reverse_end<-reverse.region(a,chrid,3,"init")
    }
    if(sum(a$orient == "+") > 0.5 * nrow(a) |sum(a$orient == "-") == 0.5 * nrow(a)){
      reverse_end<-reverse.region(a,chrid,2,"init")
    }
    if(exists("reverse_end")){
      middle<-reverse_end$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr)) 
      }
      if(sum(a$orient == "+") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"+"
      }
      if(sum(a$orient == "-") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"-"
      }
      if(sum(a$orient == "-") == 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,middle)
      minusmiddle<-reverse_end$reverse[reverse_end$reverse$ref_start<=reverse_end$reverse$ref_end & reverse_end$reverse$query_start<=reverse_end$reverse$query_end,]
      if(nrow(minusmiddle)!=0){
        minusmiddle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,minusmiddle)
      storesmall<-distinct(storesmall)
    }
  }
  if(length(aftertrans$tran)!=0){
    smalltrans<-rbind(smalltrans,aftertrans$tran)
    smalltrans <- smalltrans[smalltrans$orient=="+",]
  }
  if(!is.null(aftertrans$tran) && !isEmpty(aftertrans$tran)){
    inversion<-rbind(inversion,aftertrans$tran[aftertrans$tran$orient=="-",colnames(inversion)])
    transdup<-transduplication_extract(endcluster0,aftertrans$tran)
    if(nrow(transdup)!=0){
      duplication<-rbind(duplication,transdup[,colnames(duplication)])
    }
    
 }
  
  endcluster0<-aftertrans$endcluster0

  ## 2.extract inversion

  endcluster1<-split_region(endcluster0,1,cluster1paras)
  #invcluster <- endcluster1
  invcluster <- split_region(endcluster0,1,invparas)

  endcluster1<-minus.next(endcluster1)  ##this step??
  
  i=0
  for(miusclu in unique(endcluster1[endcluster1$orient=='-',]$cluster)){
    list<-which(endcluster1$cluster==miusclu)
    if(length(list)>=2){
      for(line in 2:length(list)){
        if(endcluster1[list[line],]$ref_start-endcluster1[list[line-1],]$ref_end >50000 && endcluster1[list[line-1],]$query_end <= endcluster1[list[line],]$query_end){
          endcluster1[list[line],]$cluster<-paste(endcluster1[list[line-1],]$cluster,i,sep="")
          i<-i+1
        }else{
          endcluster1[list[line],]$cluster<-endcluster1[list[line-1],]$cluster
        }
      }
    }
   
    endcluster1[endcluster1$cluster==miusclu,]
  }
  
  i=0
  for(miusclu in unique(invcluster[invcluster$orient=='-',]$cluster)){
    list<-which(invcluster$cluster==miusclu)
    if(length(list)>=2){
      for(line in 2:length(list)){
        if(invcluster[list[line],]$ref_start-invcluster[list[line-1],]$ref_end >50000 && invcluster[list[line-1],]$query_end <= invcluster[list[line],]$query_end){
          invcluster[list[line],]$cluster<-paste(invcluster[list[line-1],]$cluster,i,sep="")
          i<-i+1
        }else{
          invcluster[list[line],]$cluster<-invcluster[list[line-1],]$cluster
        }
      }
    }
    
    invcluster[invcluster$cluster==miusclu,]
  }

  INV_INV <- call_nested_inversion(invcluster,chrid)

  refine_breakpoints <- as.logical(args[8]) 
  inver <- inversion.extract(invcluster, chrid, refine_breakpoints = refine_breakpoints)
  if (!is.null(inver) && nrow(inver) != 0) {
    colnames(inver)[8]="orient"
  } else {
    print("Warning: 'inver' is either NULL or has no rows.")
  }
  
  if(nrow(duplication)!=1){
    duplication$cluster<-as.character(duplication$cluster)
    inver<-rbind(inver,duplication[duplication$orient=="-",colnames(inver)])
  }

  list<-which(rle(rle(endcluster1$cluster)$lengths)$values==1 &rle(rle(endcluster1$cluster)$lengths)$lengths >5)
  chaoval<-rle(rle(endcluster1$cluster)$lengths)$values
  chaolen<-rle(rle(endcluster1$cluster)$lengths)$lengths
  cross.region<-cross.calcu(endcluster1)
  initclsuer<-1000
  for(k in list){
    start<-sum(chaoval[1:k-1]*chaolen[1:k-1])+1
    merclus=endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]$cluster
    chaosdf<-endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]
    for(quevalsub in unique(chaosdf$query_chr)){
      endcluster1[endcluster1$query_chr==quevalsub & endcluster1$cluster %in% merclus,]$cluster<-as.character(initclsuer)
      chaosdfsub=chaosdf[chaosdf$query_chr==quevalsub & chaosdf$cluster %in% merclus,]
      chaosdfsub$query_start<-abs(chaosdfsub$query_start)
      chaosdfsub$query_end<-abs(chaosdfsub$query_end)
      if(nrow(chaosdfsub)!=0){
        chaosdfsub<-clusterall(chaosdfsub)
      }
      chaosdfsub$anno<-"COMPLEX"
      store<-rbind(store,chaosdfsub[,colnames(store)])
      initclsuer<-initclsuer+1
      }
  }
  
  reverse_end<-reverse.region(endcluster1,chrid,0,"init")
  duplic<-reverse_end$dup
  if(!is.character(duplic)){
    if(nrow(duplic)!=0){
      duplication<-rbind(duplication,duplic[,colnames(duplication)])
    }
    
  }
  minimap<-reverse_end$minimaploc 
  
  
  if ("delsytenic" %in% names(reverse_end)){
    for(cluster in unique(reverse_end$delsytenic$cluster)){
      datalist<-reverse_end$delsytenic[reverse_end$delsytenic$cluster==cluster,]
      datalist$cluster<-1:dim(datalist)[1]
      reverse_end1<-reverse.region(datalist,unique(datalist$ref_chr),2,"init")
      middle<-reverse_end1$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr)) 
        middle$orient<-"no"
        storesmall<-rbind(storesmall,middle)
      }
    }
  }
  
  for (value in unique(reverse_end$reverse$query_chr)) {
    datarev<-reverse_end$reverse[reverse_end$reverse$query_chr==value,]
    if(nrow(datarev)!=0){
      rows_to_remove <- which(datarev$query_chr == value)

      if (length(rows_to_remove) > 1) {
        startend<-datarev[c(rows_to_remove[1], tail(rows_to_remove, 1)), ]
        store<-rbind(store,datarev[-c(rows_to_remove[1], tail(rows_to_remove, 1)), ]) 
        startend$orient<-"+"
        storesmall<-rbind(storesmall,startend)
      }
    }
  }

  vectore<-which(store$ref_start-store$ref_end==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$ref_start<-store[l,]$ref_start-1
      store[l,]$ref_end<-store[l,]$ref_end+1
    }
  }
  vectore<-which((store$query_start-store$query_end)==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$query_start<-store[l,]$query_start-1
      store[l,]$query_end<-store[l,]$query_end+1
    }
  }
  
  if(nrow(store[store$ref_start>store$ref_end,])!=0){
    store[store$ref_start>store$ref_end,]$ref_end=store[store$ref_start>store$ref_end,]$ref_start 
  }
  initstore<-store[store$ref_chr!=0,]
  if(length(unique(store$query_chr)[-1])!=0){
    for(chr_child in unique(store$query_chr)[-1]){
      new_row<-intersect.unit(store,chr_child) 
      if(nrow(new_row)!=1){
        new_row<-complex(new_row)
      }
      percen<-max(new_row[new_row$anno=="COMPLEX",]$ref_end-new_row[new_row$anno=="COMPLEX",]$ref_start)
      if(percen>100000000){ 
        new_row<-intersect.unit(store,chr_child) 
      }
      assign(chr_child,new_row)
    }
    store<-docall(store)
    rm(list=unique(store$query_chr)) 
  }
  store<-rbind(store,anti_join(initstore,store))
  
 if(!is.null(inver) && dim(inver)[1]!=0){
    inver<-inver[,colnames(inversion)]
    inversion<-rbind(inversion,inver)
    inversion<-distinct(inversion)
  }
  
  cat("Big cluster among:cluster num",dim(store)[1],"\n")  
  cat("inversion:",dim(inversion)[1]-1,"\n") 
  cat("translocation:",dim(smalltrans)[1]-1,"\n") 

  endcluster1<-reverse_end$endcluster1
  
  m=0
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster%in% names(count_minus)])
  ##cluster
  for(clusid in names(count_minus[(count_minus / count_total) >0.6]) ){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,100000)
    minusduplic<-duplication_extract(endcluster1,miuscluster)  
    if(!isEmpty(minusduplic$dupli)){
      duplication<-rbind(duplication,minusduplic$dupli[,colnames(duplication)])
    }
    
   
    reverse_end<-reverse.region(miuscluster,chrid,3,"mius")
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
      
      local<-intersect.unit(middle,unique(middle$query_chr)) 
      changed_rows <- anti_join(local,middle)
      store<-rbind(store,inner_join(local,middle))
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
        store<-rbind(store,middle)
        store<-distinct(store)
      }}
    
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster<-reverse_end$endcluster1
      miuscluster$cluster<-paste(miuscluster$cluster, "8",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    
    m=m+1
  }
  
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster %in% names(count_minus)])

  
  SDRminudlist<-names(count_minus[(count_minus / count_total) >0.6]) ##cluster
  for(k in SDRminudlist){
    region<-endcluster1[endcluster1$cluster==k,]
    if(dim(region)[1]==1){
      next
    }
    region$cluster<-1:dim(region)[1]

    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    endcluster2before<-endcluster2$pos_end
    orientid="+"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) 
    reverse_end<-reverse.region(endcluster2,chrid,3,"init")
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){

      local<-intersect.unit(middle,unique(middle$query_chr)) 
      changed_rows <- anti_join(local,middle)
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
        store<-rbind(store,middle)
        store<-distinct(store)
      }
      
      if(length(which(middle$query_start>middle$query_end))!=0){
        middle<-middle[-(which(middle$query_start>middle$query_end)),]
        
      }
      if(nrow(middle)!=0){
        middle$orient<-"-"
        storesmall<-rbind(storesmall,middle[,colnames(storesmall)]) 
      }
      
    }
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
  }

  
  for(k in setdiff(unique(endcluster1$cluster),SDRminudlist)){
    region<-endcluster1[endcluster1$cluster==k,]
    region$cluster<-1:dim(region)[1]
    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    endcluster2before<-endcluster2$pos_end
    orientid="-"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) 
    reverse_end<-reverse.region(endcluster2,chrid,2,"init")
    middle<-reverse_end$reverse
    local<-intersect.unit(middle,unique(middle$query_chr)) 
    changed_rows <- anti_join(local,middle)
    if(nrow(changed_rows)!=0){
      changed_rows$anno<-'COMPLEX'
      store<-rbind(store,changed_rows[,colnames(store)])
      
    }

    #if(length(which(middle$query_start>middle$query_end))!=0){
    #  middle<-middle[-(which(middle$query_start>middle$query_end)),]
    #}
    
    if(length(which(middle$query_start>middle$query_end))!=0){
      subm=middle[(which(middle$query_start>middle$query_end)),]
      subm$change=subm$query_start-subm$query_end
      subm$len=subm$ref_end-subm$ref_start
      subm[subm$change/subm$len<0.05,]
      middle<-middle[-(which(middle$query_start>middle$query_end)),]
      middle<-rbind(subm[,colnames(middle)],middle)
      
    }
    
    if(nrow(middle)!=0){
      middle$orient<-"+"
      storesmall<-rbind(storesmall,middle[,colnames(storesmall)])
    }
    
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
    
  }
  cat("second cluster:SDR num",dim(store)[1],"\n")
  cat("second cluster:INV num",dim(inversion)[1],"\n")
  cat("second cluster:dup num",dim(duplication)[1],"\n")
  cat("second cluster:trans num",dim(smalltrans)[1],"\n")

	final_columns_order <- c(
		"ref_chr", "ref_start", "ref_end", 
		"query_chr", "query_start", "query_end", 
		"anno", "orient", "reflen", "querylen"
	)

	align_and_sanitize <- function(df, df_name, force_anno = NULL, force_orient = NULL, remove_first_row = FALSE) {
		if (!exists(df_name, envir = .GlobalEnv) || is.null(df) || nrow(df) == 0) {
			empty_df <- as.data.frame(matrix(ncol = length(final_columns_order), nrow = 0))
			colnames(empty_df) <- final_columns_order
			return(empty_df)
		}
		if (remove_first_row) {
			if (nrow(df) > 1) {
				df <- df[-1, , drop = FALSE]
			} else {
				empty_df <- as.data.frame(matrix(ncol = length(final_columns_order), nrow = 0))
				colnames(empty_df) <- final_columns_order
				return(empty_df)
			}
		}
		for (col in final_columns_order) {
			if (!col %in% colnames(df)) {
				df[[col]] <- NA
			}
		}
		if (!is.null(force_anno)) {
			df[["anno"]] <- force_anno
		}
		if (!is.null(force_orient)) {
			df[["orient"]] <- force_orient
		}
		df_cleaned <- df[, final_columns_order, drop = FALSE]
		return(df_cleaned)
	}

	if (exists("store") && !is.null(store) && nrow(store) > 0) {
		store <- store[store$ref_chr != "0" & !is.na(store$ref_chr), , drop = FALSE]
		store <- align_and_sanitize(store, "store", force_orient = "no")
	} else {
		store <- align_and_sanitize(NULL, "store")
	}

	inversion <- align_and_sanitize(inversion, "inversion", force_anno = "INV", force_orient = "-", remove_first_row = TRUE)

	if (exists("duplication") && !is.null(duplication) && nrow(duplication) > 0) {
		duplication_dup <- distinct(duplication)
		duplication <- align_and_sanitize(duplication_dup, "duplication_dup", force_anno = "DUP", remove_first_row = TRUE)
	} else {
		duplication <- align_and_sanitize(NULL, "duplication")
	}

	smalltrans <- align_and_sanitize(smalltrans, "smalltrans", force_anno = "TRANS", force_orient = "no", remove_first_row = TRUE)

	storesmall <- align_and_sanitize(storesmall, "storesmall", force_anno = "SDR_NM", remove_first_row = TRUE)

	if (exists("storehighdup_sdr") && !is.null(storehighdup_sdr) && nrow(storehighdup_sdr) > 0) {
		storehighdup_sdr <- storehighdup_sdr[storehighdup_sdr$ref_chr != "0" & !is.na(storehighdup_sdr$ref_chr), , drop = FALSE]
		storehighdup_sdr <- align_and_sanitize(storehighdup_sdr, "storehighdup_sdr", force_orient = "no")
	} else {
		storehighdup_sdr <- align_and_sanitize(NULL, "storehighdup_sdr")
	}

	INV_INV <- align_and_sanitize(INV_INV, "INV_INV")

	all <- rbind(store, inversion, duplication, smalltrans, storesmall, storehighdup_sdr, INV_INV)
	all <- distinct(all)

	if (nrow(all) == 0) {
		write.table(all, paste(args[4], chrid, "end.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
		next
	}

	all <- all[order(as.numeric(all$ref_start)), ]

	valid_query_chrs <- unique(all$query_chr[all$query_chr != "0" & !is.na(all$query_chr) & all$query_chr != ""])
	if (length(valid_query_chrs) > 0) {
		distance <- args[7]
		for (chr_child in valid_query_chrs) {
			sub_all <- all[all$query_chr == chr_child, , drop = FALSE]
			if (nrow(sub_all) > 0) {
				assign(chr_child, endfilter(sub_all, chrid, chr_child, distance))
			}
		}
	}

	data <- docall(all)

	if (!is.null(data) && nrow(data) > 0) {
		wrong_idx <- which(as.numeric(data$ref_start) > as.numeric(data$ref_end))
		if (length(wrong_idx) > 0) data <- data[-wrong_idx, , drop = FALSE]

		data <- distinct(data)

		zero_idx <- which(data$reflen == 0 & data$querylen == 0)
		if (length(zero_idx) > 0) data <- data[-zero_idx, , drop = FALSE]
	}

	COMPLEX <- data[data$anno == "COMPLEX" & !is.na(data$anno), , drop = FALSE]
	if (nrow(COMPLEX) > 0) {
		data <- data[data$anno != "COMPLEX" | is.na(data$anno), , drop = FALSE]
		newcomplex <- complexinte(COMPLEX)
		if (!is.null(newcomplex) && nrow(newcomplex) > 0) {
			newcomplex$anno <- "SDR_COMPLEX"
			newcomplex$orient <- "no"
			newcomplex$reflen <- as.numeric(newcomplex$ref_end) - as.numeric(newcomplex$ref_start)
			newcomplex$querylen <- as.numeric(newcomplex$query_end) - as.numeric(newcomplex$query_start)
			newcomplex_cleaned <- align_and_sanitize(newcomplex, "newcomplex")
			data <- rbind(data[, final_columns_order, drop = FALSE], newcomplex_cleaned)
		}
	}

	if (nrow(data) > 0) {
		data <- data[as.numeric(data$ref_start) <= as.numeric(data$ref_end) & 
			as.numeric(data$query_start) <= as.numeric(data$query_end), , drop = FALSE]
	}

	if (exists("highdupall") && !is.null(highdupall) && nrow(highdupall) > 0) {
		hidup <- highdupall[highdupall$ref_chr != "0" & !is.na(highdupall$ref_chr), , drop = FALSE]
		if (nrow(hidup) > 0) {
			hidup$query_chr <- ""
			hidup$query_start <- ""
			hidup$query_end <- ""
			hidup$anno <- "high-dup"
			hidup$orient <- ""
			hidup$reflen <- ""
			hidup$querylen <- ""
			hidup_cleaned <- align_and_sanitize(hidup, "hidup")
			data <- rbind(data[, final_columns_order, drop = FALSE], hidup_cleaned)
		}
	}

	if (nrow(data) > 0) {
		data <- data[!(data$ref_start == 0 & data$ref_end == 0), , drop = FALSE]
		data <- data[!(data$query_start == 0 & data$query_end == 0), , drop = FALSE]
	}

	if (nrow(data) > 0) {
		data <- data[, final_columns_order, drop = FALSE]
		numeric_cols <- c("ref_start", "ref_end", "query_start", "query_end", "reflen", "querylen")
		data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
			x_num <- as.numeric(x)
			if (any(is.na(x_num) & !is.na(x))) warning("Numeric conversion generated NAs in final printout.")
			format(x_num, scientific = FALSE, trim = TRUE, digits = 22)
		})
	} else {
		data <- as.data.frame(matrix(ncol = length(final_columns_order), nrow = 0))
		colnames(data) <- final_columns_order
	}

	write.table(data, paste(args[4], chrid, "end.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

}
})

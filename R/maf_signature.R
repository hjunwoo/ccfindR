#' Generate mutation count matrix for signature identifation
#' 
#' (Adapted from \code{TrinucleotideMatrix.R} in \code{maftools}))
#' @export
#' @param mafile MAF file to read
#' @param ref.genome Reference genome FASTA file.
#' @param prefix Chromosome name prefix to be added/substituted 
#'        to records in \code{maf}.
#' @param add Add the prefix. If \code{FALSE}, substitute instead.
#' @param renameChr A vector of two strings. All chromosome names
#'        matching the first will be substituted into second; e.g.,
#'        \code{renameChr=c('chr23','chrX')}.
mutcount <- function(mafile, ref.genome, prefix=NULL, add=TRUE,
                        kmer.size=3, skip.chr=NULL, 
                        remove.duplicates=TRUE,
                        useSyn=TRUE, rename.chr=NULL, gene.list=FALSE,
                        progress.bar=TRUE){
  
  maf <- read.table(mafile,header=TRUE,sep='\t',stringsAsFactors=FALSE)
  maf <- maf[maf$Variant_Type=='SNP',]
  if(remove.duplicates)
    maf <- maf[!duplicated(maf[,c('NCBI_Build','Chromosome','Start_Position',
                                  'End_position','Tumor_Sample_Barcode')]),]
  
  if(nrow(maf)==0) stop('No SNV present in maf file')
  
  up <- down <- 1
  
  if(!is.null(skip.chr))         # skip specified chromosomes
    maf <- maf[!maf$Chromosome %in% skip.chr,]
  if(!is.null(prefix))           # add prefix to chromosome names
    maf$Chromosome <- paste(prefix, maf$Chromosome, sep='')
    
  if(!is.null(rename.chr)){
    n <- length(rename.chr)
    for(i in seq_len(n))
      maf$Chromosome[maf$Chromosome==rename.chr[[i]][1]] <- rename.chr[[i]][2]
  }
 
  chrs <- unique(maf$Chromosome)
  
  if(class(ref.genome)!='DNAStringSet'){
    message('Reading ref. genome ...')            # read ref.genome
    if(class(ref.genome)=='character'){           # ref.genome path supplied
      ref <- Rsamtools::FaFile(file=ref.genome)
      ref <- Biostrings::getSeq(x = ref)
    }
    else if(class(ref.genome)=='BSgenome')
      ref <- Biostrings::getSeq(x = ref.genome)
    else
      stop('Inappropriate class object for ref.genome')
  }
  
  seq.lvl <- VariantAnnotation::seqlevels(ref)  # extract contigs from ref
  seq.lvl <- sapply(X=strsplit(x=seq.lvl, split=' ',fixed=TRUE), '[[',1)
  chrs.missing <- chrs[!chrs %in% seq.lvl]
  
  if(length(chrs.missing) > 0){
    warning(paste0('Chromosome ',chrs.missing,' not in ref and ignored'))
    maf <- maf[!maf$Chromosome %in% chrs.missing,]
    if(nrow(maf)==0) stop('No SNV in ref.genome')
  }
  
  dk <- as.integer((kmer.size-1)/2)
  extract.tbl <- data.frame(Chromosome=maf$Chromosome,
                    Start=maf$Start_Position-dk,
                    End=maf$Start_Position+dk,
                    Genes.affected=maf$Hugo_Symbol,
                    Reference_Allele=maf$Reference_Allele,
                    Tumor_Seq_Allele2=maf$Tumor_Seq_Allele2,
                    Tumor_Sample_Barcode=maf$Tumor_Sample_Barcode)
  
  ss <- Biostrings::subseq(x=ref[extract.tbl[,'Chromosome']],
                           start=extract.tbl[,'Start'],
                           end=extract.tbl[,'End'])
  
  polynuc <- as.character(ss)
  subst <- paste(extract.tbl$Reference_Allele, extract.tbl$Tumor_Seq_Allele2,
                        sep='>')
  subst.motif <- paste0(substr(x=polynuc,1,dk),'[',subst,']',
                        substr(polynuc,dk+2,2*dk+1))

  conv <- c('T>C','T>C','C>T','C>T','T>A','T>A','T>G','T>G','C>A',
            'C>A','C>G','C>G')
  names(conv) <- c('A>G','T>C','C>T','G>A','A>T','T>A','A>C','T>G',
                   'C>A','G>T','C>G','G>C')
  
  subst.type <- conv[subst]
  
  if(dk==1){
    complement <- c('A','C','G','T')
    names(complement) <- c('T','G','C','A')
  } else if(dk==2){
    complement <-        c('AA','AC','AG','AT',
                           'CA','CC','CG','CT',
                           'GA','GC','GG','GT',
                           'TA','TC','TG','TT')
    names(complement) <- c('TT','GT','CT','AT',
                           'TG','GG','CG','AG',
                           'TC','GC','CC','AC',
                           'TA','GA','CA','AA')
  }
  complemented.multiplets <- paste0(
    complement[substr(polynuc,dk+2,2*dk+1)],'[',subst.type,']',
    complement[substr(polynuc, 1, dk)])
  swap.ind <- which(substr(x=subst,1,1)%in% c('G','A'))
  subst.motif[swap.ind] <- complemented.multiplets[swap.ind]
  
  extract.tbl <- cbind(extract.tbl, data.frame(polynuc=polynuc,subst=subst,
                                               subst.motif=subst.motif))
  
  #-------
  subtype.levels <- mut_list(kmer.size)
  message('Creating mutation matrix..')
  
  samples <- unique(extract.tbl$Tumor_Sample_Barcode)
  mut.mat <- matrix(0, nrow=length(subtype.levels), ncol=length(samples))
  rownames(mut.mat) <- subtype.levels
  colnames(mut.mat) <- samples
  for(m in subtype.levels){
    x <- table(extract.tbl[extract.tbl$subst.motif==m,]$Tumor_Sample_Barcode)
    mut.mat[m,names(x)] <- x
  }
  mut.mat <- as(mut.mat,'dgCMatrix')
  if(!gene.list) return(mut.mat)
  
  gene.mat <- as.data.frame(
    data.table::dcast(data.table(extract.tbl),
                      formula=Tumor_Sample_Barcode~subst.motif,fill='',
                      value.var='Genes.affected', 
                      fun.aggregate=function(x){paste(x[x!='UnknownGene'],collapse=' ')}))
  rownames(gene.mat) <- gene.mat[,1]
  gene.mat <- t(gene.mat[,-1])
  gene.mat <- gene.mat[,match(colnames(mut.mat),colnames(gene.mat))]
  gene.mat <- gene.mat[match(rownames(mut.mat),rownames(gene.mat)),]

  return(list(count=mut.mat, genes=gene.mat))
}
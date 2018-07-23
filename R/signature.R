#' Generate mutation count matrix for signature identifation
#' 
#' (Adapted from \code{TrinucleotideMatrix.R} in \code{maftools})
#' @export
#' @param maf MAF object from \code{maftools}.
#' @param ref.genome Reference genome FASTA file.
#' @param prefix Chromosome name prefix to be added/substituted 
#'        to records in \code{maf}.
#' @param add Add the prefix. If \code{FALSE}, substitute instead.
#' @param renameChr A vector of two strings. All chromosome names
#'        matching the first will be substituted into second; e.g.,
#'        \code{renameChr=c('chr23','chrX')}.
kmer_matrix <- function(maf, ref.genome, prefix=NULL, add=TRUE,
                        kmer.size=3, ignoreChr=NULL, 
                        useSyn=TRUE, renameChr=NULL){
  
  maf.snp <- maftools::subsetMaf(maf=maf, includeSyn=useSyn, 
                       fields='Variant_Classification',
                       query="Variant_Type == 'SNP'")
  if(nrow(maf.snp)==0)
    stop('No snv left after filtering')
  
  up <- down <- 1
  
  if(!is.null(ignoreChr))
    maf.snp <- maf.snp[!maf.snp$Chromosome %in% ignoreChr]
  
  if(!is.null(prefix)){
    if(add)
      maf.snp$Chromosome <- paste(prefix, maf.snp$Chromosome,sep='')
    else
      maf.snp$Chromosome <- gsub(pattern=prefix, replacement='',
                                 x=maf.snp$Chromosome, fixed=TRUE)
  }
  if(!is.null(renameChr))
    maf.snp$Chromosome[maf.snp$Chromosome==renameChr[1]] <- renameChr[2]
 
  chrs <- unique(maf.snp$Chromosome)
  
  if(class(ref.genome)=='character'){       # ref.genome path supplied
    ref <- Rsamtools::FaFile(file=ref.genome)
    message(paste0('Reading ', ref.genome,'...')) # read ref.genome
    ref <- VariantAnnotation::getSeq(x = ref)
  } else
    ref <- ref.genome
  
  seq.lvl <- VariantAnnotation::seqlevels(ref)  # extract contigs from ref
  seq.lvl <- sapply(X=strsplit(x=seq.lvl, split=' ',fixed=TRUE),
                    '[[',1)
  names(ref) <- seq.lvl
  chrs.missing <- chrs[!chrs %in% seq.lvl]
  
  if(length(chrs.missing) > 0){
#    warning(paste0('Chr. in fasta: ', paste(seq.lvl, collapse=', ')))
#    warning(paste0('Chr. in input maf: ', paste(chrs, collapse=', ')))
#    warning(paste0('Chr. names in MAF must match those in ref. Ignoring ',
#                   nrow(maf.snp[Chromosome %in% chrs.missing]),
#                   'snvs from missing chr.',
#                   paste(chrs.missing, collapse=', ')))
    warning(paste0('Chromosome ',chrs.missing,' not in ref and ignored'))
    maf.snp <- maf.snp[!Chromosome %in% chrs.missing]
    if(nrow(maf.snp)==0) stop('No mutations left.')
  }
  
  dk <- as.integer((kmer.size-1)/2)
  extract.tbl <- data.table::data.table(Chromosome=maf.snp$Chromosome,
                    Start=maf.snp$Start_Position-dk,
                    End=maf.snp$Start_Position+dk,
                    Reference_Allele=maf.snp$Reference_Allele,
                    Tumor_Seq_Allele2=maf.snp$Tumor_Seq_Allele2,
                    Tumor_Sample_Barcode=maf.snp$Tumor_Sample_Barcode,
                    upstream=maf.snp$Start_Position-20,
                    downstream=maf.snp$End_Position+20)
  message("Extracting 5' and 3' adjacent bases...")
  ss <- Biostrings::subseq(x=ref[extract.tbl[,Chromosome]],
                           start=extract.tbl[,Start],
                           end=extract.tbl[,End])
  message('Extracting +/- 20bp around mutation...')
  updwn <- Biostrings::subseq(x=ref[extract.tbl[,Chromosome]],
                               start=extract.tbl[,upstream],
                               end=extract.tbl[,downstream])
  updwn.alphFreq <- 
    data.table::as.data.table(
      Biostrings::alphabetFrequency(x=updwn))[,.(A,C,G,T)]
  updwn.tnmFreq <- data.table::as.data.table(Biostrings::trinucleotideFrequency(
    x=updwn, step=1))
  
  extract.tbl[,polynucleotide:= as.character(ss)][,updwn:=as.character(updwn)]
  extract.tbl <- cbind(extract.tbl, updwn.alphFreq[,.(A,T,G,C)])
  extract.tbl <- cbind(extract.tbl, updwn.tnmFreq[,.(TCA,TCT,AGA,TGA)])
  extract.tbl[, tcw := rowSums(extract.tbl[, .(TCA,TCT)])]
  extract.tbl[, wga := rowSums(extract.tbl[, .(TGA,AGA)])]
  
  extract.tbl[,Substitution:=paste(extract.tbl$Reference_Allele,
            extract.tbl$Tumor_Seq_Allele2, sep='>')]
  extract.tbl$SubstitutionMotif <- paste(substr(
    x=as.character(extract.tbl$polynucleotide),1,dk),
    '[',extract.tbl$Substitution,']', 
    substr(as.character(extract.tbl$polynucleotide),dk+2, 2*dk+1),sep='')
  
  conv <- c('T>C','T>C','C>T','C>T','T>A','T>A','T>G','T>G','C>A','C>A','C>G','C>G')
  names(conv) <- c('A>G','T>C','C>T','G>A','A>T','T>A','A>C','T>G','C>A','G>T','C>G','G>C')
  
  extract.tbl$SubstitutionType <- conv[extract.tbl$Substitution]
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
  complemented.multiplets <- paste(
    complement[substr(x=as.character(extract.tbl$polynucleotide),dk+2,2*dk+1)],
    '[',extract.tbl$SubstitutionType,']',
    complement[substr(as.character(extract.tbl$polynucleotide), 1, dk)],
    sep='')
  swap.ind <- which(substr(x=extract.tbl$Substitution,1,1)%in% c('G','A'))
  swapSubTypeMotif <- extract.tbl$SubstitutionTypeMotif <-
    paste(substr(x=as.character(extract.tbl$polynucleotide),1,dk),'[',
          extract.tbl$SubstitutionType,']',
          substr(as.character(extract.tbl$polynucleotide), dk+2, 2*dk+1),sep='')
  swapSubTypeMotif[swap.ind] <- complemented.multiplets[swap.ind]
  
  extract.tbl$SubstitutionTypeMotif <- swapSubTypeMotif
  sub.levels <- extract.tbl[,.N, Substitution][,Substitution]

#  nt <- c('A','C','G','T')
#  stype <- c('[C>A]','[C>G]','[C>T]','[T>A]','[T>C]','[T>G]')
#  z <- vector('list',2*dk+1)
#  for(k in seq_len(2*dk)) z[[k]] <- nt
#  z[[2*dk+1]] <- stype
#  x <- expand.grid(z)
#  x <- x[,c(seq(dk+1,2*dk),2*dk+1,seq(1,dk))]
#  subtype.levels <- apply(x,1,paste, collapse='')
  
  subtype.levels <- mut_list(kmer.size)
  message('Creating mutation matrix..')
  extract.tbl.summary <- extract.tbl[,.N, by=list(Tumor_Sample_Barcode,
                                                   SubstitutionTypeMotif)]
  
  conv.mat <- as.data.frame(
    data.table::dcast(extract.tbl.summary, 
                      formula=Tumor_Sample_Barcode~SubstitutionTypeMotif,fill=0,
                      value.var='N', drop=FALSE))
  rownames(conv.mat) <- conv.mat[,1]
  conv.mat <- conv.mat[,-1]
  
  colOrder.missing <- subtype.levels[!subtype.levels %in% colnames(conv.mat)]
  
  if(length(colOrder.missing) > 0){
    zeroMat <- as.data.frame(matrix(data = 0, nrow=nrow(conv.mat),
                                    ncol = length(colOrder.missing)))
    colnames(zeroMat) <- colOrder.missing
    conv.mat <- cbind(conv.mat, zeroMat)
  }
  conv.mat <- as.matrix(conv.mat[,match(subtype.levels, colnames(conv.mat))])
  conv.mat[is.na(conv.mat)] <- 0
  conv.mat <- t(as(conv.mat, 'dgCMatrix'))
  
  message(paste('matrix of dimension ', nrow(conv.mat), 'x', ncol(conv.mat), sep=''))
  
  return(conv.mat)
}
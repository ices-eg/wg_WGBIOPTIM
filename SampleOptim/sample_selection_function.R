porSexo<-function(parte,sexR,numTotF,numTotM){
  q_s1<-NULL;
  if(dim(parte)[1]>0){
    if(sexR){ # escolher de acordo com o sexRatio
      #escolher femeas
      numF<-numTotF;
      parteF<-parte[parte$SEXO=="F",];
      if(dim(parteF)[1]>numF) q_s1_F<-sample(parteF$ID_BIO_PEIXES, numF, replace=F)
      else q_s1_F<-parteF$ID_BIO_PEIXES;
      numF<-length(q_s1_F)

      #escolher machos
      numM<-numTotM;
      parteM<-parte[parte$SEXO=="M",];
      if(dim(parteM)[1]>numM) q_s1_M<-sample(parteM$ID_BIO_PEIXES, numM, replace=F)
      else q_s1_M<-parteM$ID_BIO_PEIXES;
      numM<-length(q_s1_M)
      q_s1<-c(q_s1_F, q_s1_M);

      #escolher indeterminados para completar o numero especificado em quantidade
      numI<-numTotF+numTotM-numF-numM;
      if(numI>0){
        parteI<-parte[parte$SEXO=="I",];
        if(dim(parteI)[1]>numI) q_s1_I<-sample(parteI$ID_BIO_PEIXES, numI, replace=F)
        else q_s1_I<-parteI$ID_BIO_PEIXES;
        q_s1<-c(q_s1, q_s1_I)
      }
    }
    else{ # escolher sem especificar sexRatio
      quantidade<-numTotF+numTotM
      if(dim(parte)[1]>quantidade){
        q_s1<-sample(parte$ID_BIO_PEIXES, quantidade, replace=F);
        if(quantidade != length(q_s1)) print(paste("ERRO na quantidade:",quantidade," diferente de:",parte$ID_BIO_PEIXES));
      }else{
        q_s1<- parte$ID_BIO_PEIXES;
      }
    }
  }
  return(q_s1);
}


##
#porClassesPorSexo
#
# Escolhe os indiv?duos por classe de comprimento e por sexo (sexR=T)
#
porClassesPorSexo<-function(todos,classes,sexR,numTotF,numTotM, verbose=F){
  idRes<-NULL;
  for(cl in classes){   #por cada classe escolhida
    parte<-todos[todos$C_CLASSE == cl, ];
    q_s1<-porSexo(parte,sexR,numTotF,numTotM)
    idRes<-c(idRes, q_s1);
  }
  if(verbose) print(paste("N?mero de individuos: ",length(idRes)))
  return(idRes);
}

##
#porPortosPorClassesPorSexo
#
# Escolhe os indiv?duos por porto, por classe de comprimento e por sexo (sexR=T) (numTotF - n?mero total de f?meas, numTotM - n?mero total de machos)
# S?o escolhidos (numTotF+numTotM) de indiv?duos por porto de cada classe de comprimento (distUniPorto=FALSE) (valor m?ximo = numeroPortos*numeroClasses*quantidade)
# Caso distUniPorto=TRUE, s?o escolhidos quantidade de indiv?duos por cada classe de comprimento escolhidos  (valor m?ximo = numeroClasses*quantidade)
#  com uma distribui??o uniforme de entre os portos.
#
# todos - cont?m os indiv?duos escolhidos de acordo com o padr?o temporal seleccionado
# distUniPorto - determina como o porto influencia as amostras seleccionadas
# conjClasses - conjunto de classes
#
porPortosPorClassesPorSexo<-function(todos,classes,sexR,sexRatio,quantidade,distUniPorto, verbose=F){
  idRes<-NULL;
  conjPortos<-unique(todos$PORTO_NOME);
  numPortos<-length(conjPortos)

  if(distUniPorto){
    quantidadePorPorto <- trunc(dunif(seq(1,1, length = numPortos))/numPortos*quantidade)
    falta<-quantidade-sum(quantidadePorPorto)
    resto<-c(rep(1,falta),rep(0,(numPortos-falta)))
    quantidadePorPorto<-quantidadePorPorto+sample(resto,numPortos)
  } else quantidadePorPorto<-array(quantidade, numPortos)
  numF<-ceiling(quantidadePorPorto*sexRatio)
  numM<-quantidadePorPorto-numF
  nporto<-1
  for(p in conjPortos){
    if(verbose) print(paste("Porto: ", p));
    pPorto<-todos[todos$PORTO_NOME == p, ];
    q_s1<-porClassesPorSexo(pPorto,classes,sexR,numF[nporto],numM[nporto], verbose=verbose)
    idRes<-c(idRes, q_s1);
    nporto<-nporto+1
  }
  return(idRes);
}

##
#semestre
#
# devolve o semetre (1 ou 2) ao qual pertence o m?s
semestre<-function(mes){
  return(ceiling(mes/6));
}

##
#trimestre
#
# devolve o trimestre (1,2,3 ou 4) ao qual pertence o m?s
trimestre<-function(mes){
  return(ceiling(mes/3))
}

############################ ? L T I M A Vers?o ###############################

# Devolve um vector com os $ID_BIO_PEIXES seleccionados aleatoriamente para a quantidade referida em cada classe de comprimento
# O per?odo temporal ? o escolhido
# Podem entrar diversos anos, o resultado do per?odo escolhido ? em fun??o de cada ano
#
# sexRatio - propor??o entre f?meas e machos (ex: 1 - 100% f?meas, 0 - 100% machos, 0.2 - 20% f?meas 80% machos)
#   no caso de sexRatio <0 indica que a escolha n?o depende do sexo dos indiv?duos
# conjClasses - conjunto de classes de comprimento
# quantidade - n?mero de indiv?duos a escolher por classe de comprimento
# tm - op??o temporal
#     1S - escolher por ano e s? dentro do 1? semestre
#     2S - escolher por ano e s? dentro do 2? semestre
#     S  - escolher por ano e dentro de cada semestre
#     1T - escolher por ano e s? dentro do 1? trimestre
#     2T - escolher por ano e s? dentro do 2? trimestre
#     3T - escolher por ano e s? dentro do 3? trimestre
#     4T - escolher por ano e s? dentro do 4? trimestre
#     T  - escolher por ano e dentro de cada trimestre
#     A  - escolher por ano
#
# S?o acrescentadas duas colunas aos dados: SEMESTRE e TRIMESTRE para facilitar a escolha
#

amostraTemporal<-function(tab, quantidade, conjClasses, sexRatio=NaN, tm="T", porto=TRUE, distUniPorto=TRUE, verbose=F){
  tab$SEMESTRE<-semestre(tab$MONTH);
  tab$TRIMESTRE<-trimestre(tab$MONTH);
  original<-tab;
  idRes<-NULL
  tab<-switch(tm,
              "1S" = tab[tab$SEMESTRE==1,],   # 1? semestre
              "2S" = tab[tab$SEMESTRE==2,],   # 2? semestre
              "1T" = tab[tab$TRIMESTRE==1,],  # 1? trimestre
              "2T" = tab[tab$TRIMESTRE==2,],  # 2? trimestre
              "3T" = tab[tab$TRIMESTRE==3,],  # 3? trimestre
              "4T" = tab[tab$TRIMESTRE==4,],  # 4? trimestre
              "A"  = tab,                     # Ano
              "S"  = tab,                     # Semestre
              "T"  = tab                      # Trimestre
  );
  if(is.null(tab)){
    print(paste("ERRO: Per?odo temporal n?o foi inserido correctamente: <", tm,">"));
    return (NULL);
  }
  if(verbose) print(paste("Per?odo de tempo: ",tm));
  if(verbose) print(paste("N?mero de indiv?duos: ",dim(tab)[1]));
  sexR<-TRUE;
  numTotF<-quantidade;
  numTotM<-0;
  if(!is.nan(sexRatio) & (sexRatio>=0) & (sexRatio<=1)){ # escolher de acordo com o sexRatio
    #calcular n?mero de f?meas e macho a escolher
    numTotF<-ceiling(quantidade * sexRatio);
    numTotM<-quantidade -numTotF;
    if(verbose) print(paste("SexRatio: ", sexRatio, "Fem?as:",numTotF,"Machos:",numTotM));
  } else { sexR<-FALSE; if(verbose) print("SexRatio: Ignore"); }
  if(verbose) print(paste("N?mero de classes: ", length(conjClasses)));
  if(tm=="S"){
    ano<-unique(tab$ANO);
    for(a in ano)
      for(s in 1:2){
        if(verbose) print(paste("ano:",a," ",s,"S",sep=''))
        parte<-tab[tab$ANO==a & tab$SEMESTRE==s,]
        if(porto)
          q_s1<-porPortosPorClassesPorSexo(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
        else  q_s1<-porClassesPorSexo(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
        idRes<-c(idRes, q_s1);
      }
  }
  else if(tm=="T"){
    ano<-unique(tab$ANO);
    for(a in ano)
      for(t in 1:4){
        if(verbose) print(paste("ano:",a," ",t,"T",sep=''))
        parte<-tab[tab$ANO==a & tab$TRIMESTRE==t,]
        if(porto)
          q_s1<-porPortosPorClassesPorSexo(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
        else  q_s1<-porClassesPorSexo(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
        idRes<-c(idRes, q_s1);
      }
  }
  else if(tm=="A"){
    ano<-unique(tab$ANO);
    for(a in ano){
      if(verbose) print(paste("ano:",a,sep=''))
      parte<-tab[tab$ANO==a,]
      if(verbose) print(dim(parte))
      if(porto)
        q_s1<-porPortosPorClassesPorSexo(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
      else  q_s1<-porClassesPorSexo(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
      if(verbose) print(length(q_s1))
      idRes<-c(idRes, q_s1);
    }
  }
  else{ #tab cont?m a escolha de acordo com tm= "1S","2S","1T","2T","3T" ou "4T"
    ano<-unique(tab$ANO);
    for(a in ano){
      if(verbose) print(paste("ano:",a," ",tm,sep=''))
      parte<-tab[tab$ANO==a,]
      if(verbose) print(dim(parte))
      if(porto)
        q_s1<-porPortosPorClassesPorSexo(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
      else  q_s1<-porClassesPorSexo(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
      idRes<-c(idRes, q_s1);
    }
  }
  return(original[original$ID_BIO_PEIXES %in% idRes,])
}

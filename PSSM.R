# ATENÇÃO: A funcao de PSSM requer a instalacao do pacote do Bioconductor msa
# que realizara o multialinhamento, caso nao o possua na sua maquina, rodar as linhas abaixo no seu console:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("msa")
# Todos os outros metodos nesse script sao nativos ao R.

# Desenvolvedores: Joao Vitor Cavalcante (jvfecav@gmail.com)
#                  Vitor Gabriel Saldanha (vitor.saldanha11@gmail.com)
# Data: 09/09/2019


pssm <- function(sequence){
  #------------------------------------------ Preparacao do alinhamento ----------------------------------------
  
  #Carrega Pacote msa para fazer multialinhamento 
  library(msa)  
  
  #Alinha as proteinas com a ordem do input e algoritmo Muscle
  alignment <- msa(sequence, type="protein", order="input", "Muscle")
  
  #Transforma alinhamento contido no subobjeto unmasked em matriz 
  matriz_total <- as.matrix(alignment@unmasked)
  
  #Cria um vetor com todos os aminacidos para usar como nome da pssm
  amino <- strsplit("AGILPVFWYDERHKSTCMNQ-", "")[[1]]
  
  #Se no alinhamento nao existirem todos aminoacidos, diminui-se o escopo de aminoacidos
  names <- list(amino[amino %in% matriz_total])
  
  #Cria matriz vazia para ser adicionar a contagem a cada iterador do for abaixo
  #Essa matriz tem o mesmo tamanho que os names, mesmo numero de coluna que o alinhamento, 
  #tem as linhas com nomes dos aminoacidos e eh preenchida de 0
  #OBSERVACAO: Esse objeto se chama de pssm, mas ate la ele vai armazenar uma pfm e um ppm.
  pssm <- matrix(data=0, nrow = length(names[[1]]), ncol = ncol(matriz_total), dimnames = names)
  
  #Transforma a matriz anterior em data.frame (tabela)
  pssm <- as.data.frame(pssm)
  
  #------------------------------------------------- Calculo da PFM ----------------------------------------------------------
  #For que vai passar por cada coluna do alinhamento para fazer a contagem dos aminoacidos
  for(i in 1:ncol(matriz_total)){
    
    #A cada coluna ele usa a funcao table, a qual conta a frequencia de elementos dentro de um vetor
    #Esse objeto eh um vetor (sequencia de numeros) que cada um tem seu nome (aminoacido correspondente a contagem)
    counts <- table(matriz_total[,i])
    
    #Percorre o vetor counts a fim de criar a pfm
    for(int in 1:length(names(counts))){
      
      #Para cada elemento do vetor ele pega o nome do aminoacido
      amino_name <- names(counts[int])
      
      #Tambem para cada elemento ele pega a contagem do aminacido
      amino_count <- unname(counts[int])
      
      #Adiciona a contagem do aminoacido na linha "amino_name" (objeto que contem o nome do aminoacido)
      #E na coluna i (correspondente ao for mais exterior que passa por cada coluna)
      #No fim dos loops esse objeto chamado "pssm" que continha 0 agora eh uma pfm
      pssm[amino_name, i] <- amino_count
    }
  }
  
  #Tira a linha que contem a contagem de gaps
  pssm <- pssm[rownames(pssm)!="-",]
  
  #---------------------------------------------- Calculo da PPM ----------------------------------------------------
  
  #Divide-se a pfm (armazenada dentro do objeto chamado pssm) pelo numero de linhas do alinhamento para obter a PPM
  pssm <- pssm/nrow(matriz_total)
  
  #Atribui a todos os elementos que tem o valor de 0 a probabilidade de 0.001 (pseudocontagem)
  pssm[pssm==0] <- 0.001
  
  #--------------------------------------------- Calculo da PSSM ----------------------------------------------------
  
  #Divide todos os ementos por 0.05 e calcula o log de 2 para obter a PSSM final
  pssm <- log2(pssm/0.05)
  
  return(pssm)
}


#Sequencias de input
sequence <- c("CSEDWVGYQRKCYFISTVKRSWTSAQNACSE", "CSDDWIGHKGKYYLISKKTKNWTLAQNFCSK", "CPDDWIGYQTKCYFISKKTKNWTLAQSFCSK", "CKNEWFSYNGKCYFFSTTTKTWALAQKSCSE", "CKNEWISYKRTCYFFSTTTKSWALAQRSCSE")
pssm <- pssm(sequence)

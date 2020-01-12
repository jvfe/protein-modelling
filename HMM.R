# ATENÇÃO: A funcao de HMM para multialinhamento requer a instalacao do pacote do Bioconductor msa
# que realizara o multialinhamento, caso nao o possua na sua maquina, rodar as linhas abaixo no seu console:
# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("msa")
# Todos os outros metodos nesse script sao nativos ao R.

# Desenvolvedores: Joao Vitor Cavalcante (jvfecav@gmail.com)
#                  Vitor Gabriel Saldanha (vitor.saldanha11@gmail.com)
# Data: 09/09/2019

#--------------------------------------- Criar HMM -------------------------------------------------

createHMM <- function(proteins){
  
  # Carrega Pacote msa para fazer multialinhamento
  library(msa)  
  
  # Alinha as proteinas com a ordem do input e algoritmo Muscle
  alignment <- msa(sequence, type="protein", order="input","Muscle")
  
  # Transforma alinhamento contido no subobjeto unmasked em matriz 
  matriz_total <- as.matrix(alignment@unmasked)
  
  #------------------------------- Criar Matriz dos Estados ---------------------------------------
  
  # Cria matriz do mesmo tamanho do alinhamento (mesmo number of row e columns) preenchida com M, 
  # a qual vai ser usada posteriormente para conter o estado de cada aminoacido
  matriz_estados <- matrix("M", nrow=nrow(matriz_total), ncol=ncol(matriz_total))
  
  # Da o nome das sequencias do alinhamento a matriz de estados tambem
  row.names(matriz_estados) <- row.names(matriz_total)
  
  # Descobre os estados escondidos, o que eh insercao e delecao
  # Foi considerado que gaps mais que 20% caracterizam a coluna como estado escondido
  # Esse for abaixo procura por todas as colunas do alinhamento essa condicao
  for(i in 1:ncol(matriz_total)){
    
    # Transforma as matrizes em data.frames (tabelas), dizendo para não considerar os elementos fatores (categorias)
    matriz_total <- as.data.frame(matriz_total, stringsAsFactors= F)
    matriz_estados <- as.data.frame(matriz_estados, stringsAsFactors= F)
    
    # Nos locais que tem gap ("-" dentro da matriz do alinhamento) ele poe um D (delecao) dentro da matriz dos estados de cada aminoacido
    matriz_estados[matriz_total == "-"] <- "D"
    
    # Procura por colunas que tenham mais de 20% de gaps - se sim entra no if
    if(sum(matriz_total[,i] == "-")/nrow(matriz_total) > 0.2){
      
      # Muda o nome da coluna em que foi visto o estado escondido para H + num da coluna. Ex: coluna 2 daria um H2
      colnames(matriz_estados)[i] <- paste("H", i, sep = "")
      colnames(matriz_total)[i] <- paste("H", i, sep = "")
      
      # Ele percorre todas as linhas do estado escondido, anotando o que eh delecao e o que eh insercao
      for(iterador in 1:nrow(matriz_estados)){
        if(matriz_total[iterador,i] == "-"){
          matriz_estados[iterador,i] <- "D"
        } else {
          matriz_estados[iterador,i] <- "I"
        }
      }
    }
  }
  
  # Cria uma lista que contem duas tabelas: matriz_total (alinhamento) e matriz_estados(estado de cada aminoacido do alinhamento)
  matrizes <- list(alignment = matriz_total, states = matriz_estados)
  
  # Funcao retorna essa lista
  return(matrizes)
}


#----------------------------- Probabilidades de Emissao ----------------------------------------

emissionProb <- function(matrizes){

  # Cria lista vazia para ser adiciona um elemento a cada iterador do for abaixo
  lista <- list()
  
  # Transforma o alinhamento em data.frame(tabela)
  df <- as.data.frame(matrizes$alignment)
  
  # Para cada coluna do alinhamento calcula frequencia, contagem total e probabilidade
  for(i in 1:ncol(df)){
    
    # Calcula Frequencia da coluna - table eh uma funcao que faz a contagem de cada elemento de um vetor
    contagem_do_hmm <- as.data.frame(table(df[,i]))
    
    # Se existir um gap na coluna que lista os elementos (Var1), retira a contagem dos gaps do objeto
    if(contagem_do_hmm$Var1 == "-"){
      contagem_do_hmm <- contagem_do_hmm[!contagem_do_hmm$Var1 == "-",]
    }
    
    # Cria um nome unico para cada coluna
    nome <- paste("state", i, sep="_")
    
    # Adiciona a lista (inicialmente vazia) o elemento com o nome dentro de "name". Como nao existe o elemento,
    # ele eh criado e armazena a frequencia
    lista[[nome]] <- contagem_do_hmm
    
    # Calcula contagem total de aminoacidos naquela coluna
    n_total <- sum(lista[[nome]]$Freq)
    
    # Calcula a probabilidade de cada aminoacido - sua frequencia pela contagem total da coluna
    # e adiciona esse valor a uma coluna criada nesse momento chamda "Prob"
    lista[[nome]]$Prob <- lista[[nome]]$Freq / n_total
    
  }
  
  # Retorna essa lista que no fim vai ter todas as colunas em elementos e cada elemento vai possuir a frequencia 
  # e probabilidade de cada coluna
  return(lista)
  
}


#----------------------------- Probabilidades de Transicao ----------------------------------------

transProb <- function(matriz = matrizes){
  
  # Cria lista vazia para ser adiciona um elemento a cada iterador do for abaixo
  list_trans <- list()
  list_final <- list()
  
  # Cria dataframe(tabela) referente ao alinhamento
  df <- data.frame(matrizes$alignment, stringsAsFactors = F)
  
  # Loop para calcular frequencia e probabilidade de transicao
  # Obs: O loop vai de 1 até (ncol - 1) porque a probabilidade de transmissao nao existe na ultima coluna
  for(k in 1:(ncol(df[,-1]))) {
    
    # Cria data.frame para contar a frequencia de cada amioacido em uma coluna (funcao table)
    # A funcao começa a partir da segunda coluna porque nao existe probabilidade de transicao na 1a coluna
    contagem_transicao <- data.frame(table(df[,k+1]), stringsAsFactors = F)
    
    # Aqui foi feita uma tratativa para quando o estado anterior a um aminoacido for um gap (delecao), retirar
    # a contagem desse aminoacido que vem apos o gap.
    # Condicaçao imposta: if(coluna anterior tiver um gap & o estado dessa coluna anterior for visivel)
    # Exemplo: T - C = Se a segunda coluna for visivel, deve-se retirar a contagem de C em 1 porque a coluna 2
    #                  não existe. Entao nao eh uma tranmissao M-M, eh uma transmissao D-M.
    #                  Se a segunda coluna for invisivel, mantem-se a contagem porque C veio de T (transicao M-M).
    if(any(df[,k]=="-") & grepl("V",colnames(df)[k])){
      
      # Descobre qual aminoacido posterior ao gap (armazena uma string)
      # Subset das linhas: qual linha esta o gap; Subset da coluna: coluna posterior ao gap ( k + 1)
      idx <- df[df[,k] == "-", k + 1]
  
      # Busca na tabela de contagem qual a contagem do aminoacido armazenado em idx // Obs: 2 porque a segunda coluna eh a de count
      cont_posDel <- contagem_transicao[contagem_transicao$Var1 == idx, 2] 
      
      # Ele reatribui o calor de contagem para o valor encontrado anteriormente - 1
      contagem_transicao[contagem_transicao$Var1 == idx, 2] <- cont_posDel - 1
    }
    
    # Cria um nome unico para cada coluna do alinhamento
    # Assim: state_1 = coluna 2; state_2 = coluna 3 ....
    name <- paste("state", k, sep="_")
    
    # Cria dentro da lista_trans (inicialmente vazia) o elemento com string armazenado dentro de "name" e atribui a esse 
    # elemento a tabela de contagem
    list_trans[[name]] <- contagem_transicao
    
    # Calcula quanto eh o somatorio de todos aminoacidos (e possiveis gaps) existentes naquela coluna
    n_tot <- sum(list_trans[[name]]$Freq)
    
    # Calcula as probabilidades
    # Se há "-" na coluna, subtrai do n_tot o somatório do numero de vezes que  "-" ocorre (frequencia de gaps).
    # Se não há, realiza-se a divisão pelo n total.
    if(any(list_trans[[name]]$Var1=="-")){
      
      # Cria uma coluna nova dentro do objeto da lista chamada Prob que vai armazenar as probabilidades
      # A probabilidade eh calculada como a frequencia de cada aminoacido vididido pelo (n_tot - somatorio de gaps) 
      list_trans[[name]]$Prob <- list_trans[[name]]$Freq / (n_tot - list_trans[[name]]$Freq[list_trans[[name]]$Var1=="-"])
      
      # Caso o valor iterado se refira a um gap, conta-se como transição para um estado de deleção,
      # portanto, dividimos a frequencia de gaps pelo n_total
      list_trans[[name]]$Prob[list_trans[[name]]$Var1=="-"] <- list_trans[[name]]$Freq[list_trans[[name]]$Var1=="-"] / n_tot
    } else {
      
      # Divisão normal quando nao existe gap na coluna
      list_trans[[name]]$Prob <- list_trans[[name]]$Freq / n_tot
    }
    
    # Cria lista temporaria apenas para facilitar a sintaxe da linha posterior 
    temp_list <- list_trans[[name]]
    
    # Soma as probabilidades de todos os simbolos que nao "-" (gap) e armazena em visible.visible
    # que vai representar a transicao de estado M-M para o caso de dois estados visiveis ou ainda
    # V-I ou I-V, se existir uma coluna invisivel
    visible.visible <-  sum(temp_list[temp_list$Var1 != "-",][3])
    
    # Atribui um pseudocount a transicao de visivel para gap
    visible.gap <- -Inf
    
    # Se algum dos elementos da coluna for um gap, ele entra no if
    if(any(list_trans[[name]]$Var1=="-")){
      
      # Muda a probabilidade de -Infinito para a probabilidade de transicao do gap calculada
      # anteriormente
      visible.gap <- sum(temp_list[temp_list$Var1 == "-",][3])
      
      # Subtrai do visible.visible pela probabilidade do visible.gap (Assim a soma das probabilidade da um)
      visible.visible <- visible.visible - visible.gap
    }
    
    # Criacao da lista final, com os nomes das respectivas transicoes
    list_final[[name]] <- c(visible.visible = visible.visible, visible.gap = visible.gap)
    
  }  
  
  return(list_final)
  
}

#------------------------------------------- TESTE DO CÓDIGO -------------------------------------------

# Sequencia de input
sequence <- c("SYKHFTYL", "NYGPFTFL", "SYGPLFL", "TYEQLSFL", "KFAGNVDFL", "QYGSQVTFA", "QYKGDLSLV")

# Roda a funcao para criar o HMM
matrizes <- createHMM(sequence)
matrizes

# Roda a função para obter as probabilidades de emissao
emissao <- emissionProb(matrizes)
emissao

# Roda a funcao para obter as probabilidades de transicao
transicao <- transProb(matrizes)
transicao



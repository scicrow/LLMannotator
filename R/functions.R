
# LLM embeddings for the TOterms are already available in this github repo 
# under the 'embeddings' directory.
# Only run this funciton if you want to re-obtain the embeddings yourself. 
# 
embed_TOterms <- function()
{
  TOterms <- data.table::fread("../ontology/TO_termsdict.csv", encoding = "UTF-8")
  TOterms <- mutate(TOterms, Synonyms = str_remove_all(Synonyms, "(related)"))
  TOterms <- mutate(TOterms, bigstring = paste(`Preferred Label`, Definitions, sep=". "))
  TOterms <- dplyr::filter(TOterms, ID %in% TO$id)
  
  # embed each TO term's label and description using OpenAI
  embedded <- openai::create_embedding(model = "text-embedding-3-large",
                                       input = utf8::utf8_encode(TOterms %>% pull(bigstring)) )
  
  TOterms$embedding <- embedded$data$embedding
  return(TOterms)
}

embed_descriptors <- function()
{
  TAIR.embed <- openai::create_embedding(model = "text-embedding-3-large",
                                         input = utf8::utf8_encode(TAIR.unique %>%
                                                                   pull(PHENOTYPE)) )
  return(TAIR.embed)
}

embed_concepts <- function()
{
  TAIR.concept.embed <- openai::create_embedding(model = "text-embedding-3-large", 
                                      input = utf8::utf8_encode(concepts %>% 
                                                                  pull(embedtext)) )
  return(TAIR.concept.embed)
}


# this function reads through a data.frame of descriptors and sends each one to 
# GPT-4o to get it parsed into concepts. Each one is then written to file 
# defined by the outdir and fname parameters. 


openai_getconcepts <- function(descriptors, prompt = NULL, outdir="outputs", fname = "default", descriptorfields = c(2)) {
  for(i in 1:nrow(descriptors)) 
  {
    r <- descriptors[i,]
    print(paste0(r[[1]],": ",r[[descriptorfields]]))
    print("===========")
    
    messages = list(
      list("role"="system", "content"=prompt),
      list("role"="user", "content"= paste0(r[[1]],": ",r[[descriptorfields]]))
    )
    #print(messages)            
    out <- openai::create_chat_completion(model = "gpt-4o",
                                          temperature = 0,
                                          max_tokens = 4000,
                                          messages = messages
    )
    
    Sys.sleep(1) # this is to prevent too many OpenAI requests per minute
    print(out)
    if(!is.null(out) & !is.na(out$choices$message.content)){  # if we received concepts...
      out$choices$message.content <- paste0(out$choices$message.content, "\n")
      write_file(out$choices$message.content, paste0(DIRGPT, fname, ".", Sys.Date(), ".txt"), append = TRUE)
    } 
    else {  # ...otherwise, try again with more temperature
      print("*** TRYING AGAIN")
      out <- openai::create_chat_completion(model = "gpt-4o",
                                            temperature = 0.4,
                                            max_tokens = 1000,
                                            messages = messages
      )
      Sys.sleep(1)
      print(out)
      if(!is.null(out) & !is.na(out$choices$message.content)){
        out$choices$message.content <- paste0(out$choices$message.content, "\n")
        write_file(out$choices$message.content, paste0(outdir, fname, ".", Sys.Date(), ".txt"), append = TRUE)
      } 
    }
  }
  
}


openai_filterterms <- function(descriptors, prompt = NULL, outdir="outputs", fname = "default", descriptorfields = c(2)) {
  for(i in 1:nrow(descriptors)) 
  {
    r <- descriptors[i,]
    print(paste0(r[[1]],": ",r[[descriptorfields]]))
    print("===========")
    
    messages = list(
      list("role"="system", "content"=utf8::utf8_encode(prompt)),
      list("role"="user", "content"= utf8::utf8_encode(paste0(r[[1]],": ",r[[descriptorfields]])))
    )
    #print(messages)            
    out <- openai::create_chat_completion(model = "gpt-4o",
                                          temperature = 0,
                                          max_tokens = 4000,
                                          messages = messages
    )
    Sys.sleep(1)
    print(out)
    if(!is.null(out) & !is.na(out$choices$message.content)){  # if we received a valid response
      out$choices$message.content <- paste0(out$choices$message.content, "\n")
      write_file(out$choices$message.content, paste0(outdir, fname, ".", Sys.Date(), ".txt"), append = TRUE)
    } 
    else {  # ...otherwise, try again with more temperature
      print("*** TRYING AGAIN")
      out <- openai::create_chat_completion(model = "gpt-4o",
                                            temperature = 0.4,
                                            max_tokens = 4000,
                                            messages = messages
      )
      Sys.sleep(1)
      print(out)
      if(!is.null(out) & !is.na(out$choices$message.content)){
        out$choices$message.content <- paste0(out$choices$message.content, "\n")
        write_file(out$choices$message.content, paste0(outdir, fname, ".", Sys.Date(), ".txt"), append = TRUE)
      } 
    }
  }
  
}

getscore <- function(annotation, itemscol = NULL, preproc = NULL, annotator = NULL, dataset = NULL, ontology, ic)
{
  items <- annotation %>% pull(itemscol) %>% unique()
  
  score <- foreach(j = items, .combine = rbind) %do%
    {
      goldterms    <- filter(annotation, item_id==j) %>% 
        pull(Term.gold) %>% 
        unique()
      mappedterms  <- filter(annotation, item_id==j, !is.na(Term.auto)) %>% 
        pull(Term.auto) %>% 
        unique()
      
      if(length(mappedterms) == 0) {
        mappedterms <- NA
      }
      if(length(goldterms) == 0) {
        mappedterms <- NA
      }
      
      
      jaccard  <- length(intersect(mappedterms,goldterms)) / length(union(mappedterms,goldterms))   # jaccard 
      
      semsim = 0
      if(sum(is.na(goldterms))==1 & sum(is.na(mappedterms))==1) {
        semsim = 1
      } else if(sum(is.na(mappedterms))) {
        semsim = 0
      } else if(sum(is.na(goldterms))) {
        semsim = 0
      } else {
        semsim <- ontologySimilarity::get_sim_grid(ontology = ontology, 
                                                   information_content = ic,
                                                   term_sets = list("A" = mappedterms, 
                                                                    "B"= goldterms))[1,2]
      }
      
      recall     <- length(intersect(mappedterms,goldterms)) / length(goldterms)
      precision  <- length(intersect(mappedterms,goldterms)) / length(mappedterms)
      
      row1 <- data.frame(item_id = j, 
                         num_goldterms = length(goldterms), 
                         num_mappedterms = length(mappedterms),
                         goldterms = paste(goldterms, collapse = ","),
                         mappedterms = paste(mappedterms, collapse = ","),
                         recall = recall,
                         precision = precision,
                         jaccard = jaccard,
                         semsim = semsim,
                         preprocess = preproc,
                         annotator = annotator,
                         dataset = dataset)
      
      row1
      
    }
  
  return(score)
}

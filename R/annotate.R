library(tidyverse)
library(foreach)
library(iterators)
library(openai)
library(ontologyIndex)
library(ontologySimilarity)
library(data.table)
source("functions.R")

#########################################
#
# initialise 
#
########################################

# Need to run this line to access OpenAI API
# If you don't have a key, you need to register an account with OpenAI
Sys.setenv(OPENAI_API_KEY = "put your key here")


# load ontology graph
TO    <- ontologyIndex::get_ontology(file = "../ontology/to.obo")
# get information content for semsim
TO.ic <- ontologySimilarity::descendants_IC(TO)

# Load the pre-calculated TO embeddings. 
# You can re-calculate them using the 'embed_TOterms' function
load("../embeddings/TOterms_embedding.Rdata")





################################################
#
#     DCRAG workflow                       #####
#
################################################


#### load phenotype descriptors to be annotated ####

TAIR <- data.table::fread("../descriptors/TAIR_pheno_TO_gold100_labeled.txt",
                          header = TRUE, sep="\t", na.strings = "NA")
TAIR <- TAIR %>% 
  group_by(PHENOTYPE) %>% 
  mutate(item_id = cur_group_id(), .before = 1) %>% 
  ungroup()

TAIR.unique <- TAIR %>% 
  distinct(item_id, PHENOTYPE)




#### parse descriptors into concepts using LLM ####

# break descriptors into concepts using LLM
prompt1 <- "I will provide a plant phenotype description. It starts with a 
gene id then ':'.  You will extract the key concepts from the phenotype by 
understanding the meaning of the entire description. Separate compound concepts 
and discard concepts that say normal or unchanged from wildtype (WT). 
For each concept (C) provide three short phrases (HL1,HL2,HL3) that 
semantically explain C at a higher level, mainly using words and phrases found 
in Plant Trait Ontology labels and descriptions. There may be abbreviations for 
chemical compounds expand these where possible. e.g. SA is short for 
Salicylic Acid. Format the response as a tab-delimited file where each concept 
is a new row with columns Gene,C,HL1,HL2,HL3. If there is no concept then 
provide 'NA'. Do not output a header row."

# load the pre-parsed concepts of the TAIR descriptors
# If you want to re-parse these using GPT-4o then use the 'openai_getconcepts' function 
# which will call OpenAI and save the concepts to a file for use here.
concepts <- data.table::fread("../outputs/TAIR.concepts.2024-07-11.txt", sep="\t", header=F,
                              col.names = c("item_id", "concept","phrase1","phrase2","phrase3"), fill = TRUE)
concepts <- concepts %>% mutate(embedtext = paste(concept, phrase1, phrase2, phrase3, sep=", ")) %>% 
  mutate(concept_id = row_number(), .after=1)




#### load embeddings of full-length TAIR descriptors ####

# You can re-calculate these embeddings using the 'embed_descriptors' function
load("../embeddings/TAIR.2024-06-24_embedding.Rdata")


#### load embeddings of parsed concepts ####

# You can re-calculate these embeddings using the 'embed_concepts' function
load("../embeddings/TAIR.concepts.2024-07-11_embedding.Rdata")





#### use embedding similarity to find candidate TO terms ####
# here we use cosine similarity between the embedding of the descriptor 
# and the embedding of each TO term


# for each full descriptor, get the most similar TO terms
NTOP = 4   
best_per_item <- foreach( emb = iterators::iter(TAIR.embed$data$embedding), 
                          i   = iterators::icount()) %do% {
                            
                            cossim <- proxy::simil(x = list(emb = emb), 
                                                   y = TOterms$embedding, 
                                                   method = "cosine")
                            colnames(cossim) <- TOterms$ID
                            df <- data.frame(item_id = TAIR.unique$item_id[i], 
                                             Term    = colnames(cossim)[order(-cossim)[1:NTOP]], 
                                             cossim  = cossim[order(-cossim)[1:NTOP]])
                            df <- left_join(df, dplyr::select(TOterms, ID, bigstring), by=c("Term"="ID"))
                            df
                            
                          }
best_per_item <- data.table::rbindlist(best_per_item) %>% 
  filter(cossim >= 0.35)


# here we use cosine similarity between the embedding of each concept  
# and the embedding of each TO term

# for each concept, get the most similar TO terms
NTOP = 4  # get best terms per concept (remember, there are usually several concepts per item)
best_per_concept <- foreach( emb = iterators::iter(TAIR.concept.embed$data$embedding), 
                             i   = iterators::icount()) %do% {
                               cossim <- proxy::simil(x = list(emb = emb), 
                                                      y = TOterms$embedding, 
                                                      method = "cosine")
                               colnames(cossim) <- TOterms$ID
                               df <- data.frame(item_id    = concepts$item_id[i], 
                                                concept_id = concepts$concept_id[i],
                                                concept    = concepts$concept[i], 
                                                Term       = colnames(cossim)[order(-cossim)[1:NTOP]], 
                                                cossim     = cossim[order(-cossim)[1:NTOP]])
                               df <- left_join(df, dplyr::select(TOterms, ID, bigstring), by=c("Term"="ID"))
                               df
                             }
best_per_concept <- data.table::rbindlist(best_per_concept) %>% 
  filter(cossim >= 0.35) %>% 
  group_by(item_id) %>% 
  distinct(Term, .keep_all = TRUE) %>% 
  ungroup()



#### ask the LLM to filter the candidate TO terms ####

# Base Prompt explaining the task to the LLM

baseprompt <- "You are a plant biologist. You will be given a description (D) 
that describes observations of how a mutant plant was altered due to a treatment. 
D starts with an id (ID) then ':'. You will use your doctorate-level plant 
biology knowledge to annotate D with the most appropriate Plant Trait Ontology 
terms from the list below. Trait ontology terms start with 'TO:'. Step One is 
to find each of the phenotypic observations in D that occurred as a result of 
the treatment, but do not annotate the treatment itself. Step two is to use 
step-by-step reasoning and your understanding of plant experimental biology 
to relate the phenotypic observations to ontology terms based on their 
descriptors. Each phenotypic observation might match well with multiple 
ontology terms so report the matches that make the most sense based on plant 
biology and anatomy. Provide as many terms as necessary to address all relevant 
observations in D. If D states that no difference was observed, then do not 
provide ontology terms. An example of this is 'mutants showed no difference to 
wildtype'.\n\n Think step-by-step with reasoning about how each ontology term 
could apply to observations in D but do not show your reasoning for choosing a 
term. \n\n Format the response as a tab-delimited file where each ontology term 
is a new row with columns ID,term,label. Do not print out column headings. If 
there are no good terms then print 'NA'. \n\n"

# if you don't want to run the RAG, you can just load the pre-calculated 
# outputs (see fruther below)

# for each item (i.e. unique descriptor) ...
for(item in TAIR.unique$item_id)
{
  # combine the best TO terms for the descriptor and its concepts
  df <- rbind( best_per_concept %>% 
                 filter(item_id==item) %>% 
                 dplyr::select(Term, "description"=bigstring),
               best_per_item %>% 
                 filter(item_id==item) %>% 
                 dplyr::select(Term, "description"=bigstring)
  ) %>% distinct(Term, .keep_all = TRUE)
  
  # Augment the base prompt with candidate TO terms for this item (descriptor)
  prompt <- c(baseprompt, "```\n", knitr::kable(df, format = "simple" ), "```")
  prompt <- paste(prompt, collapse = "\n")
  cat(prompt)
  
  # ask the LLM to filter the candidate terms to a small subset
  openai_filterterms(TAIR.unique %>% 
                    filter(item_id==item) %>% 
                    select(item_id, descriptor=PHENOTYPE),
                  prompt = prompt, 
                  fname = "TAIR.RAG")
}


####  calculate score for each item  ####

# load the outputs from the RAG
TAIR.RAG <- data.table::fread("../outputs/TAIR.RAG.concepts_items.2024-07-23.txt", col.names = c("item_id","Term", "label"))
TAIR.RAG <- TAIR.RAG %>% distinct(item_id, Term, .keep_all = TRUE)
tmp <- left_join(TAIR, TAIR.RAG, by="item_id", suffix = c(".gold",".auto"))

score <- list()
score[["concepts_RAG"]] <- getscore(annotation = tmp, 
                                    itemscol  = "item_id", 
                                    preproc   = "concepts", 
                                    annotator = "RAG", 
                                    dataset   = "TAIR", 
                                    ontology  = TO, 
                                    ic        = TO.ic)



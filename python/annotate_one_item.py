################################################################################
#
# AI-driven workflows for auto-annotation of phenotypes.
# This code underpins the analysis in the manuscript titled
# "The effectiveness of Large Language Models with RAG for 
# auto-annotating phenotype descriptions"
#
# An attempt to turn this into Python code by Thomas. Probably not worth the effort.
#
#
# Author: David Kainer, Thomas Crow
# Contact: d.kainer@uq.edu.au, t.crow@uq.edu.au
#
################################################################################

import pandas as pd
import numpy as np
#import polars as pl # dataframe processing
#import openai
from owlready2 import * # Alternative for ontologyIndex
#import datatable as dt # Alternative for R's data.table


#library(tidyverse)
#library(openai)
#library(ontologyIndex)
#library(ontologySimilarity)
#library(data.table)


#########################################
#
# initialise 
#
########################################

# Need to run this line to access OpenAI API
# If you don't have a key, you need to register an account with OpenAI
#openai.api_key = "put your key here" # Maybe make variable another file?

# load ontology graph
ont_file = "/Users/thomascrow/PycharmProjects/PhD/LLMannotator/ontology/"
onto_path.append(ont_file)
ontology = get_ontology("http://purl.obolibrary.org/obo/to.owl").load()
#ontology = get_ontology("http://purl.obolibrary.org/obo/to/imports/chebi_import.owl").load()

print(list(ontology.classes()))

#TO    <- ontologyIndex::get_ontology(file = "ontology/to.obo")

# get information content for semsim
#TO.ic <- ontologySimilarity::descendants_IC(TO)


# Load the pre-calculated TO embeddings. 
# You can re-calculate them using the 'embed_TOterms' function 
#load("embeddings/TOterms_embedding.Rdata")

"""

#######################################3
#
# annotation!
#
########################################

# example phenotype descriptor to annotate
item <- "AT1G0000001: ABA hypersensitivity of guard cell anion-channel activation and stomal closing"




#### parse concepts with LLM ####

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

# set up query to LLM
messages = list(
  list("role"="system", "content"=prompt1),
  list("role"="user", "content"=item)
)
out <- openai::create_chat_completion(model = "gpt-4o",
                                      temperature = 0,
                                      max_tokens = 4000,
                                      messages = messages)

# prompt1 was designed for high-throughput usage in the study.
# i.e. The GPT output is formatted for writing to file.
# So for this single item demo we do a bit of editing of the output format
concepts <- out$choices$message.content
concepts <- strsplit(concepts,"\n" )
concepts <- lapply(concepts, str_replace_all, "\t", ",")
concepts <- lapply(concepts[[1]], FUN = function(x){ str_split(x, pattern = ",", n = 2)[[1]][2] })



#### get embeddings ####

# remove the gene ID before embedding the whole descriptor
embed.item <- openai::create_embedding(model = "text-embedding-3-large",
                                  input = utf8::utf8_encode( gsub(".*:","",item) ))

# this will embed each concept separately using one call to the LLM
embed.concepts <- openai::create_embedding(model = "text-embedding-3-large", 
                                          input = utf8::utf8_encode(unlist(concepts)))





#######################
#### DE annotation ####
#######################
NTOP = 4   

# similarity between the descriptor embedding and each TO term embedding
cossim <- proxy::simil(x = list(emb = embed.item$data$embedding[[1]]), 
                       y = TOterms$embedding, 
                       method = "cosine")

colnames(cossim) <- TOterms$ID
DEterms <- data.frame(
                 Term    = colnames(cossim)[order(-cossim)[1:NTOP]], 
                 cossim  = cossim[order(-cossim)[1:NTOP]]
                 )
DEterms <- left_join(DEterms, dplyr::select(TOterms, ID, bigstring), by=c("Term"="ID"))
DEterms
                            
# if any of the NTOP terms have cosine < 0.35 then get rid of them
DEterms <- DEterms %>%  filter(cossim >= 0.35)




########################
#### DCE annotation ####
########################

#concepts <- concepts %>% mutate(concept_id = row_number(), .after=1)
names(concepts) <- 1:length(concepts)

# for each concept, get the most similar TO terms
NTOP = 4  # get best terms per concept (remember, there are usually several concepts per item)
DCEterms <- foreach( emb = iterators::iter(embed.concepts$data$embedding), 
                             i   = iterators::icount()) %do% {
                               cossim <- proxy::simil(x = list(emb = emb), 
                                                      y = TOterms$embedding, 
                                                      method = "cosine")
                               colnames(cossim) <- TOterms$ID
                               df <- data.frame(concept_id = names(concepts)[i],
                                                concept    = concepts[[i]], 
                                                Term       = colnames(cossim)[order(-cossim)[1:NTOP]], 
                                                cossim     = cossim[order(-cossim)[1:NTOP]])
                               df <- left_join(df, dplyr::select(TOterms, ID, bigstring), by=c("Term"="ID"))
                               df
                             }
DCEterms <- data.table::rbindlist(DCEterms) %>% 
  filter(cossim >= 0.35) %>% 
  distinct(Term, .keep_all = TRUE) %>% 
  ungroup()




##########################
#### DCRAG annotation ####
##########################

prompt2 <- "You are a plant biologist. You will be given a description (D) 
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

  # combine the best TO terms for the descriptor and its concepts
  df <- rbind( DCEterms %>% 
                 dplyr::select(Term, "description"=bigstring),
               DEterms %>% 
                 dplyr::select(Term, "description"=bigstring)
  ) %>% distinct(Term, .keep_all = TRUE)
  
  # Augment the base prompt with candidate TO terms for this item (descriptor)
  prompt <- c(prompt2, "```\n", knitr::kable(df, format = "simple" ), "```")
  prompt <- paste(prompt, collapse = "\n")
  cat(prompt)
  
  # ask the LLM to filter the candidate terms to a small subset
  messages = list(
    list("role"="system", "content"=utf8::utf8_encode(prompt)),
    list("role"="user", "content"= utf8::utf8_encode(item))
  )
  #print(messages)            
  out <- openai::create_chat_completion(model = "gpt-4o",
                                        temperature = 0,
                                        max_tokens = 4000,
                                        messages = messages)

  tmp <- str_remove_all(out$choices$message.content, "```") %>% 
    str_split(.,pattern = "\n")

  DCRAGterms <- read.table(text = tmp[[1]], sep = "\t", 
                   col.names = c("id", "Term", "Label"))

"""
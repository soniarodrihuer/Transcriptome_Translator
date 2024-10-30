########################  TRANSCRIPTOME TRANSLATION #########################
########################          FUNCTION          #########################

# General information:

#### This code enables the translation of a .fasta file containing cDNA sequences derived from RNAseq.

#### The R version of this script is R 4.4.1

#### The only input needed is the fasta file path 

Transcriptome_Translator_Function <- function(route) {
  library(Biostrings)
  library(stringi)
  library(stringr)
  library(knitr)
  library(kableExtra)
  library(dplyr)
  
  
  # 1. Table with the sequences with the original labels
  # First, you have to specify a fasta. 
  
  # 1.1 Process the fasta in order to have a table with the labels and sequence
  seq <-readBStringSet(route)
  seq_strings <- toString(seq)
  Sense <- unlist(strsplit(seq_strings, split=", "))
  Label <- names(seq)
  seq_table <- data.frame(Label, Sense)
  
  
  # 2. Remove T with U (in order to have RNAm)
  # If you start with an RNA fasta you can skip this step.
  
  # 2.1 Remove T with U
  seq_table$Sense <- gsub("T", "U", seq_table$Sense)
  
  
  
  # 3. Reverse the sequences
  # As we want to observe the reverse strand reads, with this step, we reverse it (reading from right to left).
  
  # 3.1 Obtain the Reverse sequence
  seq_table$Reverse <- stri_reverse(seq_table$Sense)
  seq_table$Reverse <- chartr("UAGC", "AUCG", seq_table$Reverse)
  
  
  
  # 4. Obtain 6 orfs
  # Now the 6 possible ORFs will be generated for each sequence.
  
  # 4.1 Function for subsrings
  substrings_orf <- function(string_for_orf) {
    substrings_for_orf <- vector("character", length = 3)
    for (i in 1:3) {
      substrings_for_orf[i] <- substr(string_for_orf, i, nchar(string_for_orf))
    }
    return(substrings_for_orf)
  }
  
  # 4.2 Apply the function substrings_orf
  Seqs_Sense <- lapply(seq_table$Sense, substrings_orf)
  Seqs_Reverse <- lapply(seq_table$Reverse, substrings_orf)
  
  # 4.3 Obtain a new dataframe with the orfs
  sixorfs <- data.frame(
    ID = rep(seq_table$Label, each = 1),
    Sense_1 = as.vector(sapply(Seqs_Sense, function(x) ifelse(length(x) >= 1, x[1], NA))),
    Sense_2 = as.vector(sapply(Seqs_Sense, function(x) ifelse(length(x) >= 2, x[2], NA))),
    Sense_3 = as.vector(sapply(Seqs_Sense, function(x) ifelse(length(x) >= 3, x[3], NA))),
    Reverse_1 = as.vector(sapply(Seqs_Reverse, function(x) ifelse(length(x) >= 1, x[1], NA))),
    Reverse_2 = as.vector(sapply(Seqs_Reverse, function(x) ifelse(length(x) >= 2, x[2], NA))),
    Reverse_3 = as.vector(sapply(Seqs_Reverse, function(x) ifelse(length(x) >= 3, x[3], NA)))
  )
  
  
  
  # 5. Translate the orfs
  # Now the 6 ORFs are translated using the standard genetic code. 
  # It can be modified if necessary.
  
  # 5.1 Translation code
  translate_codon <- function(codon) {
    codon_table <- list(
      UUU = "F", UUC = "F", UUA = "L", UUG = "L",
      CUU = "L", CUC = "L", CUA = "L", CUG = "L",
      AUU = "I", AUC = "I", AUA = "I", AUG = "M",
      GUU = "V", GUC = "V", GUA = "V", GUG = "V",
      UCU = "S", UCC = "S", UCA = "S", UCG = "S",
      CCU = "P", CCC = "P", CCA = "P", CCG = "P",
      ACU = "T", ACC = "T", ACA = "T", ACG = "T",
      GCU = "A", GCC = "A", GCA = "A", GCG = "A",
      UAU = "Y", UAC = "Y", UAA = "X", UAG = "X",
      CAU = "H", CAC = "H", CAA = "Q", CAG = "Q",
      AAU = "N", AAC = "N", AAA = "K", AAG = "K",
      GAU = "D", GAC = "D", GAA = "E", GAG = "E",
      UGU = "C", UGC = "C", UGA = "X", UGG = "W",
      CGU = "R", CGC = "R", CGA = "R", CGG = "R",
      AGU = "S", AGC = "S", AGA = "R", AGG = "R",
      GGU = "G", GGC = "G", GGA = "G", GGG = "G"
    )
    
    return(codon_table[[codon]])
  }
  
  # 5.2 Obtain the sequences
  sequences <- sixorfs[, -1]      	# Exclude the Labels
  labels <- sixorfs[, 1]          	# Obtain the Labels
  translated_sequences <- matrix("", nrow = nrow(sequences), 
                                 ncol = ncol(sequences)) # Set the table
  
  # 5.3 Do the translation
  for (i in 1:nrow(sequences)) {
    for (j in 1:ncol(sequences)) {
      seq <- as.character(sequences[i, j])
      codon_number <- floor(nchar(seq) / 3)
      output <- character(0)
      for (k in 1:codon_number) {
        chunk <- substr(seq, 3 * (k - 1) + 1, 3 * k)
        translated <- translate_codon(chunk)
        output <- c(output, translated)
      }
      translated_sequences[i, j] <- paste(output, collapse = "")
    }
  }
  
  # 5.3 Obtain the table with the 6 translations of each sequence
  translated_table <- data.frame(Labels = labels,
                                 Sense1 = translated_sequences[, 1],
                                 Sense2 = translated_sequences[, 2],
                                 Sense3 = translated_sequences[, 3],
                                 Reverse1 = translated_sequences[, 4],
                                 Reverse2 = translated_sequences[, 5],
                                 Reverse3 = translated_sequences[, 6])
  
  
  
  
  # 6. Select the sequence with more letters between an M and X (with no more X)
  # This way, the longest possible sequences are selected within each sequence. 
  # - IT DOES NOT CONSIDER THE CASE WHERE 2 SEQUENCES ARE EQUAL IN LENGTH.
  # - IT CHOOSES THE LAST ONE OF GREATER LENGTH.
  
  # 6.1 Function to detect the longest seq between M and X (without no more X)
  # If there are no M a 0 is written
  find_subseq <- function(sequence) {
    matches <- str_extract_all(sequence, "M[^X]*X")[[1]]
    if (length(matches) > 0) {
      longest_match <- which.max(nchar(matches))
      return(matches[longest_match])
    } else {
      return(0)
    }
  }
  
  # 6.2 Apply function to each seq
  subseq_sense1 <- apply(translated_table[, "Sense1", drop = FALSE], 1, find_subseq)
  subseq_sense2 <- apply(translated_table[, "Sense2", drop = FALSE], 1, find_subseq)
  subseq_sense3 <- apply(translated_table[, "Sense3", drop = FALSE], 1, find_subseq)
  subseq_reverse1 <- apply(translated_table[, "Reverse1", drop = FALSE], 1, find_subseq)
  subseq_reverse2 <- apply(translated_table[, "Reverse2", drop = FALSE], 1, find_subseq)
  subseq_reverse3 <- apply(translated_table[, "Reverse3", drop = FALSE], 1, find_subseq)
  
  # 6.3 Create a new table with the results of the longest translation
  Longest_of_each_translation <- data.frame(
    Label = seq_table$Label,
    Sense_1 = subseq_sense1,
    Sense_2 = subseq_sense2,
    Sense_3 = subseq_sense3,
    Reverse_1 = subseq_reverse1,
    Reverse_2 = subseq_reverse2,
    Reverse_3 = subseq_reverse3
  )
  
  
  
  # 7. Select the longest for each reading
  # In this part the longest of the six possibilities is selected as the protein for the original sequence.
  
  
  # 7.1 Function to find the longest translation
  find_the_longest_of_the_longest <- function(sequences) {
    seq_length <- sapply(sequences, nchar)
    highest_number <- which.max(seq_length)
    return(sequences[highest_number])
  }
  
  # 7.2 Apply to the data and obtain the final result
  Longest_Protein <- apply(Longest_of_each_translation[, -1, drop = FALSE], 1,
                           find_the_longest_of_the_longest)
  RESULT <- data.frame(Label = seq_table$Label,
                       Longest_Protein = Longest_Protein)
  RESULT$Longest_Protein <- gsub("X", "", RESULT$Longest_Protein)
  
  # 7.3 Display the table with data
  RESULT %>%
    mutate(Longest_Protein = str_c(str_sub(Longest_Protein, end = 10), "...")) %>%
    kable() %>%
    kable_styling(full_width = FALSE, position = "center", 
                  font_size = 14, latex_options = c("striped", "hold_position"))
  
  # 7.4 Obtain a .fasta with the results
  Write_fasta <- function(df, file_name) {
    fileConn <- file(file_name)
    fasta_lines <- unlist(lapply(1:nrow(df), function(i) {
      c(paste0(">", df$Label[i]), df$Longest_Protein[i])
    }))
    writeLines(fasta_lines, fileConn)
    close(fileConn)
  }
  
  Write_fasta(RESULT, "RESULT.fasta")
  
  print("Done. Congratulations! You have your results downloaded as a file called <RESULT.fasta>")
  
}


########################           EXAMPLE           #########################

# Example of use
Transcriptome_Translator_Function("C:/Users/asus/Desktop/RGit/transcriptome_fern/General code and trials/Sample_for_Transcriptome_Translator.fasta")
          # or
        actual_directory <- getwd()
        path <- file.path(actual_directory, "Sample_for_Transcriptome_Translator.fasta")
        Transcriptome_Translator_Function(path)
        
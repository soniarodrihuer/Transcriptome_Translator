library(shiny)
library(bslib)
library(tidyverse)

library(Biostrings)
library(stringi)
library(stringr)

library(knitr)
library(kableExtra)
library(dplyr)




# Define UI for data upload and download app ----
ui <- fluidPage(
    theme = bs_theme(version = 4, bootswatch = "spacelab"),
  
  # App title ----
  titlePanel("TRANSCRIPTOME TRANSLATOR"),
  
  # Panel Title ----
  navlistPanel(
  
  
    
  # First Panel ----
  tabPanel("Translate your FASTA",
           
           h4(strong("LET'S TRANSLATE A TRANSCRIPTOME!")),
           p(strong("1. Upload the FASTA file with the Browse button.")),
              p("Remember: 10,000 sequences maximum! "),
           p(strong("2. Click on <Download your translated transcriptome> and wait.")),
              p("For 10,000 sequences it should be around 2 minutes maximum."),
           p(strong("3. Enjoy your translated transcriptome as a .csv or as .fasta!")),
              p("For more information about the process that has been carried out for the translation, you can refer to the tab 'Behind-the-scenes process'."),
           p("If needed, there is a FASTA file example in <Behind the scenes process>."),
           
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      width=8,
      
      # Input: Select a file ----
      fileInput("file1", "Please, upload your transcriptomic .fasta file here:",
                accept = c(".fasta")),
      
      # Horizontal line ----
      tags$hr(),
      
      
      
      # Button for downloading translated transcriptome ----
      downloadButton("downloadCSV", "CSV: Download your translated transcriptome"),
      HTML("<br><br>"),
      downloadButton("downloadFASTA", "FASTA: Download your translated transcriptome")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents"),
      
      # Output: Table of translated transcriptome ----
      tableOutput("table")
      
    )
    
  )
),

    # Second Panel ----
    tabPanel("Behind-the-scenes process",
             h4(strong("Process behind:")),
             h6(strong("1. Make a table with the sequences with the original labels")),
             p("First, you are going to upload a FASTA file with the sequences (whether they are mRNA or cDNA). 
                            In this initial step, the FASTA file is converted into a table to facilitate the translation process in the following steps."),
              HTML("<p>You can find a sample <a href='https://github.com/soniarodrihuer/Transcriptome_Translator'>here</a>!</p>"),
            
             downloadButton("downloadFile", "Download here a Sample"),
             
              h6(strong("2. The T's are replaced by U's.")),
             p("By changing T's to U's, it ensures having mRNA in case cDNA sequences were introduced."),
             h6(strong("3. The reverse sequence is obtained.")),
             p("The reading of the sequence from right to left is added, having changed U to A, A to U, G to C, and C to G."),
             h6(strong("4. The 6 open reading frames (ORFs) are obtained.")),
             p("The three possible readings for each direct sequence and three readings for each reverse sequence are obtained."),
             h6(strong("5. The possible translations for each open reading frame (ORF) are obtained.")),
             p("For this step, the standard genetic code is used."),
             h6(strong("6. In each of the 6 readings, the sequence with the most letters between a Methionine (M) and a Stop (X) is selected, without any intervening X codons. ")),
             p("This way, the longest possible sequences are selected within each sequence.
                            It is important to highlight that in case two sequences have the same length, the last one read with the greater length will be retained."),
             h6(strong("7. The longest protein sequence is selected from all the possibilities, and a .csv file is generated with the final result for each sequence.")),
             p("In this part, the longest of the six possibilities is selected as the protein for the original sequence.")
    
   
    ),


    # Third Panel ----
    tabPanel("More about TRANSCRIPTOME TRANSLATOR APP",
             h4(strong("More about the app:")),
             p("- This application was created in 2024 by Sonia María Rodríguez Huerta and José Manuel Álvarez Díaz (Vegetal Physiology, University of Oviedo)."),
             p("- R 4.4.1 was used, and the required packages are shiny, bslib, Biostrings, stringi, stringr, knitr, kableExtra, and dplyr."),
             HTML("<p>The code is available in <a href='https://github.com/soniarodrihuer/Transcriptome_Translator'>this GitHub repository</a>!</p>")
    )

)
)





# Define server logic to read selected file and handle download ----
server <- function(input, output) {
  
  # Download Sample Fasta
  output$downloadFile <- downloadHandler(
    filename = function() {
      "Sample_for_Transcriptome_Translator.fasta"
    },
    content = function(file) {
      file.copy("Sample_for_Transcriptome_Translator.fasta", file)
    }
  )
  
  
  # 1. Table with the sequences with the original labels
  datasetInput <- reactive({
    req(input$file1)
    seq <- readDNAStringSet(input$file1$datapath)
    seq_strings <- toString(seq)
    Sense <- unlist(strsplit(seq_strings, split=", "))
    Label <- names(seq)
    seq_table <- data.frame(Label, Sense)
    
    # 2. Remove T with U (in order to have RNAm)
    seq_table$Sense <- gsub("T", "U", seq_table$Sense)
    
    # 3. Reverse the sequences
    seq_table$Reverse <- stri_reverse(seq_table$Sense)
    seq_table$Reverse <- chartr("UAGC", "AUCG", seq_table$Reverse)
    
    # 4. Obtain 6 orfs
    substrings_orf <- function(string_for_orf) {
      substrings_for_orf <- vector("character", length = 3)
      for (i in 1:3) {
        substrings_for_orf[i] <- substr(string_for_orf, i, nchar(string_for_orf))
      }
      return(substrings_for_orf)
    }
    
    Seqs_Sense <- lapply(seq_table$Sense, substrings_orf)
    Seqs_Reverse <- lapply(seq_table$Reverse, substrings_orf)
    
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
    
    sequences <- sixorfs[, -1]
    labels <- sixorfs[, 1]
    translated_sequences <- matrix("", nrow = nrow(sequences), ncol = ncol(sequences))
    
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
    
    translated_table <- data.frame(Labels = labels,
                                   Sense1 = translated_sequences[, 1],
                                   Sense2 = translated_sequences[, 2],
                                   Sense3 = translated_sequences[, 3],
                                   Reverse1 = translated_sequences[, 4],
                                   Reverse2 = translated_sequences[, 5],
                                   Reverse3 = translated_sequences[, 6])
    
    find_subseq <- function(sequence) {
      matches <- str_extract_all(sequence, "M[^X]*X")[[1]]
      if (length(matches) > 0) {
        longest_match <- which.max(nchar(matches))
        return(matches[longest_match])
      } else {
        return(0)
      }
    }
    
    subseq_sense1 <- apply(translated_table[, "Sense1", drop = FALSE], 1, find_subseq)
    subseq_sense2 <- apply(translated_table[, "Sense2", drop = FALSE], 1, find_subseq)
    subseq_sense3 <- apply(translated_table[, "Sense3", drop = FALSE], 1, find_subseq)
    subseq_reverse1 <- apply(translated_table[, "Reverse1", drop = FALSE], 1, find_subseq)
    subseq_reverse2 <- apply(translated_table[, "Reverse2", drop = FALSE], 1, find_subseq)
    subseq_reverse3 <- apply(translated_table[, "Reverse3", drop = FALSE], 1, find_subseq)
    
    Longest_of_each_translation <- data.frame(
      Label = seq_table$Label,
      Sense_1 = subseq_sense1,
      Sense_2 = subseq_sense2,
      Sense_3 = subseq_sense3,
      Reverse_1 = subseq_reverse1,
      Reverse_2 = subseq_reverse2,
      Reverse_3 = subseq_reverse3
    )
    
    find_the_longest_of_the_longest <- function(sequences) {
      seq_length <- sapply(sequences, nchar)
      highest_number <- which.max(seq_length)
      return(sequences[highest_number])
    }
    
    Longest_Protein <- apply(Longest_of_each_translation[, -1, drop = FALSE], 1,
                             find_the_longest_of_the_longest)
    RESULT <- data.frame(Label = seq_table$Label,
                         Longest_Protein = Longest_Protein)
    RESULT$Longest_Protein <- gsub("X", "", RESULT$Longest_Protein)
    
    return(RESULT)
  })
  
  # Función para escribir en formato FASTA
  Write_fasta <- function(df, file_name) {
    fileConn <- file(file_name)
    fasta_lines <- unlist(lapply(1:nrow(df), function(i) {
      c(paste0(">", df$Label[i]), df$Longest_Protein[i])
    }))
    writeLines(fasta_lines, fileConn)
    close(fileConn)
  }
  
  # Download Results as CSV
  output$downloadCSV <- downloadHandler(
    filename = function() {
      "Your_translated_transcriptome_as_csv.csv"
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
  
  # ownload Results as FASTA
  output$downloadFASTA <- downloadHandler(
    filename = function() {
      "Your_translated_transcriptome_as_fasta.fasta"
    },
    content = function(file) {
      Write_fasta(datasetInput(), file)
    }
  )
  

  
}

# Run the Shiny app ----
shinyApp(ui, server)


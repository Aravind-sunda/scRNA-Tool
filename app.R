#Libraries----
library(shiny)
library(shinydashboard)
library(Seurat)
library(tidyverse)
library(shinyjs)
library(BiocManager)
library(data.table)
library(fgsea)
library(DT)
library(shinycssloaders)
library(biomaRt)
library(plotly)
library(glue)

options(repos = BiocManager::repositories())


options(shiny.maxRequestSize = 1024*1024^2)  # 30 MB, adjust as needed
#functions are named with camelCase
#ids are named with snake_case

# Global Functions----
## Page 1- Clustering Functions----
plotUMAP <- function(object = NULL, group = NULL, reduction = NULL){
  plot <- DimPlot(object = object,reduction = reduction ,label = TRUE,group.by = group)
  return(plot)
}

generateClusterCheckboxes <- function(object = NULL,group = NULL){
  Idents(object) <- group
  levels <- levels(object)
  return(levels)
}

umapSubcluster <- function(object = NULL, group = NULL, resolution = NULL,subclusters = NULL, subgroup = NULL) {
  Idents(object) <- as.character(group) 
  # object <- SetIdent(object = object,ident.use = group)
  
  # Subsetting the cluster
  subset_obj <- subset(object, idents = subclusters)
  
  #Rerunning the UMAP clustering again
  subset_obj <- RunPCA(subset_obj, verbose = FALSE)
  subset_obj <- FindNeighbors(subset_obj, dims = 1:10)
  subset_obj <- FindClusters(subset_obj, resolution = resolution)
  
  # Plotting UMAP
  plot1 <- DimPlot(subset_obj, reduction = "umap",label = TRUE,pt.size = 0.01)
  
  #Grouping by CJD13/CJD13CR
  plot2 <- DimPlot(subset_obj, reduction = "umap", group.by = subgroup,pt.size = 0.01)
  
  return(list(plot1,plot2))
}

## Page 2-GSEA functions----
rankgenes <-  function(data = NULL, rankingmetric = NULL,filter = FALSE,padj_filter = NULL){
  # Changed due to speed restrictions
  if ("avg_log2FC" %in% colnames(data)) {
    DEG <- data %>%
      na.omit() %>%
      arrange(desc(avg_log2FC))
  } else {
    DEG <- data %>%
      na.omit()
  }
  
  # distinct() %>% 
  #   group_by(gene) %>% 
  #   summarize(stat=mean(stat))%>%
  
  if(filter == TRUE){
    DEG <- DEG %>%
      arrange(desc(avg_log2FC)) %>% 
      subset(.,p_val_adj < padj_filter)  
    
  }
  
  if (rankingmetric == "1"){
    DEG <- DEG %>%
      mutate(sort_metric = -log10(p_val)) %>%
      dplyr::arrange(desc(sort_metric))
  }
  
  if(rankingmetric == "2"){
    DEG <- DEG %>%
      mutate(sort_metric = sign(avg_log2FC) * -log10(p_val)) %>%
      dplyr::arrange(desc(sort_metric))
  }
  
  if (rankingmetric == "3"){
    DEG <- DEG %>%
      mutate(sort_metric = avg_log2FC * -log10(p_val)) %>%
      dplyr::arrange(desc(sort_metric))
  }
  
  if (rankingmetric == "4"){
    DEG <- DEG %>%
      mutate(sort_metric = avg_log2FC) %>%
      dplyr::arrange(desc(sort_metric))
  }
  
  if (rankingmetric == "5") {
    DEG <- DEG %>%
      mutate(sort_metric = stat) %>%
      dplyr::arrange(desc(sort_metric))
  }
  if(rankingmetric =="6"){
    DEG <- DEG %>%
      mutate(sort_metric = rank) %>%
      dplyr::arrange(desc(sort_metric))
  }
  
  stat <- DEG$sort_metric
  names(stat) <- DEG$gene
  gc()
  return(stat)
  
}


GSEA_function <- function(dataset_path = NULL, namedlist = NULL, padj_filter = NULL, title = "") {
  set.seed(123)
  dataset <- gmtPathways(dataset_path)
  result <- fgsea::fgsea(pathways = dataset, stats = namedlist)
  gc()
  
  significant_results <- subset(result, padj < padj_filter)%>%
    arrange(desc(NES)) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", ")))
  
  # Plotting the significant results
  plot <- ggplot(significant_results, aes(x=reorder(pathway, NES), y=NES, fill=(NES > 0))) +
    geom_bar(stat="identity") +
    coord_flip() +
    labs(title = title, x = "Pathway", y = "NES",fill = "Enrichment",caption = paste("All pathways have adjusted P-values(padj) less than",padj_filter)) +
    scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "darkorange"), labels = c("TRUE"="Upregulation", "FALSE" = "Downregulation")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold.italic", size = 12))
  # Creating a more informative image in the form of a dot plot
  results <- significant_results %>%
    mutate(regulation = ifelse(NES > 0, "Up", "Down"))
  
  dotplot <-  ggplot(results, aes(x = reorder(pathway, NES), y = NES, color = regulation,size = -log10(padj))) +
    geom_point() +
    geom_text(aes(label = sprintf("%.3f", padj), hjust = ifelse(regulation == "Up", 1.5, -0.5)), vjust = -0, size = 3) + # This can be removed to remove the annotations in the image
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", color = "Regulation") +
    theme_classic() +
    coord_flip()
  
  return(list(significant_results = significant_results,filteredplot = plot, dotplot = dotplot))
}


convertion <- function(data = NULL, species = NULL, id_type = NULL){
  if (id_type == "name") {
    return(data)
  }
  
  else {
    if (species == "mouse") {
      ds <- "mmusculus_gene_ensembl"
    }
    if (species == "human") {
      ds <- "hsapiens_gene_ensembl"
    }
    if (id_type == "entrez") {
      id = "entrezgene_id"
    }
    if (id_type == "ensembl") {
      id = "ensembl_gene_id"
    }
    list <- data$gene
    ensembl <- useEnsembl(biomart = "genes", dataset = ds,mirror = 'useast')
    mapping <- getBM(attributes = c(id,"external_gene_name"),mart = ensembl,filters = id,values = list)
    
    result <- data %>%
      left_join(mapping,by = c("gene" = id),unmatched = "drop") %>%
      mutate(gene = external_gene_name) %>%
      dplyr::select(-external_gene_name) %>%
      filter(complete.cases(.))
    
    return(result)
  }
}

## Page 3- Volcano Plot Functions----
load_data <- function(path){
  file <- fread(file = path) #%>%
    # dplyr::rename("genes" = "...1")
  return(file)
}

draw_table <- function(dataf, slider) {
  filtered_table <- dataf %>% dplyr::filter(p_val_adj <= 1*10^slider) %>%
    drop_na()
  filtered_table$padj <- format(filtered_table$p_val_adj, digits = 6)
  filtered_table$pvalue <- format(filtered_table$p_val, digits = 6)
  return(filtered_table)
}

volcano_plot <-
  function(dataf, x_name, y_name, slider, color1, color2) {
    dataf <- dataf %>% drop_na()
    plot <- ggplot(data = dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)),text = gene))+
      geom_point(aes(color = if_else(!!sym(y_name) < 1*10^slider, TRUE,FALSE)))+
      labs(title="DESeq Results", x = glue('{x_name}'), y = glue('-log10({y_name})')) +
      scale_color_manual(name = glue("{x_name} < 1*10^{slider}"),
                         values = c("FALSE" = color1,
                                    "TRUE" = color2),
                         labels = c("FALSE", "TRUE"))+
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    
    return(plot)
  }


# Variables----
dataset_list <- list(
  human = c("Hallmark" = "h.all.v2023.1.Hs.symbols.gmt",
            "Chromosome Positional" = "c1.all.v2023.1.Hs.symbols.gmt",
            "Chemical Genetic Perturbation" = "c2.cgp.v2023.1.Hs.symbols.gmt",
            "BioCarta" = "c2.cp.biocarta.v2023.1.Hs.symbols.gmt",
            "KEGG" = "c2.cp.kegg.v2023.1.Hs.symbols.gmt",
            "PID" = "c2.cp.pid.v2023.1.Hs.symbols.gmt",
            "Reactome" = "c2.cp.reactome.v2023.1.Hs.symbols.gmt",
            "WikiPathways" = "c2.cp.wikipathways.v2023.1.Hs.symbols.gmt",
            "MicroRNA Targets" = "c3.mir.v2023.1.Hs.symbols.gmt",
            "Transcription Factor Targets" = "c3.tft.v2023.1.Hs.symbols.gmt",
            "GO Gene Ontology" = "c5.go.v2023.1.Hs.symbols.gmt",
            "GO Gene Ontology:BP" = "c5.go.bp.v2023.1.Hs.symbols.gmt",
            "GO Gene Ontology:CC" = "c5.go.cc.v2023.1.Hs.symbols.gmt",
            "GO Gene Ontology:MF" = "c5.go.mf.v2023.1.Hs.symbols.gmt",
            "Human Phenotype Ontology" = "c5.hpo.v2023.1.Hs.symbols.gmt",
            "ImmuneSigDB" = "c7.immunesigdb.v2023.1.Hs.symbols.gmt",
            "Vaccine Response" = "c7.vax.v2023.1.Hs.symbols.gmt",
            "Cell Type Signature" = "c8.all.v2023.1.Hs.symbols.gmt"),
  mouse = c("Hallmark" = "mh.all.v2023.1.Mm.symbols.gmt",
            "Chromosome Positional" = "m1.all.v2023.1.Mm.symbols.gmt",
            "Chemical Genetic Perturbation" = "m2.cgp.v2023.1.Mm.symbols.gmt",
            "Biocarta" = "m2.cp.biocarta.v2023.1.Mm.symbols.gmt",
            "Reactome" = "m2.cp.reactome.v2023.1.Mm.symbols.gmt",
            "Wikipathways" = "m2.cp.wikipathways.v2023.1.Mm.symbols.gmt",
            "MiRDB-miRNA pathways" = "m3.mirdb.v2023.1.Mm.symbols.gmt",
            "GTRD-TF binding Sites" = "m3.gtrd.v2023.1.Mm.symbols.gmt",
            "GO Gene Ontology" = "m5.go.v2023.1.Mm.symbols.gmt",
            "GO Gene Ontology:BP" = "m5.go.bp.v2023.1.Mm.symbols.gmt",
            "GO Gene Ontology:CC" = "m5.go.cc.v2023.1.Mm.symbols.gmt",
            "GO Gene Ontology:MF" = "m5.go.mf.v2023.1.Mm.symbols.gmt",
            "MPT - Tumour Phenotype Ontology" = "m5.mpt.v2023.1.Mm.symbols.gmt",
            "Cell Type Signature" = "m8.all.v2023.1.Mm.symbols.gmt"
  )
)


# UI----
## Header----
header <- dashboardHeader(title = "Bioinformatics Toolkit")

##Sidebar----
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("About",tabName = "About"),
    menuItem("Clustering",tabName = "Clustering"),
    menuItem("GSEA",tabName = "GSEA"),
    menuItem("Visualize DE results",tabName = "DE")
  )
)

##Body----
body <- dashboardBody(
  tabItems(
    ###About Page----
    tabItem(tabName = "About",
            h1("Bioinformatics Toolkit"),
            h3(strong("Credits")),
            p("Thanks to the Harris lab for providing the resources and feedback to make this application possible."),
            h3(strong("Summary")),
            p("This app focuses on downstream analysis of data generated for scRNA-seq and RNA-seq."),
            h3(strong("Overview")),
            p(tags$ol(
              tags$li("Clustering page: This page subclusters the idents of the UMAP given the scRNA-seq data object obtained from Seurat. The data object(RDS object) has to have been preprocessed and clusters identified for this page to work " ),
              tags$li("GSEA page: This page performs GSEA given the output list of differential expression. The input has to be the csv file output of differential expression with the column names clearly specified according to the requirements of the app. The output will be the Enrichment bar plot, datatable showing significant terms and the Enrichment dot plot"),
              tags$li("Differential Expression Data Analysis Page : This page visualizes the output of differential expression in a suitable format according to your requriements. The most common output that can be visualized is the volcano plot. The input has to be the csv file output of differential expression with the column names clearly specified according to the requirements of the app. The output will be the HTML plot with the significant genes that have been identified.")
            ))
    ),
    ###Clustering Page----
    tabItem(tabName = "Clustering",
            #box to load in the data
            fluidRow(
              box(title = "Input",width = 12,solidHeader = TRUE,collapsible = TRUE,
                  fileInput(inputId = "rdsobj","Select the Processed RDS object that you would like to analyse",)
              )
            ),
            
            fluidRow(
              box(title = "UMAP Clustering",width = 4,solidHeader = TRUE,collapsible = TRUE,
                  selectInput(inputId = "cluster_dropdown", label = "Select Option", choices = ""),
                  radioButtons(inputId = "dim_red_type",label = "Select reduction method",choices = c("UMAP" = "umap","T-SNE" = "tsne","PCA" = "pca")),
                  actionButton(inputId = "plot_umap", label  = "Generate UMAP")
              ),
              box(title = "UMAP Plot", width = 8,solidHeader = TRUE,collapsible = TRUE,
                  plotOutput("UMAP")
              )
            ),
            
            fluidRow(
              box(title = "UMAP Clustering",width = 12,solidHeader = TRUE,collapsible = TRUE,
                  checkboxGroupInput(inputId = "selected_clusters",label = "Select the clusters to be subclustered",choices = "Choices will be generated"),
                  selectInput(inputId = "cluster_dropdown_2", label = "Select Option", choices = "Choices will be generated"),
                  #cluster_dropdown2 has the same options as clusterdropdown1 but it is used for subclustering function
                  sliderInput(inputId = "umap_resolution",label = "Enter the clustering resolution for UMAp",min = 0,max = 2,step = 0.05,value = 0.25),
                  actionButton(inputId = "plot_subclustered_umap", label ="Generate Subclustered UMAPs")
              )
            ),
            
            fluidRow(
              box(title = "Subcluster UMAP plot", width = 6,solidHeader = TRUE,collapsible = TRUE,
                  plotOutput("subcluster1")
              ),
              box(title = "Subcluster UMAP plot-grouped", width = 6,solidHeader = TRUE,collapsible = TRUE,
                  plotOutput("subcluster2")
              )
            )
    ),
    
    ###GSEA Page----
    tabItem(tabName = "GSEA",
            fluidPage(
              sidebarLayout(
                sidebarPanel(
                  wellPanel(
                    fileInput("DEGsheet", "Upload the excel sheet of differentially expressed genes of the clusters to be checked", accept = c(".xlsx",".csv")),
                    helpText(HTML("<span style='color:red;'>Ensure that the table for the differentially expressed genes has the following columns named accordingly in the exactly manner given below:
                     <ul>
                        <li>gene: Gene Symbols, Entrez Ids or ENSEMBL ids(without version numbers)  being ranked. The application would convert it into the suitable common names for genes of Entrez Ids or ENSEMBL ids are given</li>
                        <li>p_val: P-values of the genes </li>
                        <li>avg_log2FC: Log 2 Fold Changes of each and every gene</li>
                        <li>p_val_adj: Adjusted p-values of the genes with a suiable correction method. Needed only if you need to filter padj values.</li>
                        <li>stat: Output of DESeq2, called the Wald statistic. It is described as the ratio of the estimated log2foldchange to the standard error of the estimated log2foldchange. If the absolute value of the 'stat' value is large, it shows that the alternative to the null hypothesis that the log2-fold change is equal to zero is true.The sign of the 'stat' value indicates the direction of the change. A positive value indicates upregulation, while a negative value indicates downregulation. Ordering the ranked gene list based on this statistic is particularly advantageous for GSEA.</li>
                        <li> Rank : Select 'Ranked or Unranked List(Gene Names + Ordering metric)' if you have only 2 columns of gene name and ranking metric. Second column has to be named 'rank'
                     </ul></span>"))
                  ),
                  radioButtons("genetype","Select the species below", choiceNames = c("Human","Mouse"), choiceValues = c("human","mouse"), selected = "Human",inline = TRUE),
                  radioButtons("geneidtype","Select one of the gene identifiers below", choiceNames = c("Gene Names","Entrez (NCBI) Gene ID","Ensembl Gene IDs"), choiceValues = c("name","entrez","ensembl"), selected = "Gene Names", inline = TRUE),
                  radioButtons("filterstatus", "Do you want to filter padj values for input DEG's?", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE", inline = TRUE),
                  conditionalPanel(
                    condition = "input.filterstatus == 'TRUE'",
                    sliderInput(inputId = "padjfilter",label = "Select the filter for the padj values for input DEG's",min = 0, max = 1,value = 1)
                  ),
                  radioButtons("rankingmetric","Select the ranking metric to use to rank the list of genes",choices = c("-log10(p_val)" = "1",
                                                                                                                        "sign of avg_log2FC * -log10(p_val)" = "2",
                                                                                                                        "avg_log2FC * -log10(p_val)" = "3",
                                                                                                                        "avg_log2FC" = "4",
                                                                                                                        "stat" = "5",
                                                                                                                        "Ranked or Unranked List(Gene Names + Ordering metric)" = "6"),selected = "2"),
                  uiOutput("dropdownUI"),
                  numericInput("GSEAfilternumeric", "Select the padj filter to apply for the output GSEA results:", value = 0.25, min = 0, max = 1,step = 1e-10),
                  actionButton(inputId = "DEGsubmit",label = "Submit")
                ),
                
                mainPanel(
                  tabsetPanel(type = "tabs",
                              tabPanel("GSEA Enrichment Plot",shiny::plotOutput("GSEAplot",width = "auto", height = "1200px"),
                                       numericInput("img_width", "Image width", 10),
                                       numericInput("img_height", "Image height", 10),
                                       numericInput("img_res", "Image resolution", 1000),
                                       downloadButton("downloadPlot", "Download GSEA Plot")),
                              tabPanel("GSEA table",
                                       DT::dataTableOutput("GSEAtable")),
                              tabPanel("Enrichment Dotplot",shiny::plotOutput("GSEAdotplot",width = "auto", height = "1200px"),
                                       numericInput("img_width_dp", "Image width", 10),
                                       numericInput("img_height_dp", "Image height", 10),
                                       numericInput("img_res_dp", "Image resolution", 1000),
                                       downloadButton("downloadPlot_dp", "Download GSEA Dot Plot"))
                  )
                )
              )
            )
    ),
    
    
    ###Volcano plot for differential expression page----
    tabItem(tabName = "DE",
            h1("DESeq Data Analysis"),
            h3(strong("Information about page:")),
            tags$p("This page takes in DESeq data and filters and displays the data in suitable ways"),
            fluidRow(
              box(title = "Input",width = 3,collapsible = TRUE,
                  fileInput("deseqfile","Choose File to read",accept = ".txt"),
                  HTML("<p>A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis.</p>"),
                  radioButtons("x_axis", "Choose the column for X-axis", choices=c("baseMean","avg_log2FC","lfcSE","stat","p_val","p_val_adj"),selected ="avg_log2FC" ),
                  radioButtons("y_axis", "Choose the column for Y-axis", choices=c("baseMean","avg_log2FC","lfcSE","stat","p_val","p_val_adj"),selected ="p_val_adj" ),
                  colourpicker::colourInput(inputId = "base", label = "Base Point Color:", value = "black"),
                  colourpicker::colourInput(inputId = "highlight", label = "Highlight Point Color:", value = "#F56767"),
                  sliderInput(inputId = "slider",label = "Select the magnitude of the p adjusted coloring:",min = -100 ,max = 0,value = -10, step = 1),
                  div(actionButton(inputId = "desubmit", label = "Plot and Filter", style = "color: white; background-color: #0072B2;"), style = "text-align: center;")
              ),
              conditionalPanel("input.desubmit",
                               tabBox(title = "Results of DESeq",width = 9,side = "right",selected = "Filtered Table",
                                      tabPanel("Volcano Plot",withSpinner(plotlyOutput("devolcano"))),
                                      tabPanel("Filtered Table",dataTableOutput("detable"))
                               )
              )
            )
    )
    
  )
)

##page
ui <-dashboardPage(header,sidebar,body)


#Server----
server <- function(input, output, session) {
  
## Logic for Clustering page----
  rds_obj <- reactive({
    req(input$rdsobj)
    file <- input$rdsobj$datapath
    readRDS(file)
  })
  
  options_data <- reactive({
    req(rds_obj())
    meta.data <- rds_obj()@meta.data
    new_options <- meta.data %>%
      keep(is.factor) %>%
      names()
    return(new_options)
  })
  
  # Update the dropdown menu choices based on the reactive expression
  observe({
    choices <- options_data()
    updateSelectInput(session, "cluster_dropdown", choices = c(choices))
    updateSelectInput(session,inputId = "cluster_dropdown_2",choices = c(choices))
  })
  
  umap <- eventReactive(input$plot_umap,{
    req(rds_obj(),input$cluster_dropdown,input$dim_red_type)
    plotUMAP(object = rds_obj(),group = input$cluster_dropdown ,reduction = input$dim_red_type)
  })
  
  output$UMAP <- renderPlot({
    umap()
  })
  
  observeEvent(input$plot_umap,{
    req(input$cluster_dropdown)
    options <- generateClusterCheckboxes(object = rds_obj(),group = input$cluster_dropdown)
    updateCheckboxGroupInput(inputId = "selected_clusters",choices = options)
    cat(options)
  })
  
  subcluster_umaps <-eventReactive(input$plot_subclustered_umap,{
    # req()#add required inputs here
    umaps <- umapSubcluster(object = rds_obj(),group = input$cluster_dropdown,resolution = input$umap_resolution,subclusters = input$selected_clusters,subgroup = input$cluster_dropdown_2)
    return(umaps)
  }) 
  
  output$subcluster1 <- renderPlot({subcluster_umaps()[1]})
  
  output$subcluster2 <- renderPlot({subcluster_umaps()[2]})
  
  
## Logic for GSEA page----
  uploaded_data <- reactive({
    req(input$DEGsheet,input$genetype ,input$geneidtype)
    gc()
    table <- data.table::fread(input$DEGsheet$datapath)
    convertion(data = table, species = input$genetype, id_type = input$geneidtype)
  })
  
  selected_option <- reactive({
    req(input$genetype)
    input$genetype
  })
  
  output$dropdownUI <- renderUI({
    # selectInput("dropdown", "Items:", choices = server_side_list[[selected_option()]])
    selectInput("DEGdataset",label = "Select the dataset used to perform GSEA:",choices = dataset_list[[selected_option()]],selectize = TRUE)#,selected = "Hallmark")
  })
  
  ranklist <- reactive({
    req(input$rankingmetric, input$filterstatus, input$padjfilter)
    rankgenes(
      data = uploaded_data(),
      rankingmetric = input$rankingmetric,
      filter = input$filterstatus,
      padj_filter = input$padjfilter
    )
  })
  
  
  
  # Adding the input$DEGsubmit in req makes sure that this function only runs when the submit button is pressed. which in turn makes sure that the images are only generated when the button is pressed
  #Changed this from reactive to event reactive to make sure that the images are only plotted when the submit button is pressed
  GSEAoutput <- eventReactive(input$DEGsubmit,{
    req(input$DEGdataset, input$GSEAfilternumeric,input$DEGsubmit)
    GSEA_function(
      dataset_path = file.path("MSIGDB Datasets/msigdb_v2023.1.Hs_Mm_GMTs", input$DEGdataset),
      namedlist = ranklist(),
      padj_filter = input$GSEAfilternumeric,
      title = paste()
    )
  })
  
  output$GSEAplot <- renderPlot({
    req(GSEAoutput())
    GSEAoutput()$filteredplot
  })
  
  output$GSEAtable <- DT::renderDT(server = FALSE,{
    req(GSEAoutput())
    datatable(
      GSEAoutput()$significant_results,
      caption = "Page might crash while loading in large sets of data. Adjust filter in options to prevent crash.",
      extensions = 'Buttons',
      options = list(
        pageLength = 10,
        autoWidth = TRUE, 
        dom = 'Blfrtip', 
        scrollX = TRUE,
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        columnDefs = list(list(width = '400px', targets = ncol(GSEAoutput()$significant_results)))
      )
    )
  })
  
  output$GSEAdotplot <- renderPlot({
    req(GSEAoutput())
    GSEAoutput()$dotplot
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("GSEA_plot_", Sys.Date(), '.png', sep='')
    },
    content = function(file) {
      ggsave(file,plot = GSEAoutput()$filteredplot, width = input$img_width, height = input$img_height, dpi = input$img_res)
    },
    contentType = "image/png"
  )
  
  output$downloadPlot_dp <- downloadHandler(
    filename = function() {
      paste("GSEA_dotplot_",Sys.Date(), '.png', sep='')
    },
    content = function(file) {
      ggsave(file,plot = GSEAoutput()$dotplot, width = input$img_width_dp, height = input$img_height_dp, dpi = input$img_res_dp)
    }
  )
  
## Logic for DE page----
  observeEvent(input$desubmit,{
    mydata <- load_data(path = input$deseqfile$datapath)
    
    output$devolcano <- renderPlotly({
      volcano_plot(dataf = mydata,x_name = input$x_axis,y_name = input$y_axis,slider = input$slider,color1 = input$base, color2 = input$highlight)
    })
    
    output$detable <- renderDataTable({
      draw_table(dataf = mydata,slider = input$slider)
    },width = "100%", options = list(scrollX = TRUE))
  })
}

#Application----
shinyApp(ui = ui, server = server)

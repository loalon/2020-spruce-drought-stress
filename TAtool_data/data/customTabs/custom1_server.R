print("tab 2 server code")
load(file=here("data/metricsData.RData"))  

output$custom1_table <- renderDataTable( {
  
  if (length(rownames(metricsData) > 0 )) {
    metricsData
  } else {
    matrix()
  }
}, selection = 'single')

# 
# observeEvent(input$custom1_table_rows_selected, {
#   req(length(input$custom1_table_rows > 0))
#   
#   selectedGene <- rownames(fullData$networks[[input$network_conditionSelect]]$stats)[input$network_table_rows_selected]
#   
#   selectedFDN <- getGeneFDN(fullData$networks[[input$network_conditionSelect]]$edgeList, selectedGene)
#   #allGenes <<- c(selectedGene, selectedFDN)
#   session$userData$tab_network$currentFDNgenes <- c(selectedGene, selectedFDN)
#   output$network_fdnGenes <- renderText({session$userData$tab_network$currentFDNgenes})
# })
# 
# observeEvent(input$network_toGenxBtn, {
#   if(length(session$userData$tab_network$currentFDNgenes)>0) {
#     updateSelectInput(session, "genX_conditionSelect", selected = input$network_conditionSelect)
#     updateTextAreaInput(session, "genX_genes", value = paste(session$userData$tab_network$currentFDNgenes, collapse=" "))
#     click("genX_loadBtn")
#     updateTabsetPanel(session, "tabs", selected = "tab_genX")
#   }
# })

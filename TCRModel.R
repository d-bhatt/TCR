# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of required packages
required_packages <- c("shiny", "ggplot2", "umap", "dplyr", "scales", "patchwork", "shinycssloaders")

# Install and load them
install_and_load(required_packages)

# ============================= #
# Simulation Functions          #
# ============================= #

set.seed(123)

generate_barcode <- function(n, length) {
  replicate(n, paste0(sample(LETTERS, length, replace=TRUE), collapse=""))
}

check_cross_reactivity <- function(tcell_barcode, antigen_barcode, match_len) {
  L <- nchar(tcell_barcode)
  if (match_len > L) return(FALSE)
  substrings <- sapply(1:(L - match_len + 1), function(i) substr(tcell_barcode, i, i+match_len-1))
  any(substrings %in% sapply(
    1:(nchar(antigen_barcode) - match_len + 1),
    function(j) substr(antigen_barcode, j, j+match_len-1)
  ))
}

init_Tcells <- function(n, grid_size, params) {
  max_tcells <- grid_size * grid_size
  n <- min(n, max_tcells)
  
  random_bc <- generate_barcode(n, params$barcode_length)
  forced_bc <- c(params$self_TCR_barcodes, params$nonself_TCR_barcodes)
  barcodes <- c(random_bc, forced_bc)
  
  Tcells <- data.frame(
    x = runif(length(barcodes), 1, grid_size),
    y = runif(length(barcodes), 1, grid_size),
    state = "resting",
    is_memory = FALSE,
    barcode = barcodes,
    lifespan = sample(80:120, length(barcodes), replace=TRUE),
    activation_level = 0,
    cytotoxic_secretion = 0,
    stimulatory_secretion = 0,
    proliferating = 0
  )
  
  cd4_frac <- runif(1, 0.40, 0.60)
  cd8_frac <- runif(1, 0.20, 0.35)
  
  Tcells$phenotype <- sample(
    c("CD4", "CD8", "OTHER"),
    size = nrow(Tcells),
    replace = TRUE,
    prob = c(cd4_frac, cd8_frac, 1 - cd4_frac - cd8_frac)
  )
  
  if (params$delete_self_reactive) {
    is_reactive <- sapply(Tcells$barcode, function(tc) {
      any(sapply(params$self_TCR_barcodes, function(sb) {
        check_cross_reactivity(tc, sb, params$ordered_match_length)
      }))
    })
    delete_mask <- is_reactive & (runif(nrow(Tcells)) < params$self_delete_prob)
    Tcells <- Tcells[!delete_mask, ]
  }
  
  return(Tcells)
}

init_Antigens <- function(amount, grid_size, params) {
  type_pool <- c(params$self_antigen_barcodes, params$nonself_antigen_barcodes)
  
  if (length(type_pool) < params$n_antigen_types) {
    extra <- params$n_antigen_types - length(type_pool)
    type_pool <- c(type_pool, generate_barcode(extra, params$barcode_length))
  }
  
  antigen_types <- sample(type_pool, amount, replace=TRUE)
  
  data.frame(
    x = runif(amount, 1, grid_size),
    y = runif(amount, 1, grid_size),
    barcode = antigen_types,
    concentration = runif(amount, 0.5, 1)
  )
}

simulate_step <- function(Tcells, Antigens, grid_size, params, step) {
  Tcells$x <- pmin(pmax(Tcells$x + rnorm(nrow(Tcells), 0, 0.8), 1), grid_size)
  Tcells$y <- pmin(pmax(Tcells$y + rnorm(nrow(Tcells), 0, 0.8), 1), grid_size)
  
  if (step %in% params$antigen_introduce_steps)
    Antigens <- rbind(Antigens, init_Antigens(params$initial_antigen_amount, grid_size, params))
  
  if (step %in% params$antigen_reintroduce_steps)
    Antigens <- rbind(Antigens, init_Antigens(params$reintroduction_antigen_amount, grid_size, params))
  
  Antigens$concentration <- pmax(0, Antigens$concentration - params$antigen_decay)
  
  radius <- 2
  
  for (i in seq_len(nrow(Tcells))) {
    tcell <- Tcells[i, ]
    d <- sqrt((Antigens$x - tcell$x)^2 + (Antigens$y - tcell$y)^2)
    close_ag <- which(d < radius)
    
    if (length(close_ag) > 0) {
      matched <- FALSE
      exact_match <- FALSE
      matched_idx <- NA
      
      for (ag_idx in close_ag) {
        ag_bc <- Antigens$barcode[ag_idx]
        if (ag_bc %in% params$self_antigen_barcodes) next
        
        if (tcell$barcode == ag_bc) {
          matched <- TRUE
          exact_match <- TRUE
          matched_idx <- ag_idx
          break
        }
        
        if (check_cross_reactivity(tcell$barcode, ag_bc, params$ordered_match_length)) {
          matched <- TRUE
          matched_idx <- ag_idx
          break
        }
      }
      
      if (matched) {
        Tcells$state[i] <- "active"
        activation_increment <- ifelse(exact_match, 0.12, 0.06)
        Tcells$activation_level[i] <- Tcells$activation_level[i] + activation_increment
        Antigens$concentration[matched_idx] <- max(0, Antigens$concentration[matched_idx] - 0.1)
        
        if (params$study_functional_responses) {
          if (tcell$phenotype == "CD4") {
            Tcells$stimulatory_secretion[i] <- runif(1, 0.4*activation_increment, 1.0*activation_increment)
            Tcells$cytotoxic_secretion[i] <- 0
          }
          if (tcell$phenotype == "CD8") {
            Tcells$cytotoxic_secretion[i] <- runif(1, 0.3*activation_increment, 1.2*activation_increment)
            Tcells$stimulatory_secretion[i] <- runif(1, 0.1*activation_increment, 0.5*activation_increment)
          }
        }
        
        if (runif(1) < params$Tcell_prolif_active) {
          new_cell <- tcell
          new_cell$x <- min(max(tcell$x + rnorm(1, 0, 1), 1), grid_size)
          new_cell$y <- min(max(tcell$y + rnorm(1, 0, 1), 1), grid_size)
          new_cell$is_memory <- runif(1) < params$memory_prob
          Tcells <- rbind(Tcells, new_cell)
        }
      }
    }
    
    Tcells$lifespan[i] <- Tcells$lifespan[i] - 1
  }
  
  Tcells <- Tcells[Tcells$lifespan > 0 | Tcells$is_memory, ]
  Antigens <- Antigens[Antigens$concentration > 0.05, ]
  
  list(Tcells = Tcells, Antigens = Antigens)
}

run_simulation <- function(params) {
  Tcells <- init_Tcells(params$n_Tcells, params$grid_size, params)
  Antigens <- init_Antigens(params$initial_antigen_amount, params$grid_size, params)
  history <- vector("list", params$n_steps)
  
  for (t in seq_len(params$n_steps)) {
    res <- simulate_step(Tcells, Antigens, params$grid_size, params, t)
    Tcells <- res$Tcells
    Antigens <- res$Antigens
    history[[t]] <- res
  }
  
  history

}

# ============================================================
# UMAP Feature Extraction and Helpers
# ============================================================

extract_tcell_features <- function(Tcells) {
  required <- c(
    "activation_level", "cytotoxic_secretion", "stimulatory_secretion",
    "state", "is_memory", "proliferating", "phenotype"
  )
  
  missing_cols <- setdiff(required, colnames(Tcells))
  if (length(missing_cols) > 0)
    stop("Missing columns in Tcells: ", paste(missing_cols, collapse = ", "))
  
  features <- Tcells %>%
    mutate(
      state_num = as.numeric(factor(state, levels = c("resting", "active", "proliferating"))),
      memory_num = as.numeric(is_memory),
      proliferating_num = as.numeric(proliferating),
      phenotype_num = as.numeric(factor(phenotype, levels = c("CD4", "CD8", "OTHER")))
    ) %>%
    select(
      activation_level,
      cytotoxic_secretion,
      stimulatory_secretion,
      state_num,
      memory_num,
      proliferating_num,
      phenotype_num
    ) %>%
    na.omit()
  
  return(features)
}

run_umap <- function(features, n_neighbors = 50, min_dist = 0.8, metric = "manhattan") {
  n_neighbors <- min(nrow(features) - 1, n_neighbors)
  layout <- umap(features, n_neighbors = n_neighbors, min_dist = min_dist)$layout
  colnames(layout) <- c("UMAP1", "UMAP2")
  return(layout)
}

merge_umap_metadata <- function(layout, Tcells) {
  data.frame(layout, 
             barcode = Tcells$barcode,
             phenotype = Tcells$phenotype,
             state = Tcells$state,
             is_memory = Tcells$is_memory,
             proliferating = Tcells$proliferating,
             cytotoxic_secretion = Tcells$cytotoxic_secretion,
             stimulatory_secretion = Tcells$stimulatory_secretion)
}


# ============================= #
# Shiny App UI and Server       #
# ============================= #

# UI
ui <- fluidPage(
  titlePanel("T Cell–Repertoire Dynamics in Response to Antigen Stimuli"),
  tags$h5("Created by Darshak K. Bhatt, 2025"),
  sidebarLayout(
    sidebarPanel(
      width = 3, 
      tabsetPanel(
        tabPanel("Simulation",
                 numericInput("n_Tcells", "Number of T cells", 500),
                 numericInput("n_steps", "Simulation steps", 200),
                 numericInput("grid_size", "Grid size", 50),
                 numericInput("n_replicates", "Replicates", 1),
                 textInput("target_barcode", "Target barcode", "DKB"),
                 numericInput("plot_interval", "Plot interval", 10),
                 checkboxInput("do_plot", "Enable plotting", FALSE)
        ),
        tabPanel("Barcodes",
                 numericInput("barcode_length", "Barcode length", 3),
                 numericInput("ordered_match_length", "Ordered match length", 2),
                 textInput("self_TCR_barcodes", "Self TCR barcodes (comma-separated)", "PQR"),
                 textInput("nonself_TCR_barcodes", "Nonself TCR barcodes", "DKB,STR"),
                 checkboxInput("delete_self_reactive", "Delete self-reactive", TRUE),
                 numericInput("self_delete_prob", "Self delete probability", 1.0)
        ),
        tabPanel("Antigens",
                 textInput("self_antigen_barcodes", "Self antigen barcodes", "AAA,ATR"),
                 textInput("nonself_antigen_barcodes", "Nonself antigen barcodes", "UVW,DKB,STR"),
                 numericInput("n_antigen_types", "Number of antigen types", 20),
                 numericInput("initial_antigen_amount", "Initial antigen amount", 500),
                 numericInput("reintroduction_antigen_amount", "Reintroduction amount", 500),
                 textInput("antigen_introduce_steps", "Introduce steps (comma-separated)", "1"),
                 textInput("antigen_reintroduce_steps", "Reintroduce steps", "30,60"),
                 numericInput("antigen_decay", "Antigen decay rate", 0.05)
        ),
        tabPanel("T Cell Dynamics",
                 numericInput("Tcell_prolif_resting", "Prolif. resting", 0.05),
                 numericInput("Tcell_prolif_active", "Prolif. active", 0.3),
                 numericInput("memory_prob", "Memory probability", 0.1)
        ),
        tabPanel("Functional Output",
                 numericInput("min_cytotoxic_secretion", "Min cytotoxic", 0),
                 numericInput("max_cytotoxic_secretion", "Max cytotoxic", 1),
                 numericInput("min_stimulatory_secretion", "Min stimulatory", 0),
                 numericInput("max_stimulatory_secretion", "Max stimulatory", 1),
                 checkboxInput("study_functional_responses", "Study functional responses (Check for UMAP)", FALSE)
        )
      ),
      actionButton("runSim", "Run Simulation")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Realtime Viewer",
                 sliderInput("stepSlider", "Step", min = 1, max = 200, value = 1, step = 1, animate = TRUE),
                 plotOutput("tcellPlot", height = "600px")
        ),
        tabPanel("Clonality & Diversity",
                 sliderInput("fig1Step", "Dashboard Timestep", min = 1, max = 200, value = 80),
                 withSpinner(plotOutput("fig1Plot", height = "800px"), type = 6, color = "#0072B2")
        ),
        tabPanel("UMAP (Only when studying Effector functions)",
                 sliderInput("fig2Step", "UMAP Timestep", min = 1, max = 200, value = 80),
                 withSpinner(plotOutput("fig2Plot", height = "800px"), type = 6, color = "#0072B2")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  sim_data <- reactiveVal(NULL)
  
  observeEvent(input$runSim, {
    params <- list(
      grid_size = input$grid_size,
      n_Tcells = input$n_Tcells,
      n_steps = input$n_steps,
      n_replicates = input$n_replicates,
      target_barcode = input$target_barcode,
      plot_interval = input$plot_interval,
      do_plot = input$do_plot,
      
      barcode_length = input$barcode_length,
      ordered_match_length = input$ordered_match_length,
      self_TCR_barcodes = strsplit(input$self_TCR_barcodes, ",")[[1]],
      nonself_TCR_barcodes = strsplit(input$nonself_TCR_barcodes, ",")[[1]],
      delete_self_reactive = input$delete_self_reactive,
      self_delete_prob = input$self_delete_prob,
      
      self_antigen_barcodes = strsplit(input$self_antigen_barcodes, ",")[[1]],
      nonself_antigen_barcodes = strsplit(input$nonself_antigen_barcodes, ",")[[1]],
      n_antigen_types = input$n_antigen_types,
      initial_antigen_amount = input$initial_antigen_amount,
      reintroduction_antigen_amount = input$reintroduction_antigen_amount,
      antigen_introduce_steps = as.numeric(unlist(strsplit(input$antigen_introduce_steps, ","))),
      antigen_reintroduce_steps = as.numeric(unlist(strsplit(input$antigen_reintroduce_steps, ","))),
      antigen_decay = input$antigen_decay,
      
      Tcell_prolif_resting = input$Tcell_prolif_resting,
      Tcell_prolif_active = input$Tcell_prolif_active,
      memory_prob = input$memory_prob,
      
      min_cytotoxic_secretion = input$min_cytotoxic_secretion,
      max_cytotoxic_secretion = input$max_cytotoxic_secretion,
      min_stimulatory_secretion = input$min_stimulatory_secretion,
      max_stimulatory_secretion = input$max_stimulatory_secretion,
      study_functional_responses = input$study_functional_responses
    )
    
    history <- run_simulation(params)
    sim_data(history)
    updateSliderInput(session, "stepSlider", max = input$n_steps, value = 1)
  })
  
  output$tcellPlot <- renderPlot({
    req(sim_data())
    step <- input$stepSlider
    Tcells <- sim_data()[[step]]$Tcells
    Antigens <- sim_data()[[step]]$Antigens
    
    ggplot() +
      # Antigen spots (background layer)
      geom_point(
        data = Antigens,
        aes(x = x, y = y),
        color = "blue", shape = 4, size = 2, alpha = 0.8
      ) +
      
      # T cell layer
      geom_point(
        data = Tcells,
        aes(x = x, y = y, color = state, shape = phenotype),
        alpha = 0.7, size = 4
      ) +
      
      # Manual color and shape scales
      scale_color_manual(
        values = c("resting" = "gray80", "active" = "red")
      ) +
      scale_shape_manual(
        values = c("CD4" = 16, "CD8" = 17, "OTHER" = 15)
      ) +
      
      # Fixed coordinate system
      coord_fixed(xlim = c(0, 50), ylim = c(0, 50)) +
      
      # Labels
      labs(
        title = paste("T Cell and Antigen Positions at Step", step),
        x = "X Position", y = "Y Position",
        color = "T Cell State", shape = "Phenotype"
      ) +
      
      # Theme with responsive fonts
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  })
  


  output$fig1Plot <- renderPlot({
    req(sim_data())
    history <- sim_data()
    t <- input$fig1Step
    if (t > length(history)) t <- length(history)
    
    Tcells <- history[[t]]$Tcells
    params <- list(
      self_antigen_barcodes = strsplit(input$self_antigen_barcodes, ",")[[1]],
      nonself_antigen_barcodes = strsplit(input$nonself_antigen_barcodes, ",")[[1]],
      self_TCR_barcodes = strsplit(input$self_TCR_barcodes, ",")[[1]],
      nonself_TCR_barcodes = strsplit(input$nonself_TCR_barcodes, ",")[[1]],
      delete_self_reactive = input$delete_self_reactive,
      ordered_match_length = input$ordered_match_length,
      target_barcode = input$target_barcode,
      n_steps = input$n_steps,
      grid_size = input$grid_size
    )
    
    clone_counts <- sort(table(Tcells$barcode), decreasing = TRUE)
    top_clones <- names(clone_counts)[1:min(5, length(clone_counts))]
    
    # matrix of clone trajectories
    clone_traj <- sapply(history, function(h) {
      tab <- table(h$Tcells$barcode)
      sapply(top_clones, function(cl) ifelse(cl %in% names(tab), tab[cl], 0))
    })
    clone_traj <- t(clone_traj)  # rows = steps, cols = clones
    
    antigen_load <- sapply(history, function(h) sum(h$Antigens$concentration))
    diversity <- sapply(history, function(h) {
      tab <- table(h$Tcells$barcode)
      p <- tab / sum(tab)
      -sum(p * log(p + 1e-9))
    })
    
    # Convert to tidy data frames
    df_clone <- as.data.frame(clone_traj)
    df_clone$step <- seq_len(nrow(df_clone))
    df_clone <- tidyr::pivot_longer(df_clone, -step,
                                    names_to = "clone", values_to = "size")
    
    df_antigen <- data.frame(step = seq_along(history), load = antigen_load)
    df_diversity <- data.frame(step = seq_along(history), H = diversity)
    df_bar <- data.frame(clone = names(clone_counts),
                         count = as.numeric(clone_counts))
    
    # 1. Top Clone Dynamics
    p1 <- ggplot(df_clone, aes(x = step, y = size, color = clone)) +
      geom_line(linewidth = 1.2) +
      geom_vline(xintercept = t, linetype = "dashed", color = "grey40") +
      labs(title = "Top Clone Dynamics", x = "Time step", y = "Clone size") +
      theme_minimal(base_size = 14)
    
    # 2. Antigen Load
    p2 <- ggplot(df_antigen, aes(x = step, y = load)) +
      geom_line(color = "orange", linewidth = 1.2) +
      geom_vline(xintercept = t, linetype = "dashed", color = "grey40") +
      labs(title = "Antigen Load Over Time", x = "Time step", y = "Total antigen concentration") +
      theme_minimal(base_size = 14)
    
    # 3. Clonal Diversity
    p3 <- ggplot(df_diversity, aes(x = step, y = H)) +
      geom_line(color = "purple", linewidth = 1.2) +
      geom_vline(xintercept = t, linetype = "dashed", color = "grey40") +
      labs(title = "Clonal Diversity (H')", x = "Time step", y = "H'") +
      theme_minimal(base_size = 14)
    
    # 4. Barplot of Top Clones
    p4 <- ggplot(head(df_bar, 10), aes(x = reorder(clone, -count), y = count)) +
      geom_col(fill = "skyblue") +
      labs(title = paste("Top Clones at Step", t), x = "Clone", y = "Count") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 5. Pie Chart
    df_pie <- head(df_bar, min(8, nrow(df_bar)))
    p5 <- ggplot(df_pie, aes(x = "", y = count, fill = clone)) +
      geom_col(width = 1, color = "black") +
      # Add labels inside slices
      geom_text(aes(label = paste0(clone, "\n", count)),
                position = position_stack(vjust = 0.5),
                color = "white", size = 4, fontface = "bold") +
      coord_polar(theta = "y") +
      labs(title = paste("Clonal Composition — Step", t)) +
      theme_void(base_size = 14) +
      theme(legend.position = "none")
    
    # 6. Text Summary
    summary_text <- paste0(
      "Simulation setup (Step ", t, "):\n",
      "- Self antigens: ", paste(params$self_antigen_barcodes, collapse = ", "), "\n",
      "- Non-self antigens: ", paste(params$nonself_antigen_barcodes, collapse = ", "), "\n",
      "- Self-reactive TCRs: ", paste(params$self_TCR_barcodes, collapse = ", "), "\n",
      "- Target/non-self TCRs: ", paste(params$nonself_TCR_barcodes, collapse = ", "), "\n",
      "- Self-reactive TCR deletion: ", params$delete_self_reactive, "\n",
      "- Cross-reactivity threshold: ", params$ordered_match_length
    )
    
    p6 <- ggplot() +
      annotate("text", x = 0, y = 1, label = summary_text,
               hjust = 0, vjust = 1, size = 6) +
      xlim(0, 1) + ylim(0, 1) +
      theme_void()
    
    # Combine plots with patchwork
    (p1 | p2 | p3) /
      (p4 | p5 | p6)
  })
  
  
  
  #Figure 2 code
  output$fig2Plot <- renderPlot({
    req(sim_data())
    if (!input$study_functional_responses) {
      plot.new()
      title("UMAP only available when functional responses are studied (Check Functional Output tab)")
      return()
    }
    
    history <- sim_data()
    t <- input$fig2Step
    if (t > length(history)) t <- length(history)
    Tcells <- history[[t]]$Tcells
    if (nrow(Tcells) == 0) {
      plot.new()
      title("No T cells present at this timestep")
      return()
    }
    
    # Feature extraction
    feats <- tryCatch(extract_tcell_features(Tcells), error = function(e) NULL)
    if (is.null(feats) || nrow(feats) < 2) {
      plot.new()
      title("Insufficient data for UMAP")
      return()
    }
    
    layout <- run_umap(feats)
    umap_df <- merge_umap_metadata(layout, Tcells)
    
    # Top clone coloring
    clone_counts <- sort(table(umap_df$barcode), decreasing = TRUE)
    top10 <- names(clone_counts)[1:min(10, length(clone_counts))]
    umap_df$clone_group <- ifelse(umap_df$barcode %in% top10, umap_df$barcode, "Other")
    umap_df$clone_group <- factor(umap_df$clone_group, levels = c(top10, "Other"))
    clone_colors <- c(setNames(rainbow(length(top10)), top10), Other = "gray70")
    
    # K-means clustering
    set.seed(42)
    k <- min(10, nrow(umap_df))
    clusters <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = k)
    umap_df$Population <- factor(clusters$cluster)
    centroids <- umap_df %>%
      group_by(Population) %>%
      summarise(across(c(UMAP1, UMAP2), mean))
    
    # Define phenotype colors
    phenotype_colors <- c("CD4" = "blue", "CD8" = "red", "OTHER" = "gray40")
    
    # Create ggplot panels
    p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = phenotype)) +
      geom_point(alpha = 0.7, size = 1) +
      scale_color_manual(values = phenotype_colors) +
      theme_minimal(base_size = 14) + ggtitle("Phenotype")+
      coord_fixed()
    
    p2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = state)) +
      geom_point(alpha = 0.6, size = 1) +
      theme_minimal(base_size = 14) + ggtitle("Activation State")+
      coord_fixed()
    
    p3 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = is_memory)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("TRUE" = "purple", "FALSE" = "gray70")) +
      theme_minimal(base_size = 14) + ggtitle("Memory Status")+
      coord_fixed()
      
    
    p4 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = stimulatory_secretion)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis_c() +
      theme_minimal(base_size = 14) + ggtitle("Stimulatory Secretion")+
      coord_fixed()
      
    
    p5 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cytotoxic_secretion)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis_c() +
      theme_minimal(base_size = 14) + ggtitle("Cytotoxic Secretion")+
      coord_fixed()
      
    
    p6 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Population)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_text(data = centroids, aes(label = Population),
                color = "black", size = 4, fontface = "bold") +
      theme_minimal(base_size = 14) + ggtitle("K-means Clusters")+
      coord_fixed()
      
    
    # Combine with patchwork
    (p1 | p2 | p3) / (p4 | p5 | p6)
  })
  
  
  
}

# Launch app
shinyApp(ui = ui, server = server)
# render plots with Tikz
render_plots <- function(plot, save_path, file_name,
                         plot_parameters = list(
                             stand_alone = FALSE,
                             compile_pdf = TRUE,
                             compile_png = FALSE,
                             aspect_ratio = 6 / 19,
                             width = 7
                         ), ...) {
    # Variables
    name.p <- as_label(dplyr::enquo(plot))
    # final_plot_name
    plot_name <- glue::glue("{name.p}_{file_name}")
    stand_alone <- plot_parameters$stand_alone
    compile_pdf <- plot_parameters$compile_pdf
    compile_pdf <- plot_parameters$compile_pdf
    compile_png <- plot_parameters$compile_png
    aspect_ratio <- plot_parameters$aspect_ratio
    width <- plot_parameters$width

    # chekc if destination exists
    dir.create(file.path(here::here(save_path)), recursive = TRUE, showWarnings = FALSE)

    # tikzpicture
    write_tikzpicture(
        plot, plot_name, save_path,
        width = width, aspect_ratio = aspect_ratio,
        stand_alone = stand_alone)

    # render tex
    if(stand_alone == TRUE && compile_pdf == TRUE){
        tinytex::latexmk(here::here(save_path, paste0(plot_name,".tex")))
    }
    # render png
    if(stand_alone == TRUE && compile_pdf == TRUE && compile_png == TRUE){
        pdf_to_png(here::here(save_path, paste0(plot_name,".pdf")))
    }
}


pdf_to_png <- function(pdf_path, dpi = 600, density = dpi, page = 1, remove_pdf = FALSE) {
    # Output file path
    png_path <- sub("\\.pdf$", ".png", pdf_path)

    # Read the PDF (you can specify density for resolution)
    img <- image_read_pdf(pdf_path, density = density, pages = page)

    # Write to PNG
    image_write(img, path = png_path, format = "png")

    if (remove_pdf) file.remove(pdf_path)

    invisible(png_path)
}


# plotting with tikzpicutre
write_tikzpicture <- function(
        plot_obj, plot_name, save_path, width = 7,
        aspect_ratio = 9 / 16, stand_alone = FALSE, compile_pdf = TRUE, ...) {
    tikz(here::here(save_path, glue::glue("{plot_name}.tex")),
         standAlone = stand_alone, width = width, height = width * aspect_ratio)
    print(plot_obj)
    dev.off()
}

own_theme <- function(...) {
    own_theme <- ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = element_text(face = "bold",
                                      hjust = 0.5,
                                      margin = ggplot2::margin(b = 10)),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(margin = margin(t = 5)),  # x specific
            axis.text.y = element_text(margin = margin(r = 5)),  # y specific
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            panel.spacing = unit(0.5, "cm"),
            #   axis.line = element_line(colour = '#004c93'),
            panel.border = element_rect(linewidth = 0.2), # colour = '#004c93', linewidth = 0.01),
            panel.grid.major.x = element_blank(),
            #      axis.ticks = element_line(colour = '#004c93'),
            panel.grid.minor.x = element_blank(),
            strip.background = element_rect(colour = '#004c93',
                                            fill = '#004c93'),
            strip.text.x = element_text(color = 'white', size = 12,
                                        margin = margin(t = 5, r = 5, b = 5, l = 5)),
            strip.text.y = element_text(color = 'white', size = 12,
                                        margin = margin(t = 5, r = 5, b = 5, l = 5))
            )

    return(own_theme)
}

ci_plot <- function(ci_data, plot_name, standAlone = FALSE,
                    plot_width = 10, plot_height = 12, ...){
    
    
    
    # create export folderif it does not exist
    export_folder <- here::here("04_export")
    dir.create(file.path(export_folder), recursive = TRUE, showWarnings = FALSE)
    
    # empty plot
    spacer <- ggplot2::ggplot() +
        ggplot2::theme_void()
    
    # Coverage plot
    cov_plot <- ci_data %>%
        ggplot2::ggplot(aes(x = CI_method, y = Coverage, color = CI_method)) +
        geom_rect(
            data = data.frame(
                xmin = c(2.5, 3.5),
                xmax = c(3.5, 4.5),
                ymin = -Inf,
                ymax = Inf
            ),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey75",
            alpha = 0.4
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::facet_grid(noise ~ split.R, scales = "free_y") +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::labs(x = NULL, y = NULL, color = NULL, title = "Coverage") +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_manual(values = c("$\\text{BLB}^{V}$" = "#004c93", "$\\text{BLB}$" = "#004c93",
                                               "$\\text{IJ}^{V}$" = "#000000", "$\\text{IJ}$" = "#000000")) +
        own_theme()
    
    # CI width
    width_plot <- ci_data  %>%
        ggplot2::ggplot(aes(x = CI_method, y = `95\\% CI Width`, color = CI_method)) +
        geom_rect(
            data = data.frame(
                xmin = c(2.5, 3.5),
                xmax = c(3.5, 4.5),
                ymin = -Inf,
                ymax = Inf
            ),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey75",
            alpha = 0.4
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::facet_grid(noise ~ split.R, scales = "free_y") +
        ggplot2::labs(x = NULL, y = NULL, color = NULL,
                      title = "95\\% CI Width") +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_manual(values = c("$\\text{BLB}^{V}$" = "#004c93", "$\\text{BLB}$" = "#004c93",
                                               "$\\text{IJ}^{V}$" = "#000000", "$\\text{IJ}$" = "#000000")) +
        own_theme()
    
    # Grob
    g1 <- ggplotGrob(cov_plot)
    g2 <- ggplotGrob(width_plot)
    
    # Force equal widths (this is what aligns the plot panels)
    g1$widths <- g2$widths <- grid::unit.pmax(g1$widths, g2$widths)
    
    # combine plot
    combined_plot <- cowplot::plot_grid(
        ggpubr::as_ggplot(g1),
        spacer,
        ggpubr::as_ggplot(g2),
        ncol = 1,
        rel_heights = c(1, .05, 1)
    )
    
    # write tikzpicture
    write_tikzpicture(combined_plot, plot_name = plot_name, 
                      save_path = export_folder,
                      standAlone = standAlone, 
                      width = plot_width, height = plot_height)
    
}


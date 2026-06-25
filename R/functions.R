#' FeaturePlot with cluster boundaries
#'
#' Draw a Seurat FeaturePlot and overlay cluster boundaries computed from the
#' embedding using DBSCAN + concave hulls, with optional overlap resolution
#' and cluster labels.
#'
#' @param seurat_obj A Seurat object.
#' @param feature Feature to plot in `Seurat::FeaturePlot()`.
#' @param reduction Reduction name passed to `Seurat::Embeddings()` and `FeaturePlot()`.
#' @param group.by Metadata column used to define clusters.
#' @param dims Dimensions of the reduction to use.
#' @param eps DBSCAN epsilon; if `NULL`, estimated automatically from kNN distances.
#' @param minPts DBSCAN `minPts`.
#' @param min_cells_component Minimum number of cells required to keep one component.
#' @param concavity Concavity parameter for `concaveman::concaveman()`.
#' @param length_threshold Length threshold for `concaveman::concaveman()`.
#' @param trim_q Quantile used to trim sparse peripheral points before hull creation.
#' @param seed Random seed used in overlap resolution and optional label repulsion.
#' @param pt.size Point size for `Seurat::FeaturePlot()`.
#' @param border.size Boundary line width.
#' @param border.color Boundary color.
#' @param label.size Label size.
#' @param label.color Label color.
#' @param order Passed to `Seurat::FeaturePlot()`.
#' @param raster Passed to `Seurat::FeaturePlot()`.
#' @param repel Logical; use `ggrepel::geom_text_repel()` for labels.
#' @param palette Colors passed to `ggplot2::scale_color_gradientn()`.
#' @param ... Additional arguments passed to `Seurat::FeaturePlot()`.
#'
#' @return A ggplot object.
#' @export
plot_feature_with_cluster_boundaries <- function(
    seurat_obj,
    feature,
    reduction = "umap",
    group.by,
    dims = c(1, 2),
    eps = 0.35,
    minPts = 20,
    min_cells_component = 30,
    concavity = 1.6,
    length_threshold = 0,
    trim_q = 0.995,
    seed = 123,
    pt.size = 0.25,
    border.size = 0.7,
    border.color = "gray40",
    label.size = 5,
    label.color = "black",
    order = TRUE,
    raster = FALSE,
    repel = FALSE,
    palette = c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026"),
    ...
) {
  needed <- c("Seurat", "ggplot2", "dbscan", "concaveman", "FNN", "sf", "polyclip")
  if (isTRUE(repel)) {
    needed <- c(needed, "ggrepel")
  }

  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop(
      "Missing suggested packages for `plot_feature_with_cluster_boundaries()`: ",
      paste(needed[!ok], collapse = ", "),
      call. = FALSE
    )
  }

  polygon_area <- function(x, y) {
    n <- length(x)
    if (n < 3) return(0)
    abs(sum(x[-1] * y[-n] - x[-n] * y[-1]) / 2)
  }

  make_boundaries <- function(df) {
    if (nrow(df) < 3) {
      stop("Not enough cells with complete embedding/metadata values.", call. = FALSE)
    }

    eps_use <- eps
    if (is.null(eps_use)) {
      if (nrow(df) < 2) {
        stop("Cannot auto-estimate `eps` with fewer than 2 cells.", call. = FALSE)
      }
      k_auto <- min(10, nrow(df) - 1)
      nn_all <- FNN::get.knn(as.matrix(df[, c("x", "y")]), k = k_auto)
      eps_use <- stats::median(nn_all$nn.dist[, k_auto], na.rm = TRUE) * 1.5
    }

    out <- list()
    j <- 1L

    for (cl in sort(unique(df$cluster))) {
      sub <- df[df$cluster == cl, c("cell", "x", "y", "cluster"), drop = FALSE]
      if (nrow(sub) < max(3, minPts)) next

      db <- dbscan::dbscan(as.matrix(sub[, c("x", "y")]), eps = eps_use, minPts = minPts)
      sub$component <- db$cluster

      comps <- sort(unique(sub$component[sub$component > 0]))
      if (length(comps) == 0) {
        sub$component <- 1L
        comps <- 1L
      }

      for (comp_id in comps) {
        comp <- sub[sub$component == comp_id, c("x", "y", "cluster"), drop = FALSE]
        if (nrow(comp) < min_cells_component) next

        if (nrow(comp) > 15) {
          k_trim <- min(10, nrow(comp) - 1)
          nn <- FNN::get.knn(as.matrix(comp[, c("x", "y")]), k = k_trim)
          local_scale <- rowMeans(nn$nn.dist)
          keep <- local_scale <= stats::quantile(local_scale, probs = trim_q, na.rm = TRUE)
          if (sum(keep) >= 3) {
            comp <- comp[keep, , drop = FALSE]
          }
        }

        pts <- as.matrix(comp[, c("x", "y"), drop = FALSE])

        hull <- tryCatch(
          concaveman::concaveman(
            pts,
            concavity = concavity,
            length_threshold = length_threshold
          ),
          error = function(e) NULL
        )

        if (is.null(hull) || nrow(hull) < 3) {
          idx_h <- grDevices::chull(pts[, 1], pts[, 2])
          hull <- pts[c(idx_h, idx_h[1]), , drop = FALSE]
        }

        hull <- as.data.frame(hull)
        colnames(hull) <- c("x", "y")
        hull$cluster <- cl
        hull$component <- paste(cl, comp_id, sep = "__")
        hull$order <- seq_len(nrow(hull))

        out[[j]] <- hull
        j <- j + 1L
      }
    }

    if (length(out) == 0) {
      stop("No boundaries generated. Try lowering `min_cells_component` or `minPts`.", call. = FALSE)
    }

    do.call(rbind, out)
  }

  resolve_overlaps <- function(boundaries) {
    split_polys <- split(
      boundaries[, c("x", "y", "cluster", "component", "order")],
      boundaries$component
    )

    poly_list <- lapply(split_polys, function(dd) {
      dd <- dd[base::order(dd$order), , drop = FALSE]
      first_xy <- as.numeric(dd[1, c("x", "y")])
      last_xy  <- as.numeric(dd[nrow(dd), c("x", "y")])

      if (!isTRUE(all.equal(first_xy, last_xy))) {
        dd <- rbind(dd, dd[1, , drop = FALSE])
      }

      list(
        x = dd$x,
        y = dd$y,
        cluster = dd$cluster[1],
        component = dd$component[1]
      )
    })

    set.seed(seed)

    ord_idx <- base::order(
      vapply(poly_list, function(z) polygon_area(z$x, z$y), numeric(1)),
      decreasing = TRUE
    )

    kept <- list()
    meta <- list()

    as_polyclip <- function(p) list(list(x = p$x, y = p$y))

    for (ii in ord_idx) {
      cur <- poly_list[[ii]]
      cur_pc <- as_polyclip(cur)

      if (length(kept) > 0) {
        for (kk in seq_along(kept)) {
          prev_pc <- as_polyclip(kept[[kk]])
          cur_pc <- polyclip::polyclip(cur_pc, prev_pc, op = "minus")
          if (length(cur_pc) == 0) break
        }
      }

      if (length(cur_pc) == 0) next

      for (part in seq_along(cur_pc)) {
        kept[[length(kept) + 1L]] <- list(
          x = cur_pc[[part]]$x,
          y = cur_pc[[part]]$y
        )
        meta[[length(meta) + 1L]] <- data.frame(
          cluster = cur$cluster,
          component = cur$component,
          part = part,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(kept) == 0) {
      stop("All polygons were removed during overlap resolution.", call. = FALSE)
    }

    out <- vector("list", length(kept))
    for (i in seq_along(kept)) {
      out[[i]] <- data.frame(
        x = kept[[i]]$x,
        y = kept[[i]]$y,
        cluster = meta[[i]]$cluster,
        component = paste0(meta[[i]]$component, "__part", meta[[i]]$part),
        order = seq_along(kept[[i]]$x),
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, out)
  }

  get_label_positions <- function(boundaries_df) {
    split_polys <- split(
      boundaries_df[, c("x", "y", "cluster", "component", "order")],
      boundaries_df$component
    )

    poly_sf <- lapply(split_polys, function(dd) {
      dd <- dd[base::order(dd$order), , drop = FALSE]
      xy <- as.matrix(dd[, c("x", "y"), drop = FALSE])

      if (!isTRUE(all.equal(xy[1, ], xy[nrow(xy), ]))) {
        xy <- rbind(xy, xy[1, , drop = FALSE])
      }

      sf::st_sf(
        cluster = dd$cluster[1],
        component = dd$component[1],
        geometry = sf::st_sfc(sf::st_polygon(list(xy)))
      )
    })

    poly_sf <- do.call(rbind, poly_sf)
    poly_sf$area <- as.numeric(sf::st_area(poly_sf))

    keep <- unlist(tapply(seq_len(nrow(poly_sf)), poly_sf$cluster, function(ii) {
      ii[which.max(poly_sf$area[ii])]
    }))

    poly_sf <- poly_sf[keep, , drop = FALSE]
    pts <- sf::st_point_on_surface(poly_sf)
    coords <- sf::st_coordinates(pts)

    data.frame(
      x = coords[, "X"],
      y = coords[, "Y"],
      cluster = pts$cluster,
      stringsAsFactors = FALSE
    )
  }

  emb <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = reduction)[, dims, drop = FALSE])
  colnames(emb) <- c("x", "y")
  emb$cell <- rownames(emb)

  meta <- seurat_obj@meta.data
  meta$cell <- rownames(meta)

  if (!group.by %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in seurat_obj@meta.data", group.by), call. = FALSE)
  }

  idx <- match(emb$cell, meta$cell)
  plot_df <- data.frame(
    cell = emb$cell,
    x = emb$x,
    y = emb$y,
    cluster = as.character(meta[idx, group.by, drop = TRUE]),
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[stats::complete.cases(plot_df), , drop = FALSE]
  if (nrow(plot_df) == 0) {
    stop("No cells left after removing incomplete cases.", call. = FALSE)
  }

  boundaries <- make_boundaries(plot_df)
  boundaries2 <- resolve_overlaps(boundaries)
  label_df <- get_label_positions(boundaries2)

  p <- Seurat::FeaturePlot(
    object = seurat_obj,
    features = feature,
    reduction = reduction,
    dims = dims,
    pt.size = pt.size,
    order = order,
    raster = raster,
    combine = TRUE,
    ...
  ) +
    ggplot2::scale_color_gradientn(colors = palette) +
    ggplot2::geom_path(
      data = boundaries2,
      ggplot2::aes(x = x, y = y, group = component),
      inherit.aes = FALSE,
      color = border.color,
      linewidth = border.size,
      lineend = "round",
      linejoin = "round"
    )

  if (isTRUE(repel)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        color = label.color,
        size = label.size,
        fontface = "bold",
        seed = seed
      )
  } else {
    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        color = label.color,
        size = label.size,
        fontface = "bold"
      )
  }

  p
}

#' Stacked ribbon plot
#'
#' Create a stacked bar plot connected by smooth ribbons between adjacent x-axis
#' categories, with optional faceting and automatic shortening/wrapping of axis
#' labels.
#'
#' @param data A data frame.
#' @param x String. Column name mapped to the x-axis category.
#' @param y String. Column name containing numeric values to stack.
#' @param fill String. Column name used for stacked fill groups.
#' @param facet_rows Optional character vector of column names used for facet rows.
#' @param facet_cols Optional character vector of column names used for facet columns.
#' @param x_levels Optional character vector giving the order of x categories.
#' @param fill_levels Optional character vector giving the order of fill groups.
#' @param bar_width Width of each stacked bar.
#' @param ribbon_alpha Alpha transparency for ribbons.
#' @param ribbon_colour Outline colour for ribbons.
#' @param bar_colour Outline colour for bars.
#' @param bar_linewidth Outline width for bars.
#' @param ribbon_n Number of interpolation points used for each ribbon polygon.
#' @param reverse_y Logical; if `TRUE`, reverse the y-axis.
#' @param facet_scales Passed to [ggplot2::facet_grid()].
#' @param facet_space Passed to [ggplot2::facet_grid()]. If `NULL`, inferred from
#'   `facet_scales`.
#' @param auto_labels Logical; if `TRUE`, attempt to shorten and wrap x labels.
#' @param label_wrap_width Optional integer wrap width for x labels.
#' @param label_angle Optional angle for x-axis labels.
#' @param label_dodge Optional number of rows used to dodge x-axis labels.
#' @param label_check_overlap Logical; passed to [ggplot2::guide_axis()].
#'
#' @return A ggplot object.
#' @export
stacked_ribbon_plot <- function(data,
                                x,
                                y,
                                fill,
                                facet_rows = NULL,
                                facet_cols = NULL,
                                x_levels = NULL,
                                fill_levels = NULL,
                                bar_width = 0.5,
                                ribbon_alpha = 0.5,
                                ribbon_colour = "grey60",
                                bar_colour = "black",
                                bar_linewidth = 0.2,
                                ribbon_n = 200,
                                reverse_y = TRUE,
                                facet_scales = "fixed",
                                facet_space = NULL,
                                auto_labels = TRUE,
                                label_wrap_width = NULL,
                                label_angle = NULL,
                                label_dodge = NULL,
                                label_check_overlap = TRUE) {

  needed <- c("ggplot2", "rlang", "tidyr")
  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop(
      "Missing suggested packages for `stacked_ribbon_plot()`: ",
      paste(needed[!ok], collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(x) || length(x) != 1L) {
    stop("`x` must be a single character string naming a column in `data`.", call. = FALSE)
  }
  if (!is.character(y) || length(y) != 1L) {
    stop("`y` must be a single character string naming a column in `data`.", call. = FALSE)
  }
  if (!is.character(fill) || length(fill) != 1L) {
    stop("`fill` must be a single character string naming a column in `data`.", call. = FALSE)
  }

  facet_rows <- facet_rows %||% character(0)
  facet_cols <- facet_cols %||% character(0)

  if (!is.character(facet_rows)) {
    stop("`facet_rows` must be NULL or a character vector.", call. = FALSE)
  }
  if (!is.character(facet_cols)) {
    stop("`facet_cols` must be NULL or a character vector.", call. = FALSE)
  }

  needed_cols <- unique(c(x, y, fill, facet_rows, facet_cols))
  missing_cols <- setdiff(needed_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[y]])) {
    stop(sprintf("Column `%s` must be numeric.", y), call. = FALSE)
  }

  facet_vars <- c(facet_rows, facet_cols)

  if (is.null(x_levels)) {
    x_levels <- if (is.factor(data[[x]])) levels(data[[x]]) else unique(as.character(data[[x]]))
  } else {
    x_levels <- as.character(x_levels)
  }

  if (is.null(fill_levels)) {
    fill_levels <- if (is.factor(data[[fill]])) levels(data[[fill]]) else unique(as.character(data[[fill]]))
  } else {
    fill_levels <- as.character(fill_levels)
  }

  if (is.null(facet_space)) {
    facet_space <- if (facet_scales %in% c("free_x", "free", "free_y")) facet_scales else "fixed"
  }

  make_ribbon <- function(ymin1, ymax1, ymin2, ymax2, x1, x2, n = 200) {
    xseq <- seq(x1, x2, length.out = n)
    t <- (xseq - x1) / (x2 - x1)
    s <- 3 * t^2 - 2 * t^3

    y_top <- ymax1 + (ymax2 - ymax1) * s
    y_bot <- ymin1 + (ymin2 - ymin1) * s

    data.frame(
      x = c(xseq, rev(xseq)),
      y = c(y_top, rev(y_bot))
    )
  }

  escape_regex <- function(z) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", z)
  }

  wrap_one <- function(z, width) {
    paste(base::strwrap(z, width = width), collapse = "\n")
  }

  shorten_labels <- function(lbls, data, facet_vars) {
    out <- as.character(lbls)

    if (length(facet_vars) > 0) {
      facet_vals <- unique(unlist(lapply(facet_vars, function(v) as.character(unique(data[[v]])))))
      facet_vals <- facet_vals[!is.na(facet_vals) & nzchar(facet_vals)]

      for (tok in facet_vals[order(nchar(facet_vals), decreasing = TRUE)]) {
        tok_esc <- escape_regex(tok)
        out <- gsub(paste0("(^|[_ .-])", tok_esc, "($|[_ .-])"), "\\1\\2", out, perl = TRUE)
      }

      out <- gsub("^[_ .-]+|[_ .-]+$", "", out, perl = TRUE)
      out <- gsub("[_ .-]{2,}", "_", out, perl = TRUE)
    }

    out
  }

  df0 <- dplyr::ungroup(data)
  df0 <- dplyr::select(df0, dplyr::all_of(c(facet_vars, x, fill, y)))

  df0$x_cat <- factor(df0[[x]], levels = x_levels)
  df0$fill_cat <- factor(df0[[fill]], levels = fill_levels)
  df0$y_val <- df0[[y]]

  df0 <- dplyr::select(df0, dplyr::all_of(facet_vars), x_cat, fill_cat, y_val)

  df1 <- df0 |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, "x_cat", "fill_cat")))) |>
    dplyr::summarise(y_val = sum(y_val, na.rm = TRUE), .groups = "drop")

  if (length(facet_vars) > 0) {
    present_x <- df1 |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, "x_cat")))) |>
      dplyr::summarise(.groups = "drop")

    template <- tidyr::crossing(
      present_x,
      fill_cat = factor(fill_levels, levels = fill_levels)
    )

    join_by <- c(facet_vars, "x_cat", "fill_cat")
  } else {
    present_x <- df1 |>
      dplyr::distinct(x_cat)

    template <- tidyr::crossing(
      present_x,
      fill_cat = factor(fill_levels, levels = fill_levels)
    )

    join_by <- c("x_cat", "fill_cat")
  }

  df2 <- dplyr::left_join(template, df1, by = join_by)
  df2$y_val[is.na(df2$y_val)] <- 0

  df2 <- df2 |>
    dplyr::mutate(
      x_id = match(as.character(x_cat), x_levels)
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(facet_vars, "x_id"))), fill_cat) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, "x_id")))) |>
    dplyr::mutate(
      ymax = cumsum(y_val),
      ymin = dplyr::lag(ymax, default = 0)
    ) |>
    dplyr::ungroup()

  segs <- df2 |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(facet_vars, "fill_cat", "x_id")))) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, "fill_cat")))) |>
    dplyr::mutate(
      x_next = dplyr::lead(x_id),
      ymin_next = dplyr::lead(ymin),
      ymax_next = dplyr::lead(ymax)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(x_next))

  ribbon_list <- lapply(seq_len(nrow(segs)), function(i) {
    row_i <- segs[i, , drop = FALSE]

    rib <- make_ribbon(
      ymin1 = row_i$ymin,
      ymax1 = row_i$ymax,
      ymin2 = row_i$ymin_next,
      ymax2 = row_i$ymax_next,
      x1 = row_i$x_id + bar_width / 2,
      x2 = row_i$x_next - bar_width / 2,
      n = ribbon_n
    )

    rib$fill_cat <- row_i$fill_cat[1]
    rib$seg_id <- i

    if (length(facet_vars) > 0) {
      rib[facet_vars] <- row_i[rep(1, nrow(rib)), facet_vars, drop = FALSE]
    }

    rib
  })

  ribbons <- dplyr::bind_rows(ribbon_list)

  axis_labels <- x_levels

  if (isTRUE(auto_labels)) {
    axis_labels <- shorten_labels(axis_labels, data, facet_vars)

    if (is.null(label_wrap_width)) {
      max_chars <- max(nchar(axis_labels), na.rm = TRUE)
      label_wrap_width <- if (max_chars > 18) 12 else if (max_chars > 12) 16 else NA_integer_
    }

    if (!is.na(label_wrap_width)) {
      axis_labels <- vapply(axis_labels, wrap_one, character(1), width = label_wrap_width)
    }

    max_label_lines <- max(vapply(strsplit(axis_labels, "\n", fixed = TRUE), length, integer(1)))
    max_label_chars <- max(nchar(gsub("\n", "", axis_labels, fixed = TRUE)), na.rm = TRUE)

    if (is.null(label_angle)) {
      label_angle <- if (max_label_lines > 1) 0 else if (max_label_chars > 16) 45 else 0
    }

    if (is.null(label_dodge)) {
      label_dodge <- if (max_label_lines > 1) 1 else if (max_label_chars > 10) 2 else 1
    }
  } else {
    if (is.null(label_angle)) label_angle <- 0
    if (is.null(label_dodge)) label_dodge <- 1
  }

  x_guide <- ggplot2::guide_axis(
    angle = label_angle,
    n.dodge = label_dodge,
    check.overlap = label_check_overlap
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbons,
      ggplot2::aes(x = x, y = y, group = seg_id, fill = fill_cat),
      alpha = ribbon_alpha,
      colour = ribbon_colour
    ) +
    ggplot2::geom_rect(
      data = df2,
      ggplot2::aes(
        xmin = x_id - bar_width / 2,
        xmax = x_id + bar_width / 2,
        ymin = ymin,
        ymax = ymax,
        fill = fill_cat
      ),
      colour = bar_colour,
      linewidth = bar_linewidth
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = axis_labels,
      guide = x_guide
    ) +
    ggplot2::labs(x = NULL, y = y, fill = fill) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        hjust = if (isTRUE(label_angle > 0)) 1 else 0.5,
        vjust = if (isTRUE(label_angle > 0)) 1 else 0.5
      )
    )

  if (reverse_y) {
    p <- p + ggplot2::scale_y_reverse()
  } else {
    p <- p + ggplot2::scale_y_continuous() + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
  }

  if (length(facet_rows) > 0 || length(facet_cols) > 0) {
    row_vars <- if (length(facet_rows) > 0) ggplot2::vars(!!!rlang::syms(facet_rows)) else ggplot2::vars()
    col_vars <- if (length(facet_cols) > 0) ggplot2::vars(!!!rlang::syms(facet_cols)) else ggplot2::vars()

    p <- p + ggplot2::facet_grid(
      rows = row_vars,
      cols = col_vars,
      scales = facet_scales,
      space = facet_space
    )
  }
  p
}

# internal helper for NULL defaulting
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

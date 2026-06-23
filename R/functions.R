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

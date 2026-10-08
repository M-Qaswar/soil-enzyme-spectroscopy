# RRSEL_pipeline.R
# Implements Sections 2.3-2.9 of the revised manuscript.
# R 4.5.1
#
# Core principles:
# 1. External validation samples are never used for model selection.
# 2. Cross-sample transformations (centering, scaling, PCA) are calibration-only.
# 3. PLSR LV selection uses calibration-only LOOCV.
# 4. RRSEL uses calibration-only out-of-fold predictions.
# 5. Validation uncertainty uses paired bootstrap resampling.
#
# IMPORTANT:
# The manuscript currently contains an unresolved sample-count inconsistency:
# BG is reported as n=111, whereas Section 2.3 says all spectra were obtained
# for 92 samples. This script requires 111 BG observations if the reported
# 78/33 partition is to be reproduced; ACP and ARS require 92.

SEED_PARTITION <- 125L
SEED_RESAMPLE  <- 126L
N_BOOT_VIP <- 500L
N_BOOT_CI  <- 2000L
MAX_LV <- 20L

PARTITION <- list(
  BG  = list(frac = 0.70, pc = 11L, n_cal = 78L, n_val = 33L),
  ACP = list(frac = 0.73, pc = 10L, n_cal = 67L, n_val = 25L),
  ARS = list(frac = 0.75, pc = 11L, n_cal = 69L, n_val = 23L)
)

MIR_RANGE <- list(
  BG  = c(3954.7, 853.6),
  ACP = c(3928.6, 853.6),
  ARS = c(3954.7, 853.6)
)

# Table 2: ma = moving average width; sg = c(window, polynomial, derivative);
# gd = c(window, polynomial, derivative); fs = c(window, polynomial, derivative).
TABLE2 <- list(
  BG = list(
    vis = list(ma=3, sg=c(3,2,1), gd=c(1,7,3), fs=c(5,2,0)),
    mir = list(ma=7, sg=c(3,2,1), gd=c(1,7,3), fs=c(5,2,0))
  ),
  ACP = list(
    vis = list(ma=5, sg=c(9,2,0), gd=NULL, fs=NULL),
    mir = list(ma=5, sg=c(3,2,1), gd=c(1,9,3), fs=c(5,2,0))
  ),
  ARS = list(
    vis = list(ma=5, sg=c(3,2,1), gd=c(1,7,3), fs=c(3,2,0)),
    mir = list(ma=9, sg=c(3,2,1), gd=c(1,7,3), fs=c(5,2,0))
  )
)

required_pkgs <- c("prospectr", "caret", "pls", "glmnet")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  stop("Install required R packages first: ", paste(missing, collapse=", "))
}

to_numeric_axis <- function(x) {
  z <- suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))
  if (all(is.na(z))) stop("Spectral column names must contain numeric wavelength/wavenumber values.")
  z
}

spectral_matrix <- function(df, lo, hi, descending = FALSE) {
  nm <- names(df)
  ax <- suppressWarnings(as.numeric(nm))
  ok <- !is.na(ax) & ax >= min(lo, hi) & ax <= max(lo, hi)
  if (!any(ok)) stop("No spectral columns found in requested range ", lo, "-", hi)
  X <- as.matrix(df[, ok, drop=FALSE])
  storage.mode(X) <- "double"
  a <- ax[ok]
  ord <- if (descending) order(a, decreasing=TRUE) else order(a)
  X[, ord, drop=FALSE]
}

reflectance_to_absorbance <- function(X) {
  # Defensive numerical guard only. The manuscript states affected channels
  # are outside the retained 394-1674 nm modelling range.
  R <- X
  R[R <= 0] <- 1e-6
  -log10(R)
}

junction_correct <- function(X, axis) {
  # Detector transition 1024-1044 nm. Correct the shorter-wave segment
  # to the adjacent detector response, then linearly interpolate the gap.
  a <- as.numeric(axis)
  if (max(a) < 1024 || min(a) > 1044) return(X)
  out <- X
  left <- which(a < 1024)
  right <- which(a > 1044)
  gap <- which(a >= 1024 & a <= 1044)
  if (!length(left) || !length(right) || !length(gap)) return(out)
  for (i in seq_len(nrow(out))) {
    lv <- out[i, max(left)]
    rv <- out[i, min(right)]
    # Offset the short-wave segment to the adjacent detector level.
    # A linear bridge is then imposed through the junction.
    offset <- rv - lv
    out[i, left] <- out[i, left] + offset
    out[i, gap] <- approx(
      x=c(a[max(left)], a[min(right)]),
      y=c(out[i,max(left)], out[i,min(right)]),
      xout=a[gap]
    )$y
  }
  out
}

moving_average <- function(X, k) {
  if (is.null(k) || k <= 1) return(X)
  if (k %% 2 == 0) stop("Moving-average width must be odd.")
  w <- rep(1/k, k)
  t(apply(X, 1, function(z) stats::filter(z, w, sides=2)))
}

sg_filter <- function(X, w, p, d=0) {
  if (is.null(w)) return(X)
  if (w %% 2 == 0 || w <= p) stop("Savitzky-Golay window must be odd and > polynomial order.")
  # prospectr::savitzkyGolay returns rows x columns for matrix input.
  prospectr::savitzkyGolay(X, m=p, p=w, w=d)
}

# The prospectr API has varied across versions; this wrapper keeps the
# analytical definition explicit and stops rather than silently changing it.
apply_pipeline <- function(X, ma, sg, gd=NULL, fs=NULL) {
  Y <- moving_average(X, ma)
  if (anyNA(Y)) {
    # Edge values of centered moving average are restored by nearest available
    # values to avoid introducing sample-wise information.
    for (i in seq_len(nrow(Y))) {
      bad <- is.na(Y[i,])
      if (any(bad)) Y[i,bad] <- X[i,bad]
    }
  }
  # Row-wise vector normalization.
  nr <- sqrt(rowSums(Y^2))
  nr[nr == 0] <- 1
  Y <- Y / nr

  # prospectr convention: m = derivative order, p = polynomial order,
  # w = window size.
  Y <- prospectr::savitzkyGolay(Y, m=sg[3], p=sg[2], w=sg[1])

  if (!is.null(gd)) {
    Y <- prospectr::savitzkyGolay(Y, m=gd[3], p=gd[2], w=gd[1])
  }
  if (!is.null(fs)) {
    Y <- prospectr::savitzkyGolay(Y, m=fs[3], p=fs[2], w=fs[1])
  }
  Y
}

fit_scaler <- function(X) {
  mu <- colMeans(X)
  sdv <- apply(X, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(mean=mu, sd=sdv)
}
transform_scaler <- function(X, s) sweep(sweep(X, 2, s$mean, "-"), 2, s$sd, "/")

fit_pls_loocv <- function(X, y, max_lv=20L, seed=SEED_RESAMPLE) {
  nlv <- min(max_lv, nrow(X)-1L, ncol(X))
  if (nlv < 1) stop("Not enough calibration observations/predictors for PLSR.")
  ctrl <- caret::trainControl(method="LOOCV", savePredictions="final", allowParallel=FALSE)
  set.seed(seed)
  fit <- caret::train(
    x=X, y=y,
    method="pls",
    preProcess=c("center","scale"),
    tuneGrid=data.frame(ncomp=seq_len(nlv)),
    trControl=ctrl,
    metric="RMSE"
  )
  fit
}

predict_pls <- function(fit, X) as.numeric(predict(fit, newdata=X))

# Kennard-Stone on retained VIS-NIR PC scores. The manuscript specifies
# Mahalanobis distance calculated from the leading PCs.
kenstone <- function(X, k, pc) {
  if (k <= 0 || k >= nrow(X)) stop("Invalid calibration size.")
  pca <- prcomp(X, center=TRUE, scale.=FALSE, rank.=pc)
  Z <- pca$x[, seq_len(pc), drop=FALSE]
  # Mahalanobis distance in the retained PC space.
  covZ <- stats::cov(Z)
  invcov <- solve(covZ + diag(1e-10, nrow(covZ)))
  D <- as.matrix(dist(Z))
  # Initial pair: farthest Euclidean pair is equivalent to the usual
  # Mahalanobis construction after whitening only if covariance is identity;
  # therefore compute Mahalanobis distances explicitly.
  M <- matrix(0, nrow(Z), nrow(Z))
  for (i in seq_len(nrow(Z))) {
    d <- sweep(Z, 2, Z[i,], "-")
    M[i,] <- sqrt(rowSums((d %*% invcov) * d))
  }
  i <- which(M == max(M), arr.ind=TRUE)[1,]
  selected <- unique(i)
  while (length(selected) < k) {
    remain <- setdiff(seq_len(nrow(Z)), selected)
    mind <- sapply(remain, function(j) min(M[j, selected]))
    selected <- c(selected, remain[which.max(mind)])
  }
  sort(selected)
}

pca_scores <- function(Xcal, Xval, threshold=.95) {
  fit <- prcomp(Xcal, center=TRUE, scale.=TRUE)
  prop <- cumsum(fit$sdev^2) / sum(fit$sdev^2)
  q <- which(prop >= threshold)[1]
  list(
    cal=fit$x[,seq_len(q),drop=FALSE],
    val=scale(Xval, center=fit$center, scale=fit$scale) %*%
      fit$rotation[,seq_len(q),drop=FALSE],
    fit=fit, ncomp=q
  )
}

oof_predictions <- function(X, y, lv, seed=SEED_RESAMPLE) {
  n <- nrow(X)
  pred <- numeric(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    fit <- caret::train(
      x=X[tr,,drop=FALSE], y=y[tr],
      method="pls",
      preProcess=c("center","scale"),
      tuneGrid=data.frame(ncomp=lv),
      trControl=caret::trainControl(method="none", allowParallel=FALSE)
    )
    pred[i] <- predict(fit, X[i,,drop=FALSE])
  }
  pred
}

fit_ridge <- function(P, y, seed=SEED_RESAMPLE) {
  set.seed(seed)
  cv.glmnet <- glmnet::cv.glmnet(
    x=as.matrix(P), y=y, alpha=0, nfolds=min(10, nrow(P)),
    standardize=TRUE, type.measure="mse"
  )
  cv.glmnet
}

vip_from_pls <- function(fit) {
  # Extract predictor matrix and scores/loadings from pls object stored by caret.
  obj <- fit$finalModel
  W <- obj$loading.weights
  T <- obj$scores
  Q <- obj$Yloadings
  h <- min(ncol(W), length(Q))
  W <- W[,seq_len(h),drop=FALSE]
  T <- T[,seq_len(h),drop=FALSE]
  Q <- as.numeric(Q[seq_len(h),])
  s <- colSums((T^2) %*% t(Q^2))
  vip <- sqrt(nrow(W) * rowSums(sweep(W^2, 2, s / sum(s), "*")))
  as.numeric(vip)
}

bootstrap_vip <- function(X, y, lv, seed=SEED_RESAMPLE, B=N_BOOT_VIP) {
  set.seed(seed)
  V <- matrix(NA_real_, nrow=B, ncol=ncol(X))
  for (b in seq_len(B)) {
    ii <- sample(seq_len(nrow(X)), replace=TRUE)
    fit <- caret::train(
      x=X[ii,,drop=FALSE], y=y[ii],
      method="pls",
      preProcess=c("center","scale"),
      tuneGrid=data.frame(ncomp=lv),
      trControl=caret::trainControl(method="none", allowParallel=FALSE)
    )
    V[b,] <- vip_from_pls(fit)
  }
  list(
    mean=colMeans(V, na.rm=TRUE),
    sd=apply(V,2,sd,na.rm=TRUE),
    lo=apply(V,2,quantile,probs=.025,na.rm=TRUE),
    hi=apply(V,2,quantile,probs=.975,na.rm=TRUE)
  )
}

stable_regions <- function(lo, axis, threshold=1, min_run=5, merge_gap=15) {
  ok <- lo > threshold
  runs <- list()
  i <- 1L
  while (i <= length(ok)) {
    if (!ok[i]) { i <- i+1L; next }
    j <- i
    while (j < length(ok) && ok[j+1L]) j <- j+1L
    if (j-i+1L >= min_run) runs[[length(runs)+1L]] <- c(i,j)
    i <- j+1L
  }
  if (length(runs) <= 1) return(runs)
  merged <- list(runs[[1]])
  for (r in runs[-1]) {
    if (r[1]-merged[[length(merged)]][2]-1 <= merge_gap) {
      merged[[length(merged)]][2] <- r[2]
    } else merged[[length(merged) + 1L]] <- r
  }
  lapply(merged, function(r) c(axis[r[1]], axis[r[2]]))
}

metrics <- function(obs, pred) {
  e <- pred-obs
  rmse <- sqrt(mean(e^2))
  mae <- mean(abs(e))
  r2 <- if (sd(obs)==0) NA_real_ else cor(obs,pred)^2
  rpd <- sd(obs)/rmse
  rpiq <- IQR(obs)/rmse
  c(R2=r2, RMSE=rmse, MAE=mae, RPD=rpd, RPIQ=rpiq)
}

bootstrap_ci <- function(obs, pred, B=N_BOOT_CI, seed=SEED_RESAMPLE) {
  set.seed(seed)
  n <- length(obs)
  out <- matrix(NA_real_, B, 5)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace=TRUE)
    out[b,] <- metrics(obs[ii], pred[ii])
  }
  lo <- apply(out,2,quantile,.025,na.rm=TRUE)
  hi <- apply(out,2,quantile,.975,na.rm=TRUE)
  list(point=metrics(obs,pred), lo=lo, hi=hi, draws=out)
}

compare_models <- function(obs, pred_list) {
  E <- sapply(pred_list, function(p) abs(p-obs))
  colnames(E) <- names(pred_list)
  fr <- stats::friedman.test(E)
  pw <- matrix(NA_real_, ncol=ncol(E), nrow=ncol(E),
               dimnames=list(colnames(E),colnames(E)))
  pv <- c()
  pairs <- list()
  z <- 1L
  for (i in seq_len(ncol(E)-1L)) for (j in (i+1L):ncol(E)) {
    p <- stats::wilcox.test(E[,i],E[,j],paired=TRUE,exact=FALSE)$p.value
    pv[z] <- p
    pairs[[z]] <- c(colnames(E)[i],colnames(E)[j]); z <- z+1L
  }
  adj <- p.adjust(pv, method="holm")
  for (k in seq_along(pv)) pw[pairs[[k]][1],pairs[[k]][2]] <- adj[k]
  # Kendall's W = Friedman statistic / (n*(k-1))
  n <- nrow(E); k <- ncol(E)
  W <- as.numeric((fr$statistic) / (n*(k-1)))
  list(friedman=fr, wilcoxon_holm=pw, kendall_W=W, errors=E)
}

bland_altman <- function(obs, pred) {
  bias <- mean(pred-obs)
  sdv <- sd(pred-obs)
  data.frame(
    mean=(obs+pred)/2,
    difference=pred-obs,
    bias=bias,
    lower=bias-1.96*sdv,
    upper=bias+1.96*sdv
  )
}

prepare_blocks <- function(vis, mir, enzyme) {
  if (!enzyme %in% names(TABLE2)) stop("Unknown enzyme: ", enzyme)
  if (!all(c("sampleID","BG","ACP","ALP","ARS") %in% names(vis)))
    stop("VIS file is missing required metadata columns.")
  if (!all(c("sampleID") %in% names(mir)))
    stop("MIR file must contain sampleID.")
  if (enzyme == "BG" && nrow(vis) != 111)
    stop("BG requires 111 observations to reproduce the reported 78/33 partition.")
  if (enzyme != "BG" && nrow(vis) != 92)
    stop(enzyme, " requires 92 observations.")
  if (nrow(mir) != nrow(vis))
    stop("VIS and MIR must have identical row counts for this enzyme.")
  if (!identical(as.character(vis$sampleID), as.character(mir$sampleID)))
    stop("VIS and MIR sampleID order must be identical.")

  vx <- spectral_matrix(vis, c(394,1674), descending=FALSE)
  va <- as.numeric(names(vis)[!is.na(suppressWarnings(as.numeric(names(vis)))) &
                               suppressWarnings(as.numeric(names(vis))) >= 394 &
                               suppressWarnings(as.numeric(names(vis))) <= 1674])
  # Recover axis in sorted order.
  vnames <- names(vis)
  vax <- suppressWarnings(as.numeric(vnames))
  vok <- !is.na(vax) & vax>=394 & vax<=1674
  vax <- vax[vok]; ord <- order(vax); vax <- vax[ord]
  vx <- vx[,order(seq_along(vax)),drop=FALSE]

  # Re-extract in deterministic axis order.
  vx <- as.matrix(vis[, which(vok)[ord], drop=FALSE])
  storage.mode(vx) <- "double"
  vx <- reflectance_to_absorbance(vx)
  vx <- junction_correct(vx, vax)

  mr <- MIR_RANGE[[enzyme]]
  maxmin <- sort(mr)
  maxmin2 <- spectral_matrix(mir, mr, descending=TRUE)
  maxmin_names <- names(mir)
  maxmin_axis <- suppressWarnings(as.numeric(maxmin_names))
  mok <- !is.na(maxmin_axis) & maxmin_axis>=min(mr) & maxmin_axis<=max(mr)
  maxmin_axis <- maxmin_axis[mok]
  maxmin_axis <- sort(maxmin_axis, decreasing=TRUE)
  mx <- as.matrix(mir[, which(mok)[order(suppressWarnings(as.numeric(maxmin_names[mok])), decreasing=TRUE)], drop=FALSE])
  storage.mode(mx) <- "double"

  list(
    y=as.numeric(vis[[enzyme]]),
    vis=apply_pipeline(vx, TABLE2[[enzyme]]$vis$ma,
                       TABLE2[[enzyme]]$vis$sg,
                       TABLE2[[enzyme]]$vis$gd,
                       TABLE2[[enzyme]]$vis$fs),
    mir=apply_pipeline(mx, TABLE2[[enzyme]]$mir$ma,
                       TABLE2[[enzyme]]$mir$sg,
                       TABLE2[[enzyme]]$mir$gd,
                       TABLE2[[enzyme]]$mir$fs),
    vis_axis=vax,
    mir_axis=maxmin_axis
  )
}

run_enzyme <- function(vis, mir, enzyme) {
  set.seed(SEED_PARTITION)
  dat <- prepare_blocks(vis, mir, enzyme)
  pp <- PARTITION[[enzyme]]
  cal <- kenstone(dat$vis, pp$n_cal, pp$pc)
  val <- setdiff(seq_len(nrow(dat$vis)), cal)

  Xv <- dat$vis
  Xm <- dat$mir
  y <- dat$y

  fit1 <- fit_pls_loocv(Xv[cal,,drop=FALSE], y[cal])
  fit2 <- fit_pls_loocv(Xm[cal,,drop=FALSE], y[cal])
  Xf <- cbind(Xv, Xm)
  fit3 <- fit_pls_loocv(Xf[cal,,drop=FALSE], y[cal])

  # M4: calibration-only PCA for each block, then PLSR on concatenated PC scores.
  p1 <- pca_scores(Xv[cal,,drop=FALSE], Xv[val,,drop=FALSE])
  p2 <- pca_scores(Xm[cal,,drop=FALSE], Xm[val,,drop=FALSE])
  Pcal <- cbind(p1$cal,p2$cal)
  Pval <- cbind(p1$val,p2$val)
  fit4 <- fit_pls_loocv(Pcal, y[cal])

  # M5: base LV values selected from the calibration set.
  lv1 <- fit1$bestTune$ncomp
  lv2 <- fit2$bestTune$ncomp
  lv3 <- fit3$bestTune$ncomp
  o1 <- oof_predictions(Xv[cal,,drop=FALSE], y[cal], lv1)
  o2 <- oof_predictions(Xm[cal,,drop=FALSE], y[cal], lv2)
  o3 <- oof_predictions(Xf[cal,,drop=FALSE], y[cal], lv3)
  meta <- fit_ridge(cbind(o1,o2,o3), y[cal])

  pred1 <- predict_pls(fit1, Xv[val,,drop=FALSE])
  pred2 <- predict_pls(fit2, Xm[val,,drop=FALSE])
  pred3 <- predict_pls(fit3, Xf[val,,drop=FALSE])
  pred4 <- predict_pls(fit4, Pval)

  # Refit base learners are already full-calibration fits from fit1-3.
  base_val <- cbind(pred1,pred2,pred3)
  pred5 <- as.numeric(predict(meta, newx=base_val, s="lambda.min"))

  preds <- list(M1_VIS_NIR=pred1, M2_MIR=pred2, M3_Direct_Fusion=pred3,
                M4_Fusion_PCA=pred4, M5_RRSEL=pred5)
  perf <- do.call(rbind, lapply(preds, function(p) metrics(y[val],p)))
  cis <- lapply(preds, function(p) bootstrap_ci(y[val],p))
  vip1 <- bootstrap_vip(Xv[cal,,drop=FALSE], y[cal], lv1)
  vip2 <- bootstrap_vip(Xm[cal,,drop=FALSE], y[cal], lv2)
  vip3 <- bootstrap_vip(Xf[cal,,drop=FALSE], y[cal], lv3)

  list(
    enzyme=enzyme,
    calibration=cal,
    validation=val,
    latent_variables=c(M1=lv1,M2=lv2,M3=lv3,M4=fit4$bestTune$ncomp),
    pca_components=c(VIS=p1$ncomp,MIR=p2$ncomp),
    predictions=preds,
    observed=y[val],
    performance=perf,
    confidence_intervals=cis,
    vip=list(M1=vip1,M2=vip2,M3=vip3),
    stable_regions=list(
      M1=stable_regions(vip1$lo,dat$vis_axis),
      M2=stable_regions(vip2$lo,dat$mir_axis),
      M3=stable_regions(vip3$lo,c(dat$vis_axis,dat$mir_axis))
    ),
    comparison=compare_models(y[val],preds),
    bland_altman=lapply(preds, function(p) bland_altman(y[val],p)),
    meta=meta,
    preprocessing=TABLE2[[enzyme]]
  )
}

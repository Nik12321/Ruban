#' SAC: package for finding global extrema global optimization method based on the selective averaging coordinate  with restrictions.
#'
#' @name SAC
#' @docType package
#'
NULL

##--------------------------------------------------------------------##
##                      Global optimization method                    ##
##                              based on                              ##
##     the selective averaging coordinate  with restrictions          ##
##--------------------------------------------------------------------##

SAC <- function(type_sac = c("sacNormal", "sacExtended", "sacIterative"), sacControl_sac = sacControl(type_sac),
                x_sac,
                delta_sac,
                fitness_sac,
                lower_sac,
                upper_sac,
                n_sac = sacControl_sac$n,
                e_sac = sacControl_sac$e,
                M_sac = sacControl_sac$M,
                y_sac = sacControl_sac$y,
                q_sac = sacControl_sac$q,
                kernelType_sac = sacControl_sac$kernelType,
                r_sac = sacControl_sac$r,
                s_sac = sacControl_sac$s) {
  if (type_sac == "sacNormal" || type_sac == "sacExtended" || type_sac == "sacIterative")
    sacMode <- get(type_sac)
  else
    sacMode <- get("sacNormal")
  result <- sacMode(type = type_sac,
                 x = x_sac,
                 delta = delta_sac,
                 fitness = fitness_sac,
                 lower = lower_sac,
                 upper = upper_sac,
                 n = n_sac,
                 e = e_sac,
                 M = M_sac,
                 y = y_sac,
                 q = q_sac,
                 kernelType = kernelType_sac,
                 r = r_sac,
                 s = s_sac)
  return(result)
}

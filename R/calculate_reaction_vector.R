#' check_matrix()
#'
#' This function of the package 'ReDirection' checks whether the user-defined
#' stoichiometry matrix is suitable for further processing.
#'
#' 'ReDirection', is reaction-centric. This means that the computations and results
#' are solely based on the reactions. The number of unique reactions must exceed
#' the number of reactants (reactions >= reactants + 2).
#'
#' Additionally, 'ReDirection' must know if the reactions are represented as
#' rows (TRUE) or columns (FALSE).
#'
#' 'ReDirection', also checks for linear dependence (rows, columns) and removes
#' the same. This can perturb the indices of the matrix and users are advised to
#' check the input stoichiometry matrix a priori.
#'
#' The checked and modified matrix is returned and processed further.
#'
#' @param input_mat This is the user-defined stoichiometry matrix of a generic
#' metabolic network and is a mandatory argument.
#'
#' @param rxy A mandatory logical argument that indicates the orientation of
#' reactions as rows (TRUE) or columns (FALSE) of the user-defined stoichiometry
#' matrix.
#'
#' @return output_mat A checked and modified version of the stoichiometry matrix
#' which is returned for further computations.
#'
#' @return flag An indicator of the suitability of the user-defined input
#' stoichiometry matrix for further computations.

check_matrix <- function(input_mat,rxy)
{
        flag <- 0
        output_mat <- input_mat
        nco <- ncol(input_mat)
        nro <- nrow(input_mat)
        if(rxy == "TRUE"){
                for(iteration in 1:nrow(input_mat))
                {
                     output_mat[nro,] <- output_mat[nro,] * (-1)
                     output_mat <- output_mat[!duplicated(t(output_mat)),]
                     nro <- nro - 1
                     flag <- 1
                }
                         }
        else if(rxy == "FALSE"){
                for(iteration in 1:ncol(input_mat))
                {
                     output_mat[,nco] <- output_mat[,nco] * (-1)
                     output_mat <- output_mat[,!duplicated(t(output_mat))]
                     nco <- nco - 1
                     flag <- 1
                }
                               }
        else if(rxy == "NULL"){print("Please enter whether the reactions are present as rows (TRUE) or columns (FALSE) in your stoichiometric matrix")}

        if(flag > 0)
        {
        rm(input_mat)
        input_mat <- output_mat
        nco <- ncol(input_mat)
        nro <- nrow(input_mat)
        addcol <- colSums(input_mat)
        addrow <- rowSums(input_mat)
        if( (rxy == "TRUE") && (nro >= nco + 2) ){flag <- 1}
        else if( (rxy == "FALSE") && (nco >= nro + 2) ){flag <- 1}

        }
list_to_return <- list(output_mat,flag)
return(list_to_return)
}


#' reaction_vector()
#'
#' This function of the package 'ReDirection', combinatorially sums non-zero
#' and unique reaction vectors of the null space and the generated subspace(s).
#'
#' The resultant reaction vectors (RVs) comprise every unique reaction as a
#' sequence of identical components across  all RVs and are processed further.
#'
#' The RVs are evaluated for the presence of duplicated vectors and the same are
#' removed.
#'
#' If and only if the combined propensity of every reaction fulfills the output
#' criteria (non-zero, greater than unity) are the iterations stopped, otherwise
#' another commences.
#'
#'
#' @param xsol The null space spanning vectors of the checked and modified user-
#' defined stoichiometry matrix.
#'
#' @return mtdf Vectors are the resultant reaction vectors which are computed by
#' by combinatorially summing the non-zero and null space spanning vectors and
#' the subspace(s) generated.
#'
#'


reaction_vector <- function(xsol)
{
    rm(space,vec,list_of_reactions,dt_frame)
    space <- dim(nrow(xsol))
    vec <- array(dim=c(nrow(xsol),1))
    dt_frame <- array(dim=c(nrow(xsol),1))
    list_of_reactions <- array(dim=c(nrow(xsol),1))

for (i in 1:ncol(xsol))
           {
                 space <- cbind(space,xsol[,i])
                 list_of_reactions <- rbind(list_of_reactions,i)
           }

space_cols <- ncol(space)

for (numcols in 2:space_cols)
{
Z <- combinations(space_cols,numcols,list_of_reactions)
for (k in 1:nrow(Z))
     {
      vec[,1] <- matrix(rep(0,times=nrow(xsol)),nrow=nrow(xsol),ncol=1)
      name_of_vector <- ""
      for (m in 1:ncol(Z))
           {
              if(m>1)
                     {
                            counter1 <- nrow(space)
                            if(counter1 == nrow(space))
                                                                {
                                                                 vec[,1] <- vec[,1] + space[,m]
                                                                 name_of_vector <- paste(name_of_vector,Z[k,m],sep="")
                                                                }

                     }
               else
                      {
                          vec[,1] <- vec[,1] + space[,m]
                          name_of_vector <- paste(name_of_vector,Z[k,m],sep="")
                      }
             }

    vector_name <- paste("v",name_of_vector,sep="")
    if(all(vec[,1]>=0.0)){vector_name <- paste(vector_name,"0+",sep="")}else{vector_name <- paste(vector_name,"-",sep="")}
    dt <- data.frame(vector_name,t(vec[,1]))
    dt_frame <- cbind(data.frame(vec[,1],row.names=NULL),dt_frame)
    rm(name_of_vector,dt)
     }
}

dt_frame <- dt_frame[,-c(length(dt_frame))]
dt_frame <- cbind(xsol,dt_frame)
mtdf <- data.matrix(dt_frame)
return(mtdf)
}


#' calculate_reaction_vector()
#
#' This is the main function of the package 'ReDirection'.
#'
#' The input stoichiometry matrix of a generic biochemical network is checked
#' and modified after which the null space is computed.
#'
#' The non-zero null space spanning and the subspace-generated unique reaction
#' vectors (RVs) are recursively summed combinatorially and evaluated for
#' redundancy.
#'
#' The function compares the generated subspace vectors reaction-wise, i.e., as
#' sequences of identical terms across all RVs. These are then evaluated on the
#' basis of pre-defined criteria (descriptive statistics, linear maps, tests of
#' convergence, probability of occurrence, vector norms).
#'
#' The output for each reaction is the combined propensity (p1-norm) of the sub-
#' propensities and is strictly positive. This is the probable rate constant, is
#' used to infer the dominant direction of each reaction and annotate a reaction
#' as "Forward (f)" or "Reverse (b)" or "Equivalent (e)".
#
#' Although, there is no restriction on the number of reactions, the inherent
#' computational complexity (NP-hard) involved per iteration suggests that an
#' suitable value for the number of reactions is around 20.
#'
#' @param smat An input stoichiometry matrix of a generic biochemical network.
#' This is a mandatory argument.
#
#' @param rar A mandatory logical argument that indicates the orientation of the
#' reactions as rows (TRUE) or columns (FALSE) in the stoichiometry matrix.
#'
#' @return code A numerically encoded ('0', no success; '1', success) text
#' message to the user and indicates the outcome of utilizing "ReDirection".
#
#' @examples
#' mx <- matrix(c(1,0,0,0,0,1,0,-1,0,1,0,0,1,0,0,0,-1,1,1,0,0,0,0,0,0,-2,0,-1,
#' -1,1,0,0,0,0,1,0,-1,1,0,0,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,-1,0,1,0,0,-1,0,-1,0)
#' ,byrow=TRUE,nrow=9,ncol=7)
#' calculate_reaction_vector(mx,TRUE)
#'
#' @import MASS pracma stats
#' @importFrom gtools combinations
#'
#' @export


calculate_reaction_vector <- function(smat,rar)
  {

if( (missing(smat)) && (missing(rar))  ){smat <- "NULL"; rar  <- "NULL"}
else if( (missing(smat)) && (missing(rar))  ){smat <- "NULL"; rar  <- "NULL"}
else if(missing(smat)){smat <- "NULL"}
else if(missing(rar)){rar <- "NULL"}


XL <- check_matrix(smat,rar)
code <- XL[[2]]

if(XL[[2]] > 0)
{

print("Stoichiometric matrix checked .....[OK]")

if( nrow(XL[[1]]) > ncol(XL[[1]]) ){sol <- Null(XL[[1]])}
else if( nrow(XL[[1]]) < ncol(XL[[1]]) ){sol <- Null(t(XL[[1]]))}


print("Null space calculated .....[OK]")
print("Commencing iterations.....[OK]")

iter <- 1

repeat
{

if(iter > 1){rm(mat_dt_frame,modes_mat,reaction,propensity,direction)}
mat_dt_frame <- paste("mat_dt_frame_",iter,sep="")
suppressWarnings(mat_dt_frame <- round(reaction_vector(sol),6))
rm(sol)
modes_mat <-mat_dt_frame [,!duplicated(t(mat_dt_frame))]

reaction <- dim(c(nrow(modes_mat),1))
propensity <- dim(c(nrow(modes_mat),1))
index <- dim(c(nrow(modes_mat),1))
direction <- dim(c(nrow(modes_mat),1))
min_value_to_calc <- dim(ncol(modes_mat))
max_value_to_calc <- dim(ncol(modes_mat))
eq_value_to_calc <- dim(ncol(modes_mat))



    for(nrw in 1:nrow(modes_mat))
      {
       min_val <- max_val <- max_counter <- min_counter <- eq_counter <- trueq_counter <-  f0 <- f1 <- F0 <- F1 <-  delta_max <- delta_min <- 0
       max_value_to_calc <- min_value_to_calc <- eq_value_to_calc <- rep(1,times=ncol(modes_mat))
       corrected_max_counter <- corrected_min_counter <- corrected_eq_counter <- peak_value <- flat_value <- 0
       dD <- "";F2 <- 0
       min_val <- mean(modes_mat[nrw,]) - (2 * sd(modes_mat[nrw,]))
       max_val <- mean(modes_mat[nrw,]) + (2 * sd(modes_mat[nrw,]))
       for(nco in 1:ncol(modes_mat))
            {if( (modes_mat[nrw,nco] < min_val)  ){min_counter=min_counter+1;min_value_to_calc[min_counter] <- modes_mat[nrw,nco]}
            if( (modes_mat[nrw,nco] > max_val)  ){max_counter=max_counter+1;max_value_to_calc[max_counter] <- modes_mat[nrw,nco]}
            if(modes_mat[nrw,nco] != 0){if( (modes_mat[nrw,nco] > -1) && (modes_mat[nrw,nco] < 1) ){eq_counter=eq_counter+1;eq_value_to_calc[eq_counter] <- modes_mat[nrw,nco]}}}

            index <- rbind(index,paste("Reaction",nrw,sep="_"))
            f2 <- Norm(eq_value_to_calc,p = Inf);


            if( (max_counter > 0) )
                         {

                            if ( (f2 > 0) && (f2 < 1) ){corrected_max_counter = corrected_max_counter + max_counter}
                            delta_max = max_counter - corrected_max_counter
                            f0 <- Norm(max_value_to_calc,p = Inf)
                            if( (delta_max > 0) && (f0 != 0) ){
                                                                                         if( (f0 < -1 ) || (f0 > 1) ) {F0 <- abs(round(f0,6))}
                                                                                         else if( (f0 > -1) && (f1 < 1) ){F0 <- round(exp(1) ^ (log(abs(f0 ^ (-1)),base = exp(1))),6)}
                                                                                 }

                         }
            if(  (min_counter > 0) )
                        {
                            if ( (f2 > 0) && (f2 < 1) ){corrected_min_counter = corrected_min_counter + min_counter}
                            delta_min = min_counter - corrected_min_counter
                            f1 <- Norm(min_value_to_calc,p = Inf)
                            if( (delta_min > 0) && (f1 != 0) ){
                                                                                      if((f1 < -1)  || (f1 > 1)){F1 <- round(exp(1) ^ (log(abs(f1) ^ (-1),base = exp(1))),6)}
                                                                                      else if( (f1 > -1) && (f1 < 1) ){F1 <- abs(round(f1,6) )}
                                                                                }
                        }

           T_counts <- delta_max + delta_min; if(T_counts == 0){T_counts = T_counts + 1}
           prob_F <- (delta_max)/T_counts
           prob_R <- (delta_min)/T_counts
           prob_E <- (corrected_eq_counter)/T_counts

           if( (prob_F > prob_R) ){prob_R <- prob_E <- 0}
           else if( (prob_R > prob_F) ){prob_F <- prob_E <- 0}
           else if( (prob_F == prob_R) && ((prob_F > 0) && (prob_R > 0)) ){prob_F <- prob_R <- 0;prob_E <- 1}
           else{prob_F <- prob_R <- prob_E <- 0}

           reaction <- rbind(reaction,paste(prob_F,prob_R,prob_E,sep=","))
           norm_prob <- c((F0 * prob_F),(F1 * prob_R),prob_E)
           propensity  <- rbind(propensity,Norm(norm_prob,p=1))
           if( Norm(norm_prob,p=1) > 1){dD <- "Forward"} else if( (Norm(norm_prob,p=1) > 0) && (Norm(norm_prob,p=1) < 1)){dD <- "Reverse"} else if(Norm(norm_prob,p=1) == 1){dD <- "Equivalent"}
           direction <- rbind(direction,dD)
      }
     if( (Norm(propensity,p=Inf) > 1) && (Norm(propensity,p=-Inf) > 0)  ) {sol <- modes_mat;break}
     else {iter = iter + 1; sol <- modes_mat}
}

display <- dim(c(nrow(modes_mat),4))
display <- data.frame(cbind(index,reaction,propensity,direction))
names(display) <- c("Index", "RV (f,b,e)", "Propensity", "Direction")

print("Dominant direction search complete.....[OK]")
print(display,row.names=FALSE)

code <- 1
} else{print("Stoichiometric matrix checked..... [FAIL]")}


return(code)
}


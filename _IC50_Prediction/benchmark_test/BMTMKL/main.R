library(dplyr)
library(RhpcBLASctl)

args_def = function(args, default=NULL) {
  args_x = args
  if (toupper(args_x) %in% c("T", "TRUE")) args_x = T
  if (toupper(args_x) %in% c("F", "FALSE")) args_x = F
  args_x = ifelse(!is.na(args_x), args_x, default)
  return(args_x)
}

args = commandArgs(trailingOnly=T)
# nth = args_def(args[1], 0)
# ic50_th = args_def(args[2], 2)
# choice = args_def(args[3], 0)
# cpus = args_def(args[4], 12)

args = commandArgs(trailingOnly=T)
nth = args_def(args[1], 0)
dir_pred = args_def(args[2], "Results/ex")
cpus = args_def(args[3], 12)

blas_get_num_procs()
blas_set_num_threads(cpus)

# dir_ic50 = c("IC50_GDSC", "IC50_GDSC1", "IC50_GDSC2")[ic50_th+1]
# dir_test = c("Normal", "Cell_Blind")[choice+1]
# dir_pred = sprintf("Results/%s/%s", dir_ic50, dir_test)

file_train = sprintf("%s/ic50_train_%s.txt", dir_pred, nth)
file_test = sprintf("%s/ic50_test_%s.txt", dir_pred, nth)
file_pred = sprintf("%s/pred_test_%s.csv", dir_pred, nth)

# Input Data
IC50_Train = read.csv(file_train, sep="\t")
IC50_Test = read.csv(file_test, sep="\t")
IC50_Train = IC50_Train %>% mutate(Cell=as.character(Cell), Drug=as.character(Drug))
IC50_Test = IC50_Test %>% mutate(Cell=as.character(Cell), Drug=as.character(Drug))

file = "_data/View.RData"
load(file, verbose=T)

n_drug_tn = IC50_Train$Drug %>% unique
n_drug_tt = IC50_Test$Drug %>% unique
drugs = intersect(n_drug_tn, n_drug_tt)

cond1 = length(drugs)!=length(n_drug_tn)
cond2 = length(drugs)!=length(n_drug_tt)

if (cond1 | cond2) {
  sprintf("# The number of drugs is decreased [Train] : %s > %s", 
          length(drugs), length(n_drug_tn)) %>% print
  sprintf("# The number of drugs is decreased [Test] : %s > %s", 
          length(drugs), length(n_drug_tt)) %>% print
  
  IC50_Train = IC50_Train %>% subset(Drug %in% drugs)
  IC50_Test = IC50_Test %>% subset(Drug %in% drugs)
}

n_view = 14
cells = rownames(RNA_Array_View)
view = sprintf("View%s", 1:n_view)

n_cell = length(cells)
dim_view = c(n_cell, n_cell, n_view)

Cell_List = list(Prot_View, RNA_Array_View, RNA_Seq_View, Meth_View, CNV_View, MUT_View, 
                 EXP_C2_View, EXP_CP_View, CNV_C2_View, CNV_CP_View, MUT_C2_View, MUT_CP_View,
                 Expr_CNV_View, RNA_Seq_PA_View)

Cell_List = Cell_List %>% lapply(as.matrix)
Cell_View = array(unlist(Cell_List), dim=dim_view, dimnames=list(cells, cells, view))

# for (i in 1:length(Cell_List)) {
#   identical(Cell_List[[i]], Cell_View[, , i]) %>% print   # T
# }


source("bayesian_multitask_multiple_kernel_learning_train.R")
source("bayesian_multitask_multiple_kernel_learning_test.R")

#initalize the parameters of the algorithm
parameters <- list()

#set the hyperparameters of gamma prior used for sample weights
parameters$alpha_lambda <- 1
parameters$beta_lambda <- 1

#set the hyperparameters of gamma prior used for intermediate noise
parameters$alpha_upsilon <- 1
parameters$beta_upsilon <- 1

#set the hyperparameters of gamma prior used for bias
parameters$alpha_gamma <- 1
parameters$beta_gamma <- 1

#set the hyperparameters of gamma prior used for kernel weights
parameters$alpha_omega <- 1
parameters$beta_omega <- 1

#set the hyperparameters of gamma prior used for output noise
parameters$alpha_epsilon <- 1
parameters$beta_epsilon <- 1

### IMPORTANT ###
#For gamma priors, you can experiment with three different (alpha, beta) values
#(1, 1) => default priors
#(1e-10, 1e+10) => good for obtaining sparsity
#(1e-10, 1e-10) => good for small sample size problems (like in Nature Biotechnology paper)

#set the number of iterations
parameters$iteration <- 200

#determine whether you want to calculate and store the lower bound values
parameters$progress <- 0

#set the seed for random number generator used to initalize random variables
parameters$seed <- 1606

#set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
T <- length(drugs)
#set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
P <- dim(Cell_View)[3]

#initialize the kernels and outputs of each task for training
Ktrain <- vector("list", T)
ytrain <- vector("list", T)
cells_train = list()

for (t in 1:T) {
  IC50_Temp = IC50_Train %>% subset(Drug==drugs[t])
  cells_omics = rownames(RNA_Array_View)
  cells_ic50 = unique(IC50_Temp$Cell)
  cells = intersect(cells_omics, cells_ic50)
  cells_train[[t]] = cells
  
  idx = match(cells, IC50_Temp$Cell)
  IC50_Temp = IC50_Temp[idx, ]
  
  Ktrain[[t]] <- Cell_View[cells, cells, ]
  # should be an Ntra x Ntra x P matrix containing similarity values between training samples of task t
  ytrain[[t]] <- IC50_Temp[, "LN_IC50", drop=F] %>% as.matrix
  # should be an Ntra x 1 matrix containing target outputs of task t
}

#perform training
state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)

# #display the kernel weights
# print(state$be$mu[(T+1):(T+P)])

#initialize the kernels of each task for testing
Prediction = data.frame()
Ktest <- vector("list", T)

for (t in 1:T) {
  IC50_Temp = IC50_Test %>% subset(Drug==drugs[t])
  cells_omics = rownames(RNA_Array_View)
  cells_ic50 = unique(IC50_Temp$Cell)
  cells = intersect(cells_omics, cells_ic50)
  
  idx = match(cells, IC50_Temp$Cell)
  IC50_Temp = IC50_Temp[idx, ]
  Prediction = Prediction %>% rbind(IC50_Temp)
  
  Ktest[[t]] <- Cell_View[cells_train[[t]], cells, ]
  # should be an Ntra x Ntest x P matrix containing similarity values between training and test samples of task t
}

#perform prediction
prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)

# #display the predictions for each task
# for (t in 1:T) {
#   print(prediction$y[[t]]$mu)
# }

y_pred = c()
for (t in 1:T) y_pred = y_pred %>% c(prediction$y[[t]]$mu)

Prediction$Prediction = y_pred
write.csv(Prediction, file=file_pred, row.names=F)

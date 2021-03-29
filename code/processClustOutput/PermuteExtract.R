library(batchtools)

loadRegistry('DSGEs_update/SWpermutations2019', writeable = TRUE)
res = reduceResultsList(fun=function(x) x)
res = do.call(rbind, res) # list to matrix
train_mse = res[,1:7]
test_mse = res[,8:14]
par_ests = res[,15:ncol(res)]

save(train_mse, test_mse, par_ests, file="DSGEs_Update/perm_results.Rdata")



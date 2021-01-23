mut <- read.csv(file = "features/Mutation.csv")
mut <- t(mut)
dump_features(R.object = )
test <- read.csv(file = "dream_data/dnaseq.csv")
test$lab_id %>% factor %>% levels %>% length

mut <- mut[intersect(row.names(mut), row.names(auc)),]
auc <- auc[intersect(row.names(mut), row.names(auc)),]
dump_features(R.object = mut, "features/mut/alex_mut.RData")
dump_features(R.object = auc, )
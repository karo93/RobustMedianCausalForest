Test_simulation2 <- read_excel("~/Documents/Test_simulation2.xlsx",  col_types = c( "numeric", "numeric", "numeric", "numeric", "numeric",  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric" ), skip = 1)
library(readxl)
#Test with small Dataset and prints for visualization
set.seed(12345)
set.seed(12)
tree <- causalTree::causalTree(y~ ., data = simulation.2, treatment = simulation.2$treatment,
                               split.Rule = "user", cv.option = "matching", split.Honest = F, cv.Honest = T, split.Bucket = F, 
                               xval = 5, cp = 0, minsize = 1, propensity = 0.5, split.alpha = 1)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]

opfit <- prune(tree, opcp)
rpart.plot(opfit, main="USER")

set.seed(12345)
tree1 <- causalTree::causalTree(y~ ., data = simulation.2, treatment = simulation.2$treatment,
                               split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                               xval = 5, cp = 0, minsize = 1, propensity = 0.5, split.alpha = 1)

opcp1 <- tree1$cptable[,1][which.min(tree1$cptable[,4])]

opfit <- prune(tree1, opcp1)
rpart.plot(opfit, main="CT")

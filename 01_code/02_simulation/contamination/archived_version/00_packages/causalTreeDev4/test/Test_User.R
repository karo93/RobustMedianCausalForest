set.seed(1)
tree <- causalTree::causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                   split.Rule = "user", cv.option = "matching", split.Honest = T, cv.Honest = F, split.Bucket = F, 
                   xval = 10)

tree$cptable
#rpart.plot(tree)
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
optree <- prune(tree, cp = opcp)
opfit <- prune(tree, opcp)
rpart.plot(optree, main="User_user")


tree1 <- causalTree::causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                               split.Rule = "user", cv.option = "CT", split.Honest = T, cv.Honest = F, split.Bucket = F, 
                               xval = 10)
tree1$cptable
rpart.plot(tree1)

opcp1<- tree1$cptable[,1][which.min(tree1$cptable[,4])]

opfit1 <- prune(tree1, opcp1)
rpart.plot(opfit1,main="user_CT")

tree2 <- causalTree::causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                                split.Rule = "CT", cv.option = "matching", split.Honest = F, cv.Honest = T, split.Bucket = F, 
                                xval = 10)

opcp2<- tree2$cptable[,1][which.min(tree2$cptable[,4])]

opfit2 <- prune(tree2, opcp2)
rpart.plot(opfit2,main="CT_matching")

tree3 <- causalTree::causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                                split.Rule = "user", cv.option = "fit", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                                xval = 5, cp = 0, minsize = 20, propensity = 0.5, split.alpha = 1)

opcp3<- tree3$cptable[,1][which.min(tree3$cptable[,4])]

opfit3 <- prune(tree3, opcp3)
rpart.plot(opfit3,main="User_fit")

tree4 <- causalTree::causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                                split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                                xval = 5, cp = 0, minsize = 20, propensity = 0.5, split.alpha = 1)

opcp4<- tree4$cptable[,1][which.min(tree4$cptable[,4])]

opfit4 <- prune(tree4, opcp4)
rpart.plot(opfit4,main="User_ToT")

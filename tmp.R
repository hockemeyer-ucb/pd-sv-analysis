library(ape)

myTree <- ape::read.tree(text = "((A, B), ((C, D), (E, F)));")

plot(myTree)

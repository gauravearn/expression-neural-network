deepNeuralExpression <- function(expression_file, epochs) {
    library(deepnet)
    # a deep neural expression based classifier demonstrated
    # to fit a unbalanced dataset, expression datasets across the 
    # samples where summed and then divided by the length of the total
    # replicates and then i classify them as those which have expression 
    # greater than 100 as 1 and those which have lower than 100 as 0 and build 
    # this neural network. Another optimization would be to the dropoutlayer 
    # and also present a normalizing approach, which can make the expression 
    # unbiasness looks straight. I am working on the same so that the deviations
    # caused by the change in the wide expression is not a problem. 
    expressionfile <- expression_file
    epoch = epochs
    expression <- read.csv(expressionfile,
        stringsAsFactors = FALSE, sep = "\t"
    )
    filtered_expression_complete <- expression[which(complete.cases(expression == TRUE)), ]
    expression_subset <- subset(filtered_expression_complete, select = -c(Gene.ID, Gene.Name))
    e.class <- vector(length = length(unlist(expression_subset["X7.day"], use.names = FALSE)))
    for (i in 1:length(e.class)) {
        e.class[i] <- sum(expression_subset[i, ]) / length(expression_subset[i, ])
    }
    final_expression <- cbind(expression_subset, e.class)
    final_round <- sapply(final_expression, as.integer)
    x <- as.numeric(as.matrix(final_round[, 1:10]))
    x <- matrix(as.numeric(x), ncol = 10)
    y <- final_round[, 11]
    y <- as.numeric(ifelse(y > 100, 1, 0))
    expression_sigm <- nn.train(x, y,
        hidden = c(14), activationfun = "sigm",
                     numepochs = as.integer(epoch), learningrate = 0.8
    )
    expression_sigm_predict <- nn.predict(expression_sigm, x)
    expression_sigm_y <- matrix(0, length(expression_sigm_predict), 1)
    expression_sigm_y[which(expression_sigm_predict > mean(expression_sigm_predict))] = 1
    expression_sigm_y[which(expression_sigm_predict <= mean(expression_sigm_predict))] = 0
    expression_sigm_table <- table(y, expression_sigm_y)
    expression_tanh <- nn.train(x, y,
        hidden = c(14), activationfun = "tanh",
                           numepochs = as.integer(epoch), learningrate = 0.8
    )
    expression_tanh_predict <- nn.predict(expression_tanh, x)
    expression_tanh_y <- matrix(0, length(expression_tanh_predict), 1)
    expression_tanh_y[which(expression_sigm_predict >
                                  mean(expression_tanh_predict))] = 1
    expression_tanh_y[which(expression_sigm_predict <=
                                    mean(expression_tanh_predict))] = 0
    expression_tanh_table <- table(y, expression_tanh_y)
    print(sum(diag(expression_sigm_table)) / sum(expression_sigm_table))
    print(sum(diag(expression_tanh_table)) / sum(expression_tanh_table))
}

### Workshop on Machine Learning in R
### author: Dmytro Fishman
### date: 29.09.2016

# Load corresponding library
library(caret)

# Load data into R
training <- read.csv('data/traindata.csv')

# Take a look at your data
head(training)
table(training$y)

# Enable reproducibility
set.seed(999)

# Parameters of cross validation 
fitControl <- trainControl(
  method = "cv",
  number = 4, # 4-fold CV, you can tune this number
  classProbs = TRUE
)

# Train model, 
# here you could change method to some other method
# from https://topepo.github.io/caret/modelList.html
training$y <- factor(ifelse(training$y==0, "Zero", "One"))
modelFit <- train(y~., 
                  data =  training,
                  method = "knn", 
                  trControl = fitControl
)

modelFit

# Build a confusion matrix
confusionMatrix(predict(modelFit, training),  training$y)

# Generate a submission
testing <- read.csv('data/testdata.csv')

# Generate predictions for testing data
predictions <- predict(modelFit, newdata = testing, type = 'prob')[,1]
submission <- data.frame(id = 1:length(predictions), y = predictions)

# Save into a file
write.csv(x = submission,
          file = 'sample_submission.csv', row.names = FALSE)